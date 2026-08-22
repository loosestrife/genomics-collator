from collections import defaultdict
from numpy import object_
import os
import json
import traceback
from typing import List, Optional
import duckdb
from fastapi import FastAPI, HTTPException, status
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field
import sqlite3
import ollama
import asyncio

app = FastAPI(title="NCBI Orthologs FastAPI Server")

NCBI_DB_PATH="../ncbi.duckdb"
ORTHODB_DB_PATH="../orthodb.duckdb"

from collections import defaultdict
import duckdb

from collections import defaultdict
import duckdb

def odb_mrca(
    conn: duckdb.DuckDBPyConnection, 
    species1: str, 
    species2: str, 
    sdata_lineage: list[tuple[int, str]] | list[str],
    tid2name: dict[int, str]
) -> tuple[int, str]:
    """
    Finds the most specific OrthoDB clade TaxID shared by two species by 
    matching against the ordered NCBI lineage tree.
    """
    # 1. Option A: Fast check using the `lineage` array in level2species
    l2s_array_query = """
        SELECT DISTINCT UNNEST(l1.lineage) AS level_tax_id
        FROM level2species l1
        JOIN species s1 ON l1.species_code = s1.species_code
        JOIN level2species l2 ON l1.species_code = s2.species_code
        JOIN species s2 ON s1.species_code != s2.species_code
        WHERE s1.species_name = $species1
          AND s2.species_name = $species2;
    """
    
    # 2. Option B: Check og2genes (Fixed query alias: og1.og_id instead of og.og_id)
    og_levels_query = """
        SELECT DISTINCT 
            CAST(regexp_extract(og1.og_id, 'at([0-9]+)$', 1) AS INTEGER) AS level_tax_id
        FROM species s1
        JOIN genes g1 ON s1.species_code = g1.species_code
        JOIN og2genes og1 ON g1.odb_gene_id = og1.protein_id
        JOIN og2genes og2 ON og1.og_id = og2.og_id
        JOIN genes g2 ON og2.protein_id = g2.odb_gene_id
        JOIN species s2 ON g2.species_code = s2.species_code
        WHERE s1.species_name = $species1
          AND s2.species_name = $species2;
    """

    shared_odb_levels = None
    # try:
    #     results = conn.execute(og_levels_query, {"species1": species1, "species2": species2}).fetchall()
    #     shared_odb_levels = {row[0] for row in results if row[0] is not None}
    # except Exception as e:
    #     print(f"og_levels_query fallback to level2species: {e}")
    #     shared_odb_levels = set()

    # Fallback to level2species if og2genes is empty
    if not shared_odb_levels:
        l2s_query = """
            SELECT DISTINCT UNNEST(list_intersect(l1.lineage, l2.lineage)) AS level_tax_id
FROM level2species l1
JOIN species s1 ON l1.species_code = s1.species_code
CROSS JOIN level2species l2
JOIN species s2 ON l2.species_code = s2.species_code
WHERE s1.species_name = $species1
  AND s2.species_name = $species2;
        """
        l2s_results = conn.execute(l2s_query, {"species1": species1, "species2": species2}).fetchall()
        print('lineage query results', l2s_results)
        shared_odb_levels = {row[0] for row in l2s_results}

    print('shared odb levels:', species1, species2, [(level, tid2name[level]) for level in shared_odb_levels])

    if not shared_odb_levels:
        raise ValueError(f"No common OrthoDB clade found for '{species1}' and '{species2}'.")

    # Walk NCBI lineage from most specific (index 0) to root
    for node in sdata_lineage:
        tax_id = node[0] if isinstance(node, (tuple, list)) else node
        if tax_id in shared_odb_levels:
            node_name = tid2name.get(
                tax_id, 
                str(node[1] if isinstance(node, (tuple, list)) else tax_id)
            )
            return (tax_id, node_name)

def odb_mrcall(
    conn: duckdb.DuckDBPyConnection,
    species_list: list[str],
    sdata_lineage: list[tuple[int, str]] | list[str],
    tid2name: dict[int, str]
) -> tuple[int, str]:
    """
    Finds the most specific OrthoDB clade TaxID shared by a list of 3+ species by
    intersecting the `lineage` INTEGER[] arrays across all species in level2species.
    """
    if not species_list:
        raise ValueError("species_list cannot be empty.")

    # Dynamically build JOINs and list_intersect calls for N species
    joins = []
    intersect_expr = "l0.lineage"

    for i in range(1, len(species_list)):
        joins.append(f"JOIN level2species l{i} ON l0.species_code != l{i}.species_code")
        joins.append(f"JOIN species s{i} ON l{i}.species_code = s{i}.species_code")
        intersect_expr = f"list_intersect({intersect_expr}, l{i}.lineage)"

    where_clauses = [f"s{i}.species_name = $species{i}" for i in range(len(species_list))]

    query = f"""
        SELECT DISTINCT UNNEST({intersect_expr}) AS level_tax_id
        FROM level2species l0
        JOIN species s0 ON l0.species_code = s0.species_code
        {" ".join(joins)}
        WHERE {" AND ".join(where_clauses)};
    """

    params = {f"species{i}": name for i, name in enumerate(species_list)}
    
    results = conn.execute(query, params).fetchall()
    shared_odb_levels = {row[0] for row in results if row[0] is not None}
    print('shared odb levels across all species:', species_list, shared_odb_levels)

    if not shared_odb_levels:
        raise ValueError(f"No common OrthoDB clade found across species: {species_list}")

    # Walk NCBI lineage from most specific (index 0) to root
    for node in sdata_lineage:
        tax_id = node[0] if isinstance(node, (tuple, list)) else node
        if tax_id in shared_odb_levels:
            node_name = tid2name.get(
                tax_id, 
                str(node[1] if isinstance(node, (tuple, list)) else tax_id)
            )
            return (tax_id, node_name)


def getOrthologs(
    species1: str, 
    species2: str, 
    mrca_tid: int, 
    search_query: str, 
    limit: int, 
    conn: duckdb.DuckDBPyConnection,
    taxid_to_name: dict[int, str]
) -> dict:
    
    BASE_QUERY = """
    SET enable_progress_bar = true;
    SET memory_limit = '1GB';
    WITH species_a_genes AS (
    SELECT
      g.protein_id AS protein_id_a, 
      g.gene_symbol AS gene_symbol_a,
      g.locus_tag AS locus_tag_a,
      g.replicon AS replicon_a,
      g.coordinates AS coordinates_a,
      g.description AS description_a,
      og.og_id
    FROM species s
    JOIN genes g ON s.species_code = g.species_code
    JOIN og2genes og ON g.odb_gene_id = og.protein_id
    WHERE s.species_name = $species1
      AND og.og_id LIKE $mrca_pattern
    QUALIFY ROW_NUMBER() OVER (
        PARTITION BY 
            og.og_id, 
            COALESCE(
                g.gene_symbol, 
                g.replicon || ':' || g.coordinates, 
                g.odb_gene_id
            )
        ORDER BY 
            g.description IS NOT NULL DESC,  -- Prefer annotated transcripts
            g.protein_id ASC                 -- Deterministic tie-breaker
    ) = 1
),
     species_b_genes AS (
        SELECT
        g.protein_id AS protein_id_b, 
        g.gene_symbol AS gene_symbol_b,
        g.locus_tag AS locus_tag_b,
        g.replicon AS replicon_b,
        g.coordinates AS coordinates_b,
        g.description AS description_b,
        og.og_id
        FROM species s
        JOIN genes g ON s.species_code = g.species_code
        JOIN og2genes og ON g.odb_gene_id = og.protein_id
        WHERE s.species_name = $species2
        AND og.og_id LIKE $mrca_pattern
    QUALIFY ROW_NUMBER() OVER (
        PARTITION BY 
            og.og_id, 
            COALESCE(
                g.gene_symbol, 
                g.replicon || ':' || g.coordinates, 
                g.odb_gene_id
            )
        ORDER BY 
            g.description IS NOT NULL DESC,  -- Prefer annotated transcripts
            g.protein_id ASC                 -- Deterministic tie-breaker
    ) = 1
    ),
    paired_orthologs AS (
        SELECT 
          $species1 AS species_a,
          ap.gene_symbol_a,
          ap.locus_tag_a,
          ap.protein_id_a,
          ap.replicon_a,
          ap.coordinates_a,
          ap.description_a,
          ap.og_id AS orthologous_group,
          bp.gene_symbol_b,
          bp.locus_tag_b,
          bp.protein_id_b,
          bp.replicon_b,
          bp.coordinates_b,
          bp.description_b,
          $species2 AS species_b
        FROM species_a_genes ap
        JOIN species_b_genes bp ON ap.og_id = bp.og_id
        WHERE 
            $search_query = '' OR 
            ap.gene_symbol_a ILIKE $search_pattern OR
            ap.locus_tag_a ILIKE $search_pattern OR
            ap.protein_id_a ILIKE $search_pattern OR
            ap.description_a ILIKE $search_pattern OR
            ap.og_id ILIKE $search_pattern OR
            bp.gene_symbol_b ILIKE $search_pattern OR
            bp.locus_tag_b ILIKE $search_pattern OR
            bp.protein_id_b ILIKE $search_pattern OR
            bp.description_b ILIKE $search_pattern
    ),
    target_ogs AS (
        -- Limit by Ortholog Group count rather than pairwise row count
        SELECT DISTINCT orthologous_group 
        FROM paired_orthologs
        ORDER BY orthologous_group
        LIMIT $limit
    )
    SELECT p.*
    FROM paired_orthologs p
    JOIN target_ogs t ON p.orthologous_group = t.orthologous_group
    ORDER BY p.orthologous_group;
    """

    params = {
        "species1": species1,
        "species2": species2,
        "mrca_pattern": f"%at{mrca_tid}",
        "search_query": search_query,
        "search_pattern": f"%{search_query}%",
        "limit": limit
    }

    cursor = conn.execute(BASE_QUERY, params)
    columns = [desc[0] for desc in cursor.description]
    db_results = [dict(zip(columns, row)) for row in cursor.fetchall()]

    og_dict = defaultdict(lambda: {
        species1: [],
        species2: []
    })

    for row in db_results:
        og_id = row["orthologous_group"]
        og_entry = og_dict[og_id]

        gene_a = {
            "gene_symbol": row.get("gene_symbol_a"),
            "locus_tag": row.get("locus_tag_a"),
            "replicon:coords": (row.get("replicon_a") or '')+":"+row.get("coordinates_a"),
            "description": row.get("description_a")
        }
        if not any(g["locus_tag"] == gene_a["locus_tag"] for g in og_entry[species1]):
            og_entry[species1].append(gene_a)

        gene_b = {
            "gene_symbol": row.get("gene_symbol_b"),
            "locus_tag": row.get("locus_tag_b"),
            "replicon:coords": (row.get("replicon_b") or '')+':'+row.get("coordinates_b"),
            "description": row.get("description_b")
        }
        if not any(g["locus_tag"] == gene_b["locus_tag"] for g in og_entry[species2]):
            og_entry[species2].append(gene_b)

    return dict(og_dict)


def getLocalSimilarities(
    species1: str, 
    species2: str, 
    speciesX: str,
    mrcall: int, 
    search: str, 
    limit: int, 
    conn: duckdb.DuckDBPyConnection
) -> dict:
    print('getLocalSimilarities', species1, species2, speciesX, mrcall)
    """
    Finds orthologs shared between species1 and species2, excluding those 
    present in speciesX, and applies an optional search filter and limit.
    Returns nested dictionary matching getOrthologs output structure.
    """
    subtractive_query = """
    WITH species_a_genes AS (
        SELECT
            g.protein_id AS protein_id_a, 
            g.gene_symbol AS gene_symbol_a,
            g.locus_tag AS locus_tag_a,
            g.replicon AS replicon_a,
            g.coordinates AS coordinates_a,
            g.description AS description_a,
            og.og_id
        FROM species s
        JOIN genes g ON s.species_code = g.species_code
        JOIN og2genes og ON g.odb_gene_id = og.protein_id
        WHERE s.species_name = $species1
          AND og.og_id LIKE $target_node_pat
        QUALIFY ROW_NUMBER() OVER (
            PARTITION BY 
                og.og_id, 
                COALESCE(
                    g.gene_symbol, 
                    g.replicon || ':' || g.coordinates, 
                    g.odb_gene_id
                )
            ORDER BY 
                g.description IS NOT NULL DESC,
                g.protein_id ASC
        ) = 1
    ),
    species_b_genes AS (
        SELECT
            g.protein_id AS protein_id_b, 
            g.gene_symbol AS gene_symbol_b,
            g.locus_tag AS locus_tag_b,
            g.replicon AS replicon_b,
            g.coordinates AS coordinates_b,
            g.description AS description_b,
            og.og_id
        FROM species s
        JOIN genes g ON s.species_code = g.species_code
        JOIN og2genes og ON g.odb_gene_id = og.protein_id
        WHERE s.species_name = $species2
          AND og.og_id LIKE $target_node_pat
        QUALIFY ROW_NUMBER() OVER (
            PARTITION BY 
                og.og_id, 
                COALESCE(
                    g.gene_symbol, 
                    g.replicon || ':' || g.coordinates, 
                    g.odb_gene_id
                )
            ORDER BY 
                g.description IS NOT NULL DESC,
                g.protein_id ASC
        ) = 1
    ),
    exclude_species_ogs AS (
        SELECT DISTINCT og.og_id
        FROM species s_exc
        JOIN genes g_exc ON s_exc.species_code = g_exc.species_code
        JOIN og2genes og ON g_exc.odb_gene_id = og.protein_id
        WHERE s_exc.species_name = $speciesX
    ),
    paired_orthologs AS (
        SELECT 
            $species1 AS species_a,
            ap.gene_symbol_a,
            ap.locus_tag_a,
            ap.protein_id_a,
            ap.replicon_a,
            ap.coordinates_a,
            ap.description_a,
            ap.og_id AS orthologous_group,
            bp.gene_symbol_b,
            bp.locus_tag_b,
            bp.protein_id_b,
            bp.replicon_b,
            bp.coordinates_b,
            bp.description_b,
            $species2 AS species_b
        FROM species_a_genes ap
        JOIN species_b_genes bp ON ap.og_id = bp.og_id
        WHERE ap.og_id NOT IN (SELECT og_id FROM exclude_species_ogs)
          AND (
            $search_query = '' OR 
            ap.gene_symbol_a ILIKE $search_pattern OR
            ap.locus_tag_a ILIKE $search_pattern OR
            ap.protein_id_a ILIKE $search_pattern OR
            ap.description_a ILIKE $search_pattern OR
            ap.og_id ILIKE $search_pattern OR
            bp.gene_symbol_b ILIKE $search_pattern OR
            bp.locus_tag_b ILIKE $search_pattern OR
            bp.protein_id_b ILIKE $search_pattern OR
            bp.description_b ILIKE $search_pattern
          )
    ),
    target_ogs AS (
        SELECT DISTINCT orthologous_group 
        FROM paired_orthologs
        ORDER BY orthologous_group
        LIMIT $limit
    )
    SELECT p.*
    FROM paired_orthologs p
    JOIN target_ogs t ON p.orthologous_group = t.orthologous_group
    ORDER BY p.orthologous_group;
    """

    target_node_pat = "%at" + str(mrcall)
    search_term = search if search is not None else ""
    search_pattern = f"%{search_term}%"

    params = {
        "species1": species1,
        "species2": species2,
        "speciesX": speciesX,
        "target_node_pat": target_node_pat,
        "search_query": search_term,
        "search_pattern": search_pattern,
        "limit": limit
    }

    try:
        cursor = conn.execute(subtractive_query, params)
        columns = [desc[0] for desc in cursor.description]
        db_results = [dict(zip(columns, row)) for row in cursor.fetchall()]

        og_dict = defaultdict(lambda: {
            species1: [],
            species2: []
        })

        for row in db_results:
            og_id = row["orthologous_group"]
            og_entry = og_dict[og_id]

            gene_a = {
                "gene_symbol": row.get("gene_symbol_a"),
                "locus_tag": row.get("locus_tag_a"),
                "replicon:coords": (row.get("replicon_a") or '') + ":" + (row.get("coordinates_a") or ''),
                "description": row.get("description_a")
            }
            if not any(g["locus_tag"] == gene_a["locus_tag"] for g in og_entry[species1]):
                og_entry[species1].append(gene_a)

            gene_b = {
                "gene_symbol": row.get("gene_symbol_b"),
                "locus_tag": row.get("locus_tag_b"),
                "replicon:coords": (row.get("replicon_b") or '') + ":" + (row.get("coordinates_b") or ''),
                "description": row.get("description_b")
            }
            if not any(g["locus_tag"] == gene_b["locus_tag"] for g in og_entry[species2]):
                og_entry[species2].append(gene_b)

        return dict(og_dict)

    except Exception as e:
        print(f"Error executing getLocalSimilarities: {e}")
        return {}

def getLineage(species: str, conn: duckdb.DuckDBPyConnection) -> List[tuple]:
    # Force lowercase comparison and restrict to actual scientific names
    tax_query = """
        SELECT tax_id 
        FROM names 
        WHERE LOWER(name) = LOWER(?) 
          AND name_class = 'scientific name' 
        LIMIT 1;
    """
    print(f"DEBUG: Looking up tax_id for species: '{species}'")
    
    # Pass the raw species string directly without manual string chopping
    res = conn.execute(tax_query, [species]).fetchone()
    
    if not res:
        print(f"DEBUG: No scientific tax_id found for '{species}' in names table!")
        return []
    
    start_tax_id = res[0]
    print(f"DEBUG: Found scientific start_tax_id {start_tax_id} for {species}")
    
    # 2. Walk up the tree using a recursive CTE 
    # (Assuming you also loaded nodes.dmp into a 'nodes' table with tax_id and parent_tax_id)
    lineage_query = """
        WITH RECURSIVE lineage_tree AS (
            -- Base case: start at the given tax_id, seed the path array
            SELECT 
                tax_id, 
                parent_tax_id, 
                0 AS depth, 
                [tax_id] AS path
            FROM nodes
            WHERE tax_id = ?
            
            UNION ALL
            
            -- Recursive step: check if parent_tax_id is already in our path array
            SELECT 
                n.tax_id, 
                n.parent_tax_id, 
                lt.depth + 1,
                list_prepend(n.tax_id, lt.path)
            FROM nodes n
            JOIN lineage_tree lt ON lt.parent_tax_id = n.tax_id
            WHERE NOT list_contains(lt.path, n.tax_id)
              AND lt.depth < 50
        )
        SELECT lt.tax_id, nm.name
        FROM lineage_tree lt
        JOIN names nm ON lt.tax_id = nm.tax_id
        WHERE nm.name_class = 'scientific name'
        ORDER BY lt.depth ASC;
    """
    
    results = conn.execute(lineage_query, [start_tax_id]).fetchall()
    lineage = results
    print("Lineage:", lineage)
    return lineage

def species_gene_count(species: str, conn: duckdb.DuckDBPyConnection) -> int:
    """
    Returns the total number of genes mapped to a given species name 
    in the OrthoDB database, supporting clean integration into the UI.
    """
    count_query = """
        SELECT COUNT(g.locus_tag) 
        FROM species s
        LEFT JOIN genes g ON s.species_code = g.species_code
        WHERE LOWER(s.species_name) = LOWER(?);
    """
    try:
        result = conn.execute(count_query, [species]).fetchone()
        # Return the count integer, defaulting to 0 if none are found or row is null
        return result[0] if result and result[0] is not None else 0
    except Exception as e:
        print(f"Error fetching gene count for '{species}': {e}")
        return 0


def inspect_duckdb(conn):
    """Inspects tables, schemas, and row counts for DuckDB."""
    print("DuckDB Database Schema & Row Counts:")
    tables_query = """
        SELECT table_name, table_type 
        FROM information_schema.tables 
        WHERE table_schema = 'main';
    """
    tables = conn.execute(tables_query).fetchall()
    
    for table_name, table_type in tables:
        # Safely fetch row count for each table
        try:
            row_count = conn.execute(f"SELECT COUNT(*) FROM {table_name};").fetchone()[0]
        except Exception:
            row_count = "N/A (View or virtual table)"
            
        print(f"Table: {table_name} ({table_type}) | Rows: {row_count}")

# Enable CORS matching express cors()
app.middleware("http")
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Database Setup & Query Preparation
NCBI_SQLITE_DB_PATH = "../ncbi_orthologs.db"

ORTHOLOGS_QUERY_TEXT = """
WITH variables AS (
  SELECT 
    (SELECT tax_id FROM names WHERE name LIKE ? LIMIT 1) AS t1,
    (SELECT tax_id FROM names WHERE name LIKE ? LIMIT 1) AS t2
)
SELECT 
  gi1.symbol AS symbol1, 
  gi2.symbol AS symbol2 
FROM orthologs o
CROSS JOIN variables v
JOIN genes AS gi1 ON (o.gene_id = gi1.gene_id AND gi1.tax_id = v.t1)
                  OR (o.gene_id_2 = gi1.gene_id AND gi1.tax_id = v.t1)
JOIN genes AS gi2 ON (o.gene_id_2 = gi2.gene_id AND gi2.tax_id = v.t2)
                  OR (o.gene_id = gi2.gene_id AND gi2.tax_id = v.t2)
WHERE (o.tax_id = v.t1 AND o.tax_id_2 = v.t2)
  OR (o.tax_id = v.t2 AND o.tax_id_2 = v.t1);
"""

# Pydantic models for incoming payload validation
class ImageFile(BaseModel):
    fileData: Optional[str] = None
    species: Optional[str] = None

class OrthologsRequest(BaseModel):
    image1: ImageFile
    image2: ImageFile
    imageX: Optional[ImageFile] = None
    database: Optional[str] = "orthodb"
    search: Optional[str] = Field(default="", description="Search string to filter rows")

def calculate_base64_byte_size(data_url_string: str) -> Optional[float]:
    """Calculates the exact byte size of a file from its Base64 Data URL string."""
    base64_marker = ";base64,"
    marker_index = data_url_string.find(base64_marker)
    
    if marker_index == -1:
        return None

    base64_content = data_url_string[marker_index + len(base64_marker):]
    
    padding_bytes = 0
    if base64_content.endswith("=="):
        padding_bytes = 2
    elif base64_content.endswith("="):
        padding_bytes = 1

    return (len(base64_content) * 3) / 4 - padding_bytes

def extract_base64_file_data(data_url_string: str) -> str:
    """Extracts the raw, unadorned Base64 data string from a Data URL."""
    base64_marker = ";base64,"
    marker_index = data_url_string.find(base64_marker)
    
    if marker_index == -1:
        raise ValueError("Invalid base64 data URL format")

    return data_url_string[marker_index + len(base64_marker):]

@app.post("/orthologs")
async def process_orthologs(payload: OrthologsRequest):
    async def event_generator():
        organisms = [payload.image1, payload.image2]
        if payload.imageX:
            organisms.append(payload.imageX)

        sdatas = []
        # Phase 1 Progress: Vision model / Taxonomy lookup
        yield json.dumps({'step': 'species_identification'}) + "\n\n"

        for organism in organisms:
            if organism.species:
                sdatas.append({"scientificName": organism.species})
            else:
                raw_base64 = extract_base64_file_data(organism.fileData)
                prompt = "Analyze this image. Break down your response into your freeform impressions, and your final specific scientific taxonomy classification (genus species)."
                
                # Run sync Ollama call in thread so loop stays alive
                response = await asyncio.to_thread(
                    ollama.generate,
                    model='qwen2.5vl',
                    prompt=prompt,
                    images=[raw_base64],
                    format={
                        "type": "object",
                        "properties": {
                            "scientificName": {"type": "string"},
                            "impressions": {"type": "string"}
                        },
                        "required": ["scientificName", "impressions"]
                    }
                )
                sdatas.append(json.loads(response.response))

        for sdata in sdatas:
            sdata['scientificName'] = sdata['scientificName'].capitalize()

        yield json.dumps({'step': 'lineage_resolution'}) + "\n\n"

        extras = {}
        odb_name_mapping = {
            "Dictyostelium discoideum": "Dictyostelium discoideum AX4",
            "Gorilla gorilla": "Gorilla gorilla gorilla"
        }
        def odbName(name):
            return odb_name_mapping.get(name, name)

        # Execute DuckDB tasks with active progress polling
        if payload.database == 'ncbi':
            # SQLite non-progress fallback
            yield json.dumps({'step': 'querying_sqlite'})+"\n\n"
            with sqlite3.connect(NCBI_DB_PATH) as conn:
                conn.row_factory = sqlite3.Row  
                cursor = conn.cursor()
                cursor.execute(ORTHOLOGS_QUERY_TEXT, (f"%{sdatas[0]['scientificName']}%", f"%{sdatas[0]['scientificName']}%"))
                db_results = [dict(row) for row in cursor.fetchall()]
        else:
            with duckdb.connect(NCBI_DB_PATH, read_only=True) as conn:
                for sdata in sdatas:
                    sdata['lineage'] = getLineage(sdata['scientificName'], conn)
                    yield json.dumps({"step": "lineage", "lineage": sdata['lineage']}, indent=2) + "\n\n"
                lin1 = sdatas[0]['lineage'][::-1]
                lin2 = sdatas[1]['lineage'][::-1]
                shared_nodes = [a for a, b in zip(lin1, lin2) if a == b]
                extras['mrca'] = shared_nodes[-1] if shared_nodes else None
                yield json.dumps({"step": "mrca", "mrca": extras['mrca']})

            # 2. Main Query Execution with Progress Tracking
            with duckdb.connect(ORTHODB_DB_PATH, read_only=True) as conn:
                taxid_to_name = {}
                for sdata in sdatas:
                    sdata['gene_count'] = species_gene_count(odbName(sdata['scientificName']), conn)
                    for n, s in sdata['lineage']:
                        taxid_to_name[n] = s

                extras['odb_mrca'] = odb_mrca(
                    conn,
                    odbName(sdatas[0]['scientificName']), 
                    odbName(sdatas[1]['scientificName']),
                    sdatas[0]['lineage'],
                    taxid_to_name
                )

                if len(sdatas) > 2:
                    extras['odb_mrcall'] = odb_mrcall(
                        conn,
                        [odbName(s['scientificName']) for s in sdatas],
                        sdatas[0]['lineage'],
                        taxid_to_name
                    )
                    target_fn = lambda: getLocalSimilarities(
                        odbName(sdatas[0]['scientificName']),
                        odbName(sdatas[1]['scientificName']),
                        odbName(sdatas[2]['scientificName']),
                        extras['odb_mrcall'][0],
                        payload.search or "",
                        20000,
                        conn
                    )
                else:
                    target_fn = lambda: getOrthologs(
                        odbName(sdatas[0]['scientificName']),
                        odbName(sdatas[1]['scientificName']),
                        extras['odb_mrca'][0],
                        payload.search or "",
                        20000,
                        conn,
                        taxid_to_name
                    )

                task = asyncio.create_task(asyncio.to_thread(target_fn))
                while not task.done():
                    yield json.dumps({'step': 'get_orthologs', 'progress': conn.query_progress()})+"\n\n"
                    await asyncio.sleep(1)  # poll every 200ms
                db_results = await task

        for sdata in sdatas:
            sdata['lineage'] = [line[1] for line in sdata['lineage']]

        # Yield Final Payload
        final_payload = {
            "type": "result",
            "message": f"Analysis complete using database: {payload.database}",
            "sdatas": sdatas,
            "extras": extras,
            "dbResults": db_results
        }
        yield json.dumps(final_payload, indent=2)

    return StreamingResponse(event_generator(), media_type="text/event-stream")

@app.get("/")
async def serve_frontend():
    html_path = "../barnyard.html"
    if os.path.exists(html_path):
        return FileResponse(html_path)
    return {"error": "barnyard.html not found"}

@app.get("/GBbackgroundimage.jpg")
async def serve_frontend():
    file_path = "../GBbackgroundimage.jpg"
    if os.path.exists(file_path):
        return FileResponse(file_path)
    return {"error": "GBbackgroundimage.jpg not found"}

@app.get("/barn.jpg")
async def serve_frontend():
    file_path = "../barn.jpg"
    if os.path.exists(file_path):
        return FileResponse(file_path)
    return {"error": "barn.jpg not found"}

PICTURES_DIR = "../pictures"
app.mount("/pictures", StaticFiles(directory=PICTURES_DIR), name="picture")
@app.get("/pictures.json")
async def list_pictures():
    return [f for f in os.listdir(PICTURES_DIR)
        if os.path.isfile(os.path.join(PICTURES_DIR, f))
    ]

if __name__ == "__main__":
    import uvicorn
    uvicorn.run("main:app", host="0.0.0.0", port=3000, reload=True)