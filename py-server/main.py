import os
import traceback
from typing import List, Optional
import duckdb
from fastapi import FastAPI, HTTPException, status
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field
import sqlite3
import ollama

app = FastAPI(title="NCBI Orthologs FastAPI Server")

NCBI_DB_PATH="../ncbi.duckdb"
ORTHODB_DB_PATH="../orthodb.duckdb"

def getOrthologs(species1: str, species2: str, search_query: str, limit: int, conn: duckdb.DuckDBPyConnection):
    print('getOrthologs', species1, species2)
    ORTHODB_DB_QUERY = """
    WITH species_a_proteins AS (
        SELECT 
          g.og_level, 
          g.locus_tag AS locus_tag_a,
          g.replicon AS replicon_a,
          og.og_id
        FROM species s
        JOIN genes g ON s.species_code = g.species_code
        JOIN og2genes og ON g.og_level = og.protein_id
        WHERE s.species_name = ?
    ),
    species_b_proteins AS (
        SELECT 
          g.og_level, 
          g.locus_tag AS locus_tag_b,
          g.replicon AS replicon_b,
          g.description,
          og.og_id,
          s.species_name
        FROM species s
        JOIN genes g ON s.species_code = g.species_code
        JOIN og2genes og ON g.og_level = og.protein_id
        WHERE s.species_name = ?
    ),
    base_query AS (
        SELECT DISTINCT ON (ap.locus_tag_a, bp.locus_tag_b)
          ? AS species_a,
          ap.locus_tag_a,
          ap.og_level AS protein_id_a,
          ap.replicon_a,
          ap.og_id AS orthologous_group,
          bp.locus_tag_b,
          bp.og_level AS protein_id_b,
          bp.replicon_b,
          bp.description AS description_b,
          bp.species_name AS species_b
        FROM species_a_proteins ap
        JOIN species_b_proteins bp ON ap.og_id = bp.og_id
    )
    SELECT * FROM base_query
    WHERE 
        ? = '' OR 
        CAST(species_a AS VARCHAR) ILIKE ? OR
        CAST(locus_tag_a AS VARCHAR) ILIKE ? OR
        CAST(protein_id_a AS VARCHAR) ILIKE ? OR
        CAST(replicon_a AS VARCHAR) ILIKE ? OR
        CAST(orthologous_group AS VARCHAR) ILIKE ? OR
        CAST(locus_tag_b AS VARCHAR) ILIKE ? OR
        CAST(protein_id_b AS VARCHAR) ILIKE ? OR
        CAST(replicon_b AS VARCHAR) ILIKE ? OR
        CAST(description_b AS VARCHAR) ILIKE ? OR
        CAST(species_b AS VARCHAR) ILIKE ?
    ORDER BY replicon_a
    LIMIT ?;
    """
    # Parameters mapping:
    # 1, 2: Species names for CTEs
    # 3: Species A label for SELECT
    # 4: Empty check for search (? = '')
    # 5-14: ILIKE patterns for each of the 10 output columns
    # 15: LIMIT value
    search_pattern = "%"+search_query+"%"
    params = [
        species1, species2, species1, 
        search_query, 
        search_pattern, search_pattern, search_pattern, search_pattern, 
        search_pattern, search_pattern, search_pattern, search_pattern, 
        search_pattern, search_pattern, 
        limit
    ]
    cursor = conn.execute(ORTHODB_DB_QUERY, params)
    columns = [desc[0] for desc in cursor.description]
    db_results = [dict(zip(columns, row)) for row in cursor.fetchall()]
    return db_results

def getLocalSimilarities(
    species1: str, 
    species2: str, 
    speciesX: str,
    mrclist: int, 
    search: str, 
    limit: int, 
    conn: duckdb.DuckDBPyConnection
) -> List[any]:
    """
    Finds orthologs shared between species1 and species2, excluding those 
    present in speciesX, and applies an optional search filter and limit.
    """
    print("getLocalSimilarities", species1, species2, speciesX, mrclist, search, limit)
    subtractive_query = """
WITH species_a_proteins AS (
        SELECT 
          g.og_level, 
          g.locus_tag AS locus_tag_a,
          g.replicon AS replicon_a,
          og.og_id
        FROM species s
        JOIN genes g ON s.species_code = g.species_code
        JOIN og2genes og ON g.og_level = og.protein_id
        WHERE s.species_name = ?
    ),
    species_b_proteins AS (
        SELECT 
          g.og_level, 
          g.locus_tag AS locus_tag_b,
          g.replicon AS replicon_b,
          g.description,
          og.og_id,
          s.species_name
        FROM species s
        JOIN genes g ON s.species_code = g.species_code
        JOIN og2genes og ON g.og_level = og.protein_id
        WHERE s.species_name = ?
    ),
    exclude_species_proteins AS (
        SELECT DISTINCT og.og_id
        FROM species s_exc
        JOIN genes g_exc ON s_exc.species_code = g_exc.species_code
        JOIN og2genes og ON g_exc.og_level = og.protein_id
        WHERE s_exc.species_name = ?
    ),
    base_query AS (
        SELECT DISTINCT ON (ap.locus_tag_a, bp.locus_tag_b) 
          ap.og_id AS orthologous_group,
          ? AS species_a,
          ap.locus_tag_a,
          ap.og_level AS protein_id_a,
          ap.replicon_a,
          bp.locus_tag_b,
          bp.og_level AS protein_id_b,
          bp.replicon_b,
          bp.description AS description_b,
          bp.species_name AS species_b
        FROM species_a_proteins ap
        JOIN species_b_proteins bp ON ap.og_id = bp.og_id
        WHERE ap.og_id NOT IN (SELECT og_id FROM exclude_species_proteins)
          -- Updated to match a comma-separated list of target taxonomic node IDs
          -- e.g., 'at314145,at33554,at32523'
          AND regexp_matches(ap.og_id, ?)
    )
    SELECT * FROM base_query
    WHERE 
        ? = '' OR 
        CAST(species_a AS VARCHAR) ILIKE ? OR
        CAST(locus_tag_a AS VARCHAR) ILIKE ? OR
        CAST(protein_id_a AS VARCHAR) ILIKE ? OR
        CAST(replicon_a AS VARCHAR) ILIKE ? OR
        CAST(orthologous_group AS VARCHAR) ILIKE ? OR
        CAST(locus_tag_b AS VARCHAR) ILIKE ? OR
        CAST(protein_id_b AS VARCHAR) ILIKE ? OR
        CAST(replicon_b AS VARCHAR) ILIKE ? OR
        CAST(description_b AS VARCHAR) ILIKE ? OR
        CAST(species_b AS VARCHAR) ILIKE ?
    ORDER BY replicon_a
    LIMIT ?;
    """
    
    target_nodes_regex = "at(" + '|'.join([str(x) for x in mrclist]) + ")"
    search_term = search if search is not None else ""
    search_pattern = f"%{search_term}%"

    params = [
        species1,            # 1. species_a_proteins
        species2,            # 2. species_b_proteins
        speciesX,            # 3. exclude_species_proteins (subtraction target)
        species1,            # 4. base_query SELECT species_a
        target_nodes_regex,        # 5. ap.og_id LIKE ?
        search_term,         # 5. WHERE ? = ''
        search_pattern,      # 6. species_a ILIKE
        search_pattern,      # 7. locus_tag_a ILIKE
        search_pattern,      # 8. protein_id_a ILIKE
        search_pattern,      # 9. replicon_a ILIKE
        search_pattern,      # 10. orthologous_group ILIKE
        search_pattern,      # 11. locus_tag_b ILIKE
        search_pattern,      # 12. protein_id_b ILIKE
        search_pattern,      # 13. replicon_b ILIKE
        search_pattern,      # 14. description_b ILIKE
        search_pattern,      # 15. species_b ILIKE
        limit                # 16. LIMIT
    ]

    try:
        results = conn.execute(subtractive_query, params).fetchall()
        return results
    except Exception as e:
        print(f"Error executing getLocalSimilarities: {e}")
        return []

def getLineage(species: str, conn: duckdb.DuckDBPyConnection) -> List[str]:
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
    
    try:
        results = conn.execute(lineage_query, [start_tax_id]).fetchall()
        lineage = results
        print("Lineage:", lineage)
        return lineage
    except Exception as e:
        print(f"Lineage traversal error (is nodes table missing?): {e}")
        return [None, species] # Fallback gracefully

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
    print("received call to /orthologs endpoint. Selected database:", payload.database)
    organisms = [payload.image1, payload.image2]
    if payload.imageX:
        organisms.append(payload.imageX)
    try:
        sdatas = []
        for organism in organisms:
            if organism.species:
                sdatas.append({"scientificName": organism.species})
            else:
                print('determining species from image... ... ...')
                ollama_response_schema = {
                    "type": "object",
                    "properties": {
                        "scientificName": { 
                            "type": "string", 
                            "description": 'The clean Genus species name only, formatted as "Genus species".' 
                        },
                        "impressions": { 
                            "type": "string", 
                            "description": "Your freeform thoughts, observations, and context about the specimen." 
                        }
                    },
                    "required": ["scientificName", "impressions"]
                }

                raw_base64 = extract_base64_file_data(organism.fileData)
                prompt = "Analyze this image. Break down your response into your freeform impressions, and your final specific scientific taxonomy classification (genus species)."
                
                response = ollama.generate(
                    model='qwen2.5vl',
                    prompt=prompt,
                    images=[raw_base64],
                    format=ollama_response_schema,
                )
                print(response)
                
                import json
                sdatas.append(json.loads(response.response))

        for sdata in sdatas:
            sdata['scientificName'] = sdata['scientificName'].capitalize()
        print("got names", sdatas)

        db_results = []
        extras = {}
        # Explicit routing based on frontend payload choice
        if payload.database == 'ncbi':
            with sqlite3.connect(NCBI_DB_PATH) as conn:
                conn.row_factory = sqlite3.Row  
                cursor = conn.cursor()
                cursor.execute(ORTHOLOGS_QUERY_TEXT, (f"%{raw_name1}%", f"%{raw_name2}%"))
                rows = cursor.fetchall()
                db_results = [dict(row) for row in rows]
        else:
            search_query = payload.search or ""
            search_pattern = f"%{search_query}%"
            
            with duckdb.connect(NCBI_DB_PATH, read_only=True) as conn:
                inspect_duckdb(conn)
                for sdata in sdatas:
                    sdata['lineage'] = getLineage(sdata['scientificName'], conn)
                matching_nodes = [pair[0] for pair in zip(
                    sdatas[0]['lineage'][::-1], sdatas[1]['lineage'][::-1]
                ) if pair[0] == pair[1]]
                extras['mrca'] = matching_nodes[-1] if matching_nodes else None                    
                print("MCRA calculation", sdatas[0]['scientificName'], sdatas[1]['scientificName'], extras['mrca'])
                if len(sdatas)>2:
                    matching_nodes = [trip[0] for trip in zip(
                        sdatas[0]['lineage'][::-1], sdatas[1]['lineage'][::-1], sdatas[2]['lineage'][::-1]
                    ) if trip[0] == trip[1] and trip[1] == trip[2]]
                    extras['mrcall'] = matching_nodes[-1] if matching_nodes else None
                    print("MRCA of all species", extras['mrcall'])
                    extras['mrclist'] = [trip[0] for trip in zip(
                        sdatas[0]['lineage'][::-1], sdatas[1]['lineage'][::-1], sdatas[2]['lineage'][::-1]
                    ) if trip[0] == trip[1] and trip[1] != trip[2]] + [extras['mrcall']]
                    



            with duckdb.connect(ORTHODB_DB_PATH, read_only=True) as conn:
                inspect_duckdb(conn)
                for sdata in sdatas:
                    sdata['gene_count'] = species_gene_count(sdata['scientificName'], conn)
                if(len(sdatas) > 2):
                    db_results = getLocalSimilarities(
                        sdatas[0]['scientificName'],
                        sdatas[1]['scientificName'],
                        sdatas[2]['scientificName'],
                        [x[0] for x in extras['mrclist']],
                        search_query,
                        20000,
                        conn)
                else: 
                    db_results = getOrthologs(
                        sdatas[0]['scientificName'],
                        sdatas[1]['scientificName'],
                        search_query,
                        20000,
                        conn)


        print(f"got {len(db_results)} results from {payload.database}")
        extras['dbResultsLen'] = len(db_results)

        for sdata in sdatas:
            sdata['lineage'] = [line[1] for line in sdata['lineage']]
        if 'mrclist' in extras:
            del extras['mrclist']

        return {
            "message": f"Analysis complete using database: {payload.database}",
            "sdatas": sdatas,
            "extras": extras,
            "dbResults": db_results
        }

    except Exception as error:
        traceback.print_exc()
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail={"error": str(error), "sdatas": sdatas}
        )

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