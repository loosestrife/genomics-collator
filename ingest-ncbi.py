import sqlite3
import gzip

# Connect to your SQLite database file
conn = sqlite3.connect("ncbi_orthologs.db")
cursor = conn.cursor()

# 1. Execute schema creation
cursor.executescript("""
CREATE TABLE IF NOT EXISTS genes (
    tax_id INTEGER,
    gene_id INTEGER PRIMARY KEY,
    symbol TEXT,
    synonyms TEXT,
    description TEXT
);
CREATE TABLE IF NOT EXISTS names (
    tax_id INTEGER,
    name TEXT
);
CREATE TABLE IF NOT EXISTS orthologs (
    tax_id INTEGER,
    gene_id INTEGER,
    tax_id_2 INTEGER,
    gene_id_2 INTEGER
);
""")

# 2. ingestion
batch_size = 100000

# Ingest gene_orthologs.gz (NCBI gene orthology mapping)
# Format columns: 0:tax_id_1, 1:gene_id_1, 2:tax_id_2, 3:gene_id_2
print("Ingesting organism names")
names_batch = []

with open("names.dmp") as f:
    for i, line in enumerate(f):
        if i%1000000 == 0:
            print("line ", i)
        if line.startswith("#"):
            continue
        parts = [part.strip() for part in line.strip().split("|")]
        if len(parts) >= 2:
            names_batch.append(
                (int(parts[0]), parts[1])
            )

            if len(names_batch) >= batch_size:
                cursor.executemany(
                    "INSERT INTO names VALUES (?, ?)", names_batch
                )
                conn.commit()
                names_batch = []
        else:
            print("malformed line ", line)

if names_batch:
    cursor.executemany("INSERT INTO names VALUES (?, ?)", names_batch)
    conn.commit()
# 4. Create Indexes post-insertion (significantly speeds up the initial load)
print("Creating indexes...")
cursor.executescript("""
CREATE INDEX IF NOT EXISTS idx_genes_symbol ON genes(symbol);
CREATE INDEX IF NOT EXISTS idx_genes_tax ON genes(tax_id);
CREATE INDEX IF NOT EXISTS idx_ortho_g1 ON orthologs(gene_id);
CREATE INDEX IF NOT EXISTS idx_ortho_g2 ON orthologs(gene_id_2);
""")

# Ingest gene_info.gz (NCBI gene metadata)
# Format columns: 0:tax_id, 1:GeneID, 2:Symbol, 4:synonyms, 8:description
print("Ingesting gene_info...")
gene_batch = []

with gzip.open("gene_info.gz", "rt") as f:
    for i, line in enumerate(f):
        if i%1000000 == 0:
            print("line ", i)
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) > 8:
            tax_id, gene_id, symbol, synonyms, description = (
                parts[0],
                parts[1],
                parts[2],
                parts[4],
                parts[8],
            )
            gene_batch.append((int(tax_id), int(gene_id), symbol, synonyms, description))

            if len(gene_batch) >= batch_size:
                cursor.executemany(
                    "INSERT OR IGNORE INTO genes VALUES (?, ?, ?, ?, ?)", gene_batch
                )
                conn.commit()
                gene_batch = []

if gene_batch:
    cursor.executemany("INSERT OR IGNORE INTO genes VALUES (?, ?, ?, ?, ?)", gene_batch)
    conn.commit()

# 3. Ingest gene_orthologs.gz (NCBI gene orthology mapping)
# Format columns: 0:tax_id_1, 1:gene_id_1, 2:tax_id_2, 3:gene_id_2
print("Ingesting gene_orthologs...")
ortho_batch = []

with gzip.open("gene_orthologs.gz", "rt") as f:
    for i, line in enumerate(f):
        if i%1000000 == 0:
            print("line ", i)
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) == 5:
            ortho_batch.append(
                (int(parts[0]), int(parts[1]), int(parts[3]), int(parts[4]))
            )

            if len(ortho_batch) >= batch_size:
                cursor.executemany(
                    "INSERT INTO orthologs VALUES (?, ?, ?, ?)", ortho_batch
                )
                conn.commit()
                ortho_batch = []
        else:
            print("malformed line ", line)

if ortho_batch:
    cursor.executemany("INSERT INTO orthologs VALUES (?, ?, ?, ?)", ortho_batch)
    conn.commit()


conn.close()
print("Database build complete!")