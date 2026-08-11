-- 1. Load Species Table
CREATE TABLE species AS 
SELECT * FROM read_csv(
    'odb12v2_species.tab.gz', 
    delim = '\t', 
    header = false, 
    columns = {
        'tax_id': 'INTEGER', 
        'species_code': 'VARCHAR', 
        'species_name': 'VARCHAR', 
        'assembly_acc': 'VARCHAR',
        'column4': 'VARCHAR',
        'column5': 'VARCHAR',
        'column6': 'VARCHAR'
    }
);

-- 2. Load Genes Table
CREATE TABLE genes AS 
SELECT * FROM read_csv(
    'odb12v2_genes.tab.gz', 
    delim = '\t', 
    header = false, 
    columns = {
        'og_level': 'VARCHAR', 
        'species_code': 'VARCHAR', 
        'protein_id': 'VARCHAR', 
        'locus_tag': 'VARCHAR', 
        'uniprot_id': 'VARCHAR', 
        'description': 'VARCHAR', 
        'coordinates': 'VARCHAR', 
        'replicon': 'VARCHAR',
        'extra': 'VARCHAR',
        'column9': 'VARCHAR',
        'column10': 'VARCHAR',
    }
);

-- 3. Load OG-to-Genes Table
CREATE TABLE og2genes AS 
SELECT * FROM read_csv(
    'odb12v2_OG2genes.tab.gz', 
    delim = '\t', 
    header = false, 
    columns = {
        'og_id': 'VARCHAR', 
        'protein_id': 'VARCHAR'
    }
);
