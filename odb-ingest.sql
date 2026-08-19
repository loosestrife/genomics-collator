-- 1. Load Species Table
CREATE TABLE species AS 
SELECT * FROM read_csv(
    'data/odb12v2_species.tab.gz', 
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
    'data/odb12v2_genes.tab.gz', 
    delim = '\t', 
    header = false, 
    columns = {
        'odb_gene_id':  'VARCHAR',  -- Col 1: Internal OrthoDB Gene ID (e.g., 100_0:000000)
        'species_code': 'VARCHAR',  -- Col 2: Species Code (e.g., 100_0)
        'protein_id':   'VARCHAR',  -- Col 3: RefSeq/NCBI Accession (e.g., WP_131833652.1)
        'locus_tag':    'VARCHAR',  -- Col 4: Locus Tag
        'uniprot_id':   'VARCHAR',  -- Col 5: UniProt ID
        'gene_symbol':  'VARCHAR',  -- Col 6: Public Gene Name / Symbol
        'alt_id':       'VARCHAR',  -- Col 7: Alternate / Secondary Identifier
        'description':  'VARCHAR',  -- Col 8: Functional Description
        'coordinates':  'VARCHAR',  -- Col 9: Genomic Coordinates [start:end](strand)
        'replicon':     'VARCHAR',  -- Col 10: Contig / Chromosome Accession
        'source':       'VARCHAR'   -- Col 11: Data Source / Status
    }
);

-- 3. Load OG-to-Genes Table
CREATE TABLE og2genes AS 
SELECT * FROM read_csv(
    'data/odb12v2_OG2genes.tab.gz', 
    delim = '\t', 
    header = false, 
    columns = {
        'og_id': 'VARCHAR', 
        'protein_id': 'VARCHAR'
    }
);

CREATE TABLE level2species AS 
SELECT 
    column0::INTEGER AS level_tax_id,
    column1 AS species_code,
    column2::INTEGER AS status,
    string_split(trim(both '{}' from column3), ',')::INTEGER[] AS lineage
FROM read_csv(
    'data/odb12v2_level2species.tab.gz', 
    delim = '\t', 
    header = false
);