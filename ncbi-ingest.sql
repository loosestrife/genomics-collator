-- Connect or create your DuckDB file
-- ATTACH 'ncbi_orthologs.duckdb' AS ncbi;

-- 1. Load Names Table (from pipe-delimited names.dmp)
CREATE TABLE names AS 
SELECT 
    CAST(column0 AS INTEGER) AS tax_id, 
    REGEXP_REPLACE(column1, '^[[:space:]]+|[[:space:]]+$', '', 'g') AS name,
    REGEXP_REPLACE(column3, '^[[:space:]]+|[[:space:]]+$', '', 'g') AS name_class
FROM read_csv(
    'names.dmp', 
    header = false, 
    delim = '|', 
    quote = '', 
    all_varchar = true
) 
WHERE column0 NOT LIKE '#%';


-- 2. Load Genes Table (from gzipped tab-delimited gene_info.gz)
CREATE TABLE genes AS 
SELECT 
    CAST(column00 AS INTEGER) AS tax_id,
    CAST(column01 AS UBIGINT) AS gene_id,
    column02 AS symbol,
    column04 AS synonyms,
    column08 AS description
FROM (
    SELECT * FROM read_csv(
        'gene_info.gz', 
        header = false, 
        delim = '\t', 
        quote = '', 
        nullstr = '-', 
        all_varchar = true,
        ignore_errors = true
    )
) 
WHERE column00 NOT LIKE '#%';

-- 3. Load Orthologs Table (from gzipped tab-delimited gene_orthologs.gz)
CREATE TABLE orthologs AS 
SELECT 
    CAST(column0 AS INTEGER) AS tax_id,
    CAST(column1 AS UBIGINT) AS gene_id,
    CAST(column3 AS INTEGER) AS tax_id_2,
    CAST(column4 AS UBIGINT) AS gene_id_2
FROM (
    SELECT * FROM read_csv(
        'gene_orthologs.gz', 
        header = false, 
        delim = '\t', 
        quote = '', 
        all_varchar = true,
        ignore_errors = true
    )
) 
WHERE column0 NOT LIKE '#%';


CREATE TABLE nodes AS 
SELECT 
    CAST(column00 AS INTEGER) AS tax_id, 
    CAST(column01 AS INTEGER) AS parent_tax_id,
    REGEXP_REPLACE(column02, '^[[:space:]]+|[[:space:]]+$', '', 'g') AS rank
FROM read_csv(
    'nodes.dmp', 
    header = false, 
    delim = '|', 
    quote = '', 
    all_varchar = true
);