-- Step 1: Extract your t1 and t2 variables cleanly at the top
WITH variables AS (
    SELECT 
        (SELECT tax_id FROM names WHERE name LIKE '%Homo sapiens%' LIMIT 1) AS t1,
        (SELECT tax_id FROM names WHERE name LIKE '%Mus musculus%' LIMIT 1) AS t2
)
-- Step 2: Run your exact mental model using those variables
SELECT 
    gi1.symbol, 
    gi2.symbol 
FROM orthologs o
CROSS JOIN variables v
JOIN genes AS gi1 ON o.gene_id = gi1.gene_id AND gi1.tax_id = v.t1
JOIN genes AS gi2 ON o.gene_id_2 = gi2.gene_id AND gi2.tax_id = v.t2
WHERE o.tax_id = v.t1 AND o.tax_id_2 = v.t2;