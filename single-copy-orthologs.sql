-- Save as extract_single_copy.sql
COPY (
WITH mammal_spp AS (
    SELECT DISTINCT species_code 
    FROM level2species 
    WHERE level_tax_id = 40674 OR list_contains(lineage, 40674)
),
og_species_counts AS (
    SELECT 
        o.og_id,
        g.species_code,
        COUNT(*) AS copy_count
    FROM og2genes o
    JOIN genes g ON o.protein_id = g.odb_gene_id  -- Join on internal odb_gene_id
    WHERE g.species_code IN (SELECT species_code FROM mammal_spp)
    GROUP BY o.og_id, g.species_code
),
valid_single_copy_ogs AS (
    SELECT og_id
    FROM og_species_counts
    GROUP BY og_id
    -- Must be present in at least 30 mammals with strictly 1 copy per species
    HAVING COUNT(species_code) >= 30 
       AND MAX(copy_count) = 1
    LIMIT 200
)
SELECT 
    o.og_id,
    g.species_code,
    s.species_name,
    g.protein_id
FROM og2genes o
JOIN genes g ON o.protein_id = g.odb_gene_id  -- Join on internal odb_gene_id
JOIN species s ON g.species_code = s.species_code
WHERE o.og_id IN (SELECT og_id FROM valid_single_copy_ogs)
  AND g.species_code IN (SELECT species_code FROM mammal_spp)
) TO 'single_copy_orthologs.csv' (HEADER, DELIMITER ',');