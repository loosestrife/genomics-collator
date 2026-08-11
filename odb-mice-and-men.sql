WITH species_a_proteins AS (
  SELECT 
    g.og_level, 
    g.locus_tag AS locus_tag_a,
    g.replicon AS replicon_a,
    og.og_id
  FROM species s
  JOIN genes g ON s.species_code = g.species_code
  JOIN og2genes og ON g.og_level = og.protein_id
  WHERE s.species_name = 'Homo sapiens'
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
  WHERE s.species_name = 'Mus musculus'
)
SELECT 
  'Homo sapiens' AS species_a,
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
LIMIT 50;