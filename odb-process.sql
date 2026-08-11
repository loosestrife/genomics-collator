
-- 4. Create Indexes
CREATE INDEX idx_genes_species ON genes(species_code);
CREATE INDEX idx_genes_prot ON genes(protein_id);
CREATE INDEX idx_og2genes_og ON og2genes(og_id);
CREATE INDEX idx_og2genes_prot ON og2genes(protein_id);
CREATE INDEX idx_species_code ON species(species_code);