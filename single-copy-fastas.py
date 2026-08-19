import os
import time
from io import StringIO
from Bio import Entrez, SeqIO
import pandas as pd

Entrez.email = "your_email@domain.com"
df = pd.read_csv("single_copy_orthologs.csv")

os.makedirs("fastas", exist_ok=True)


def fetch_records_validated(protein_ids):
    """Fetches FASTA records from NCBI.

    If NCBI returns an inline error message in the batch, it falls back to
    1-by-1 fetching to isolate and skip bad IDs.
    """
    try:
        handle = Entrez.efetch(
            db="protein", id=protein_ids, rettype="fasta", retmode="text"
        )
        raw_text = handle.read()
        handle.close()
    except Exception as e:
        print(
            f"  [HTTP Error] Batch fetch failed: {e}. Switching to individual fetch..."
        )
        raw_text = ""

    # Check for NCBI error text embedded in the stream
    if (
        not raw_text
        or "CEFetchPApplication" in raw_text
        or "Error:" in raw_text
    ):
        print(
            "  [Warning] Batch contains invalid IDs. Isolating valid IDs individually..."
        )
        valid_records = []

        for pid in protein_ids:
            time.sleep(0.2)  # Respect NCBI rate limits
            try:
                h = Entrez.efetch(
                    db="protein", id=pid, rettype="fasta", retmode="text"
                )
                single_txt = h.read()
                h.close()

                if (
                    "CEFetchPApplication" in single_txt
                    or "Error:" in single_txt
                    or not single_txt.strip()
                ):
                    print(f"  [Skipped Bad ID] {pid} is invalid or deprecated.")
                    continue

                # Parse validated single record
                recs = list(SeqIO.parse(StringIO(single_txt), "fasta"))
                valid_records.extend(recs)
            except Exception as err:
                print(f"  [Failed] Could not fetch {pid}: {err}")

        return valid_records

    # Parse clean batch text directly from memory stream
    return list(SeqIO.parse(StringIO(raw_text), "fasta"))


# Main Loop
for og_id, group in df.groupby("og_id"):
    fasta_path = f"fastas/{og_id}.fa"
    if os.path.exists(fasta_path):
        continue

    protein_ids = group["protein_id"].tolist()
    species_map = dict(zip(group["protein_id"], group["species_code"]))

    print(f"Fetching {len(protein_ids)} sequences for {og_id}...")
    records = fetch_records_validated(protein_ids)

    if not records:
        print(f"  [Warning] No valid records retrieved for {og_id}. Skipping.")
        time.sleep(0.3)
        continue

    clean_records = []
    for rec in records:
        # 1. Sanitize any +Error string artifacts attached to sequence ends
        seq_str = str(rec.seq)
        if "+Error" in seq_str:
            seq_str = seq_str.split("+Error")[0]
        rec.seq = rec.seq.__class__(seq_str)

        # 2. Match ID (handles both versioned 'NP_1234.1' and unversioned 'NP_1234')
        acc = rec.id.split()[0]
        acc_no_version = acc.split(".")[0]

        species_code = species_map.get(acc) or species_map.get(acc_no_version)

        if species_code:
            rec.id = species_code
            rec.description = ""
            clean_records.append(rec)
        else:
            print(
                f"  [Warning] Accession '{acc}' not found in species_map for {og_id}"
            )

    # Only write file if valid sequences exist
    if clean_records:
        SeqIO.write(clean_records, fasta_path, "fasta")

    time.sleep(0.3)