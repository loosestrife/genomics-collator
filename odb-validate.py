import gzip
import sys

files = [
    ("odb12v2_species.tab.gz", 7),
    ("odb12v2_genes.tab.gz", 11),
    ("odb12v2_OG2genes.tab.gz", 2),
]

for filepath, expected_cols in files:
    print(f"Scanning {filepath} for structural anomalies...")
    line_count = 0
    malformed_count = 0

    try:
        with gzip.open(filepath, "rt", encoding="utf-8", errors="replace") as f:
            for line_num, line in enumerate(f, 1):
                line_count += 1
                # Split strictly on tabs
                cols = line.rstrip("\n\r").split("\t")

                if len(cols) != expected_cols:
                    malformed_count += 1
                    if malformed_count <= 5:  # Print first 5 errors as samples
                        print(
                            f"  [Line {line_num}] Expected {expected_cols} columns, found {len(cols)}. Content preview: {line[:100].strip()}"
                        )

        print(
            f"Finished {filepath}: {line_count} total lines, {malformed_count} malformed rows found.\n"
        )
    except Exception as e:
        print(f"Error reading {filepath}: {e}", file=sys.stderr)