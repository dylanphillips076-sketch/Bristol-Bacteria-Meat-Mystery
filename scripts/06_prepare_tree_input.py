# =============================================================================
# 06_prepare_tree_input.py
# Bristol Bacteria Meat Mystery - Step 6 (part 1): prepare sequences for the
# protein alignment
#
#   - downloads the protein sequences of the chosen references and outgroup
#     from NCBI (by accession, using Biopython's Entrez module)
#   - gives every sequence a short, informative label for the tree
#   - combines the references, outgroup and the 5 representative proteins
#     into one FASTA file, ready for alignment with Clustal Omega
#
# Input : results/selected_references.csv      (chosen by hand in Step 5)
#         results/representative_proteins.fasta (from Step 4)
# Output: results/references_raw.fasta          downloaded reference proteins
#         results/tree_input_proteins.fasta     all sequences, labelled
#         results/tree_labels.csv               label -> group, species, role
#
# Downloaded references are saved, so the download is not repeated if the
# script is run again (delete references_raw.fasta to download again).
#
# Run from the repository root:  python scripts/06_prepare_tree_input.py
# =============================================================================

# %% 1. Packages and file paths
import re
from pathlib import Path

import pandas as pd
from Bio import Entrez, SeqIO

ROOT = Path.cwd()
while not (ROOT / "results").exists() and ROOT != ROOT.parent:
    ROOT = ROOT.parent
RESULTS = ROOT / "results"

SELECTED = RESULTS / "selected_references.csv"
REPRESENTATIVES = RESULTS / "representative_proteins.fasta"
REFERENCES_RAW = RESULTS / "references_raw.fasta"
TREE_INPUT = RESULTS / "tree_input_proteins.fasta"
LABELS = RESULTS / "tree_labels.csv"

Entrez.email = "your.email@example.com"   # NCBI asks for an email address with downloads
print("Repository root:", ROOT)


# %% 2. Read the chosen references
refs = pd.read_csv(SELECTED)
refs["accession"] = refs["accession"].str.strip()
print(refs[["group", "accession", "species", "role"]].to_string(index=False))


# %% 3. Download the reference proteins from NCBI (only if not already saved)
def base_accession(acc):
    """Accession without the version number, e.g. WP_013582405.1 -> WP_013582405."""
    return acc.split(".")[0]


if REFERENCES_RAW.exists():
    print(f"{REFERENCES_RAW.name} already exists - not downloading again")
else:
    handle = Entrez.efetch(db="protein", id=",".join(refs["accession"]),
                           rettype="fasta", retmode="text")
    REFERENCES_RAW.write_text(handle.read())
    handle.close()
    print("Downloaded", len(refs), "sequences to", REFERENCES_RAW.name)

downloaded = {base_accession(r.id): r for r in SeqIO.parse(REFERENCES_RAW, "fasta")}

missing = [a for a in refs["accession"] if base_accession(a) not in downloaded]
if missing:
    raise SystemExit(f"Not downloaded - check these accessions: {missing}")


# %% 4. Label every sequence for the tree
def clean(text):
    """Make text safe for alignment and tree software: letters, numbers and _ only."""
    text = re.sub(r"[^A-Za-z0-9.]+", "_", text).strip("_")
    return text.replace(".", "")


records = []
label_rows = []

# references and outgroup: Genus_species_accession, e.g. Staphylococcus_aureus_YP_499096
for _, row in refs.iterrows():
    record = downloaded[base_accession(row["accession"])]
    species = " ".join(row["species"].split()[:2])          # Genus species only
    label = f"{clean(species)}_{base_accession(row['accession'])}"
    if row["group"] == "outgroup":
        label = "OUTGROUP_" + label
    record.id, record.name, record.description = label, label, ""
    records.append(record)
    label_rows.append({"label": label, "type": "outgroup" if row["group"] == "outgroup" else "reference",
                       "group": row["group"], "accession": row["accession"],
                       "species": row["species"], "role": row["role"]})

# representatives: GROUP1_outlet02_seq053 (capitals so they stand out in the tree)
for record in SeqIO.parse(REPRESENTATIVES, "fasta"):
    match = re.match(r"group(\d+)_outlet_(\d+)_(seq\d+)", record.id)
    label = f"GROUP{match.group(1)}_outlet{match.group(2)}_{match.group(3)}"
    original_id = record.id
    record.id, record.name, record.description = label, label, ""
    records.append(record)
    label_rows.append({"label": label, "type": "representative",
                       "group": f"group{match.group(1)}", "accession": original_id,
                       "species": "", "role": "representative"})


# %% 5. Save the combined FASTA file and the label table
for record in records:
    record.seq = record.seq.rstrip("*")   # remove any stop symbol

SeqIO.write(records, TREE_INPUT, "fasta")
labels = pd.DataFrame(label_rows)
labels.to_csv(LABELS, index=False)

print(f"\nSaved {len(records)} sequences to {TREE_INPUT.name}:")
for record in records:
    print(f"  {record.id:<45} {len(record.seq):>5} aa")

# %%