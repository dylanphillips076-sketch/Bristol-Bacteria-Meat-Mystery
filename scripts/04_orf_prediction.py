# =============================================================================
# 04_orf_prediction.py
# Bristol Bacteria Meat Mystery - Step 4: ORF prediction and translation
#
# Adapted and extended from the Biopython practical (Exercise 8.3, longest ORF).
# Extensions to the practical code:
#   - searches all SIX reading frames (3 forward + 3 on the reverse strand),
#     not just the 3 forward frames
#   - an ORF must begin with a bacterial start codon (ATG, GTG or TTG), not
#     just be a stretch without stop codons
#   - uses the bacterial genetic code (NCBI table 11)
#   - records the strand, frame and nucleotide coordinates of each ORF and
#     whether it ends in a stop codon
#   - saves the candidate ORFs so the choice of ORF can be justified
#
# Input : results/representatives.fasta (made by scripts/03_sequence_groups.R)
# Output: results/orf_candidates.csv
#         results/representative_orfs_dna.fasta
#         results/representative_proteins.fasta   <- use this for BLASTP
#         results/protein_stats.csv
#
# Written with "# %%" cells: in VS Code click "Run Cell" above each cell,
# or run the whole file from the repository root with:
#   python scripts/04_orf_prediction.py
# =============================================================================

# %% 1. Packages and file paths
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# find the repository root (the folder containing "results"), so the script
# works whether it is run from the root or cell by cell from scripts/
ROOT = Path.cwd()
while not (ROOT / "results").exists() and ROOT != ROOT.parent:
    ROOT = ROOT.parent
RESULTS = ROOT / "results"
print("Repository root:", ROOT)


# %% 2. Settings
GENETIC_CODE = 11   # NCBI table 11: bacterial, archaeal and plant plastid code
# start codons used by bacteria (table 11). ATG is the most common, but GTG and
# TTG are also used; all are translated as methionine (M) when used as a start
START_CODONS = ("ATG", "GTG", "TTG")
MIN_ORF_AA = 50     # ignore very short ORFs (in amino acids) when listing candidates
N_CANDIDATES = 5    # number of candidate ORFs saved per sequence


# %% 3. Function to find ORFs in all six reading frames
def find_orfs(dna_seq, table=GENETIC_CODE, min_aa=MIN_ORF_AA):
    """Return all ORFs (start codon ATG -> stop codon) in the six reading
    frames of dna_seq, longest first."""
    orfs = []
    seq_len = len(dna_seq)

    # +1 = the sequence as given, -1 = its reverse complement
    for strand, nuc in [(+1, dna_seq), (-1, dna_seq.reverse_complement())]:
        for frame in range(3):
            # trim the end so the length is a multiple of 3
            end = frame + ((seq_len - frame) // 3) * 3
            protein = str(nuc[frame:end].translate(table=table, to_stop=False))

            # as in the practical: split the translation at stop codons ("*")
            aa_pos = 0   # position of the current fragment in the translation
            for fragment in protein.split("*"):
                # extension: the ORF starts at the first start codon in the fragment
                m = -1
                for i in range(len(fragment)):
                    codon_start = frame + 3 * (aa_pos + i)
                    if str(nuc[codon_start:codon_start + 3]) in START_CODONS:
                        m = i
                        break
                if m != -1 and len(fragment) - m >= min_aa:
                    start_codon = str(nuc[frame + 3 * (aa_pos + m):frame + 3 * (aa_pos + m) + 3])
                    # a start codon is always translated as methionine
                    orf_protein = "M" + fragment[m + 1:]
                    # the fragment ends in a stop codon unless it runs off the end
                    has_stop = aa_pos + len(fragment) < len(protein)

                    # nucleotide coordinates on this strand (0-based, end exclusive,
                    # including the stop codon if there is one)
                    nt_start = frame + 3 * (aa_pos + m)
                    nt_end = frame + 3 * (aa_pos + len(fragment)) + (3 if has_stop else 0)

                    # convert to 1-based coordinates on the original sequence
                    if strand == 1:
                        start, stop = nt_start + 1, nt_end
                    else:
                        start, stop = seq_len - nt_end + 1, seq_len - nt_start

                    orfs.append({
                        "strand": "+" if strand == 1 else "-",
                        "frame": frame + 1,
                        "start_codon": start_codon,
                        "start_nt": start,
                        "end_nt": stop,
                        "length_aa": len(orf_protein),
                        "length_nt": nt_end - nt_start,
                        "ends_with_stop": has_stop,
                        "protein": orf_protein,
                        "orf_dna": str(nuc[nt_start:nt_end]),
                    })
                aa_pos += len(fragment) + 1   # +1 for the "*"

    return sorted(orfs, key=lambda o: o["length_aa"], reverse=True)


# %% 4. Read the representative sequences
records = list(SeqIO.parse(RESULTS / "representatives.fasta", "fasta"))
print("Number of representative sequences:", len(records))

for record in records:
    # make sure the sequence is upper case DNA with no gaps
    record.seq = Seq(str(record.seq).upper().replace("-", ""))
    print(f"{record.id}: {len(record.seq)} bp")


# %% 5. Find the ORFs in each sequence
all_candidates = []
best_orfs = {}

for record in records:
    orfs = find_orfs(record.seq)
    if not orfs:
        print(f"{record.id}: no ORF of at least {MIN_ORF_AA} aa found")
        continue

    best_orfs[record.id] = orfs[0]

    # keep the top candidates for the table
    for rank, orf in enumerate(orfs[:N_CANDIDATES], start=1):
        row = {"sequence": record.id, "rank": rank,
               "seq_length_nt": len(record.seq)}
        row.update({k: v for k, v in orf.items() if k not in ("protein", "orf_dna")})
        all_candidates.append(row)

candidates = pd.DataFrame(all_candidates)
candidates.to_csv(RESULTS / "orf_candidates.csv", index=False)
print(candidates.to_string(index=False))


# %% 6. Check the chosen (longest) ORF for each sequence
summary = []
for record in records:
    if record.id not in best_orfs:
        continue
    best = best_orfs[record.id]
    others = candidates[(candidates["sequence"] == record.id) & (candidates["rank"] == 2)]
    second = int(others["length_aa"].iloc[0]) if len(others) else 0

    summary.append({
        "sequence": record.id,
        "strand": best["strand"],
        "frame": best["frame"],
        "start_nt": best["start_nt"],
        "end_nt": best["end_nt"],
        "length_aa": best["length_aa"],
        "second_longest_aa": second,
        "start_codon": best["start_codon"],
        "ends_with_stop": best["ends_with_stop"],
        "length_nt_multiple_of_3": best["length_nt"] % 3 == 0,
        "percent_of_sequence": round(100 * best["length_nt"] / len(record.seq), 1),
    })

orf_summary = pd.DataFrame(summary)
print(orf_summary.to_string(index=False))


# %% 7. Save the chosen ORFs (DNA) and their translations (protein)
dna_records = []
protein_records = []

for record in records:
    if record.id not in best_orfs:
        continue
    best = best_orfs[record.id]
    info = (f"longest ORF strand {best['strand']} frame {best['frame']} "
            f"start {best['start_codon']} "
            f"{best['start_nt']}-{best['end_nt']} "
            f"{best['length_aa']} aa genetic code {GENETIC_CODE}")

    dna_records.append(SeqRecord(Seq(best["orf_dna"]), id=record.id, description=info))
    protein_records.append(SeqRecord(Seq(best["protein"]), id=record.id, description=info))

SeqIO.write(dna_records, RESULTS / "representative_orfs_dna.fasta", "fasta")
SeqIO.write(protein_records, RESULTS / "representative_proteins.fasta", "fasta")
print("Saved", len(protein_records), "protein sequences to",
      RESULTS / "representative_proteins.fasta")


# %% 8. Protein statistics (as in practical Exercise 8.4)
ambiguous_aas = set("XBZJUO")
stats = []

for rec in protein_records:
    seq_str = str(rec.seq)
    ambiguous_count = sum(aa in ambiguous_aas for aa in seq_str)
    row = {"sequence": rec.id, "length_aa": len(seq_str),
           "ambiguous_count": ambiguous_count}

    # ProteinAnalysis cannot handle ambiguous residues
    if ambiguous_count == 0:
        prot = ProteinAnalysis(seq_str)
        row.update({
            "molecular_weight": round(prot.molecular_weight(), 1),
            "isoelectric_point": round(prot.isoelectric_point(), 2),
            "gravy": round(prot.gravy(), 3),
        })
    stats.append(row)

protein_stats = pd.DataFrame(stats)
protein_stats.to_csv(RESULTS / "protein_stats.csv", index=False)
print(protein_stats.to_string(index=False))

# %%