# =============================================================================
# 05_blastp.py
# Bristol Bacteria Meat Mystery - Step 5: summarise BLASTP results
#
# The BLASTP search is run on the NCBI BLAST website (see README), with all
# predicted proteins in one search, and the results downloaded with
# "Download All -> XML". This script reads that file and summarises it to help:
#   - identify the protein encoded by the sequences
#   - identify the closest bacterial taxa for each sequence group
#   - choose reference sequences (one best hit per named species)
#
# Input : results/blast/all_queries_blastp.xml  (downloaded from NCBI BLAST)
# Output: results/blast/blast_hits.csv           every hit, all queries
#         results/blast/blast_summary.csv        top hit per query
#         results/blast/reference_candidates.csv best hit per named species
#         results/blast/blast_settings.txt       details of the search
#
# Run from the repository root:  python scripts/05_blastp.py
# or cell by cell in VS Code ("Run Cell").
# =============================================================================

# %% 1. Packages and file paths
import re
from pathlib import Path

import pandas as pd
from Bio.Blast import NCBIXML

ROOT = Path.cwd()
while not (ROOT / "results").exists() and ROOT != ROOT.parent:
    ROOT = ROOT.parent
BLAST_DIR = ROOT / "results" / "blast"
BLAST_XML = BLAST_DIR / "all_queries_blastp.xml"
print("Repository root:", ROOT)


# %% 2. Details of the web BLAST search (recorded for the README)
SEARCH_DETAILS = {
    "Source": "NCBI BLAST website, one search containing all predicted proteins",
    "RID": "GNHDXHWY014",
    "Date searched": "2026-11-28",
    "Program": "blastp",
    "Database": "RefSeq Select proteins (refseq_select_prot)",
    "Max target sequences": 100,
    "Expect threshold": 0.05,
    "Word size": 5,
    "Matrix": "BLOSUM62",
    "Gap costs": "existence 11, extension 1",
    "Compositional adjustments": "conditional compositional score matrix adjustment",
    "Low-complexity filter": "off",
}

MIN_QUERY_COVER = 90   # % query cover required for a reference candidate

settings_text = "\n".join(f"{k}: {v}" for k, v in SEARCH_DETAILS.items())
settings_text += f"\nReference candidates: query cover >= {MIN_QUERY_COVER}%\n"
(BLAST_DIR / "blast_settings.txt").write_text(settings_text)
print(settings_text)


# %% 3. Read the BLAST results into one table
def species_from(description):
    """Species name in the square brackets at the end of a description."""
    match = re.search(r"\[([^\[\]]+)\]\s*$", description)
    return match.group(1) if match else ""


def protein_from(description):
    """Protein name: the description without the [species] part."""
    name = re.sub(r"\s*\[[^\[\]]+\]\s*$", "", description)
    return name.replace("MULTISPECIES: ", "")


rows = []
with open(BLAST_XML) as fh:
    for blast_record in NCBIXML.parse(fh):
        # query name = first word of the FASTA header pasted into BLAST
        query_name = blast_record.query.split()[0]
        query_length = blast_record.query_length

        for rank, alignment in enumerate(blast_record.alignments, start=1):
            hsp = alignment.hsps[0]                          # best-scoring local alignment
            description = alignment.hit_def.split(" >")[0]   # first title if several
            rows.append({
                "query": query_name,
                "rank": rank,
                "accession": alignment.accession,
                "protein": protein_from(description),
                "species": species_from(description),
                "bit_score": round(hsp.bits, 1),
                "e_value": hsp.expect,
                "percent_identity": round(100 * hsp.identities / hsp.align_length, 2),
                "query_cover": round(100 * (hsp.query_end - hsp.query_start + 1) / query_length, 1),
                "hit_length": alignment.length,
            })

hits = pd.DataFrame(rows)
hits.to_csv(BLAST_DIR / "blast_hits.csv", index=False)
print("Queries:", hits["query"].nunique(), "| Total hits:", len(hits))


# %% 4. Summary: what protein is it, and what is the closest species?
summary_rows = []
for query, q_hits in hits.groupby("query", sort=False):
    top = q_hits.iloc[0]
    top20 = q_hits.head(20)
    summary_rows.append({
        "query": query,
        "most_common_protein_top20": top20["protein"].mode().iloc[0],
        "top_hit_species": top["species"],
        "top_hit_accession": top["accession"],
        "top_hit_identity": top["percent_identity"],
        "top_hit_query_cover": top["query_cover"],
        "top_hit_e_value": top["e_value"],
    })

summary = pd.DataFrame(summary_rows)
summary.to_csv(BLAST_DIR / "blast_summary.csv", index=False)
print(summary.to_string(index=False))


# %% 5. Reference candidates: best hit per named species
# Remove hits that would make poor references:
#   - unnamed or uncultured organisms (no clear species label for the tree)
#   - hypothetical or partial proteins (poorly annotated or incomplete)
#   - low query cover (would align over only part of the protein)
bad_words = r"uncultured|unclassified|hypothetical|partial|metagenome"
named_species = hits["species"].str.split().str.len() >= 2   # "Genus species"
not_sp = ~hits["species"].str.contains(r"\bsp\.", regex=True)

candidates = hits[
    named_species
    & not_sp
    & ~hits["species"].str.contains(bad_words, case=False, regex=True)
    & ~hits["protein"].str.contains(bad_words, case=False, regex=True)
    & (hits["query_cover"] >= MIN_QUERY_COVER)
].copy()

candidates["genus"] = candidates["species"].str.split().str[0]

# hits are already ordered best first, so keep the first hit for each species
candidates = candidates.drop_duplicates(subset=["query", "species"], keep="first")
candidates = candidates.groupby("query", sort=False).head(15)

candidates.to_csv(BLAST_DIR / "reference_candidates.csv", index=False)

for query, q_cand in candidates.groupby("query", sort=False):
    print(f"\n=== {query} ===")
    print(q_cand[["accession", "species", "genus", "percent_identity",
                  "query_cover", "e_value"]].to_string(index=False))

# %%