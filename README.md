______      _     _        _  ___  ___           _    ___  ___          _
| ___ \    (_)   | |      | | |  \/  |          | |   |  \/  |         | |
| |_/ /_ __ _ ___| |_ ___ | | | .  . | ___  __ _| |_  | .  . |_   _ ___| |_ ___ _ __ _   _
| ___ \ '__| / __| __/ _ \| | | |\/| |/ _ \/ _` | __| | |\/| | | | / __| __/ _ \ '__| | | |
| |_/ / |  | \__ \ || (_) | | | |  | |  __/ (_| | |_  | |  | | |_| \__ \ ||  __/ |  | |_| |
\____/|_|  |_|___/\__\___/|_| \_|  |_/\___|\__,_|\__| \_|  |_/\__, |___/\__\___|_|   \__, |
                                                               __/ |                  __/ |
___________________________________________________________________________________________
-------------------------------------------------------------------------------------------

Bioinformatic investigation of 100 bacterial DNA sequences isolated from meat
samples at 8 food outlets (simulated outbreak exercise). This repository
contains all code needed to reproduce the analysis in the accompanying report.






**> Items in [square brackets] are placeholders - fill them in from your own results.
**
---




## Repository structure

```
data/                          raw FASTA files, one folder per outlet (never modified)
  outlet_01/ ... outlet_08/
scripts/
  01_prepare_data.sh           Step 1: combine and label sequences (Bash)
  02_align_clustalo.sh         Step 2a: Slurm job - align sequences with Clustal Omega (Blue Pebble HPC)
  03_sequence_groups.R         Steps 2b-3: distances, PCoA, clustering, representatives (R)
combined/                      copies of each FASTA file with the outlet added to the header
all_sequences.fasta            all 100 labelled sequences in one multi-FASTA file (unaligned)
results/
  all_sequences_aligned.fasta  Clustal Omega alignment of all 100 sequences
  slurm-20382747.out           log from the Blue Pebble alignment job
  silhouette_widths.csv        mean silhouette width for each k
  groups_by_outlet.csv         number of sequences from each group at each outlet
  sequence_groups.csv          group assigned to every sequence
  representatives.fasta        one representative sequence per group (ungapped)
  sessionInfo_R.txt            R and package versions
figures/
  fig1_pcoa_by_outlet.png      PCoA coloured by food outlet
  fig2_pcoa_by_group.png       PCoA coloured by inferred sequence group
  silhouette_widths.png        mean silhouette width against k
```

---

## Software and dependencies

| Step | Software | Version | Where it was run |
|---|---|---|---|
| 1 | Bash, `sed`, `cat` | [bash --version] | GitHub Codespaces (Ubuntu Linux) |
| 2a | Git Bash (`ssh`, `scp`) | [ ] | Local Windows PC |
| 2a | Slurm job scheduler | - | Blue Pebble HPC, University of Bristol |
| 2a | Clustal Omega | 1.2.4 (module `clustal-omega/1.2.4-skoy`) | Blue Pebble HPC |
| 2b-3 | R | [R.version.string] | RStudio on a local Windows PC |
| 2b-3 | ape | [see results/sessionInfo_R.txt] | |
| 2b-3 | cluster | [ ] | |
| 2b-3 | ggplot2 | [ ] | |
| 2b-3 | here | [ ] | |

Exact R package versions are saved automatically to `results/sessionInfo_R.txt`
when `03_sequence_groups.R` is run.

---

## How to reproduce the analysis

1. Clone the repository:
   `git clone https://github.com/dylanphillips076-sketch/Bristol-Bacteria-Meat-Mystery.git`
2. **Step 1** - from the repository root in a Bash terminal:
   `bash scripts/01_prepare_data.sh`
   Check: `grep -c ">" all_sequences.fasta` should return `100`.
3. **Step 2a** - align the sequences on Blue Pebble (manual commands below,
   see *Step 2a*). The resulting alignment is already included as
   `results/all_sequences_aligned.fasta`, so this step can be skipped to
   reproduce the downstream results.
4. **Steps 2b-3** - open `Bristol-Bacteria-Meat-Mystery.Rproj` in RStudio and run
   `scripts/03_sequence_groups.R` from top to bottom (missing packages are
   installed automatically).

---

## Step 1: Data preparation (Bash)

**Script:** `scripts/01_prepare_data.sh`

**What it does**

1. Creates a folder `combined/`.
2. Loops over each outlet folder in `data/` (`outlet_01` to `outlet_08`).
3. For every FASTA file in that folder, uses `sed "s/^>/>${outlet}_/"` to add the
   outlet name to the start of the sequence header
   (e.g. `>seq014` becomes `>outlet_03_seq014`), and writes the edited copy to
   `combined/<outlet>_<filename>`.
4. Concatenates all labelled files into a single multi-FASTA file,
   `all_sequences.fasta` (100 sequences).

**Rationale**

- In the raw data the outlet is recorded only by the folder a file sits in. Once
  sequences are loaded into other software that information would be lost, so the
  outlet is written into each header and stays with the sequence through every
  later step (outlet colours in Figure 1, tip labels in the tree).
- The raw files in `data/` are never edited; all changes are written to new files.
  This keeps the original data intact and means the step can be rerun from scratch.
- The outlet is also added to each file name so files with the same name in
  different outlets cannot overwrite each other.
- A single multi-FASTA file is the standard input for alignment and BLAST.

**Manual step**

- The eight `outlet_` folders were uploaded into `data/` by hand (drag and drop into
  the Codespace), keeping the folder structure so the outlet information was
  preserved. A macOS system file (`.DS_Store`) in the original download was not
  included.

---

## Step 2a: Multiple sequence alignment (Clustal Omega on Blue Pebble)

**Script:** `scripts/02_align_clustalo.sh` (Slurm job script)

**Why an alignment was needed**

The sequences varied in length (3,126-3,615 bp), so positions could not be compared
directly. Pairwise distances require homologous positions to be lined up in the
same column, so all 100 sequences were aligned first.

**What the job does**

- Requests 1 node, 1 task, 1000 MB memory and up to 1 hour on the `teach_cpu`
  partition, under the course HPC account.
- Loads Clustal Omega with `module load clustal-omega/1.2.4-skoy`.
- Runs:
  `clustalo -i ~/job_data/all_sequences.fasta -o ~/job_data/all_sequences_aligned.fasta --outfmt=fasta --force`
  with all other settings at their defaults.
- Prints job details (Clustal Omega version, host, job ID, node, end time) to the
  Slurm log, `results/slurm-20382747.out`.

**Rationale**

- **Clustal Omega** is a widely used, fast and accurate multiple sequence aligner.
  It builds a guide tree from k-mer similarities and then aligns sequences
  progressively, which scales well to around 100 sequences of a few kb.
- **Blue Pebble HPC** was used so that the alignment ran as a recorded, reproducible
  batch job (job script plus log) through the Slurm scheduler, rather than
  interactively, and so that the same workflow would scale to larger datasets.
- **Aligned FASTA output** (`--outfmt=fasta`) can be read directly into R.

**Manual steps (run in Git Bash on a local PC on the university network/VPN)**

Replace `<username>` with your University of Bristol username.

```bash
# from the repository root: stage the input files
mkdir job_data
cp all_sequences.fasta scripts/02_align_clustalo.sh job_data/
scp -r job_data <username>@bp1-login.acrc.bris.ac.uk:

# log in to Blue Pebble and submit the job
ssh -X <username>@bp1-login.acrc.bris.ac.uk
module spider clustal                         # confirm module name
sbatch job_data/02_align_clustalo.sh
squeue --me                                   # job leaves the queue when finished
grep -c ">" job_data/all_sequences_aligned.fasta   # check: 100
exit

# back on the local PC: retrieve the results
mkdir -p results
scp <username>@bp1-login.acrc.bris.ac.uk:job_data/all_sequences_aligned.fasta results/
scp <username>@bp1-login.acrc.bris.ac.uk:slurm-<jobID>.out results/
rm -r job_data
```

**Notes**

- The account code in the course worksheet was not valid for this unit. The correct
  account was found with
  `sacctmgr show associations user=<username> format=account,partition`
  and set in the job script (`#SBATCH --account=bisc040764`). Other users must
  replace this with their own account code.
- `zip` is not available in Git Bash, so the folder was transferred with
  `scp -r` instead of zipping it first.
- The job (ID 20382747) ran on node `bp1-compute189` and finished on 24 Nov 2026.
  The output was checked to contain all 100 sequences.

---

## Step 2b: Identifying sequence groups (R)

**Script:** `scripts/03_sequence_groups.R` (sections 1-6)

### Reading the alignment

- `results/all_sequences_aligned.fasta` is read with `ape::read.FASTA()`.
- Headers are trimmed to their first word, and the outlet is extracted from each
  header with a regular expression (`^(outlet_[0-9]+)_`).
- The script checks that there are 100 sequences and that all have the same
  (alignment) length ([N] columns).

### Pairwise distances

- **Function:** `ape::dist.dna(model = "JC69", pairwise.deletion = TRUE)`
- **Why JC69:** raw (p) distances underestimate divergence because a site can
  change more than once, or change back, over time. The Jukes-Cantor (JC69) model
  corrects for these hidden multiple substitutions. It is the simplest correction
  and is appropriate for grouping sequences; more complex models are used later
  for the phylogeny.
- **Why pairwise deletion:** the sequences start and end at different positions, so
  the alignment contains gaps, especially at the ends. Pairwise deletion ignores a
  site only for the pair of sequences where one has a gap, so gaps are not counted
  as differences and as much data as possible is used for each comparison.

### Principal coordinates analysis (PCoA)

- **Function:** `cmdscale(jc_dist, k = 10, eig = TRUE)`
- **Why PCoA:** the data are a distance matrix rather than a table of variables, so
  PCoA (classical multidimensional scaling) is used to place each sequence in a
  low-dimensional space while preserving the pairwise distances as closely as
  possible. This makes clusters of similar sequences visible.
- The percentage of variation explained by each axis is calculated from the
  positive eigenvalues. PCoA1 and PCoA2 explain [x]% and [y]% respectively.
- **Figure 1** (`figures/fig1_pcoa_by_outlet.png`): PCoA1 vs PCoA2, points coloured
  by food outlet.

### k-means clustering and silhouette analysis

- **Clustering:** `kmeans(centers = k, nstart = 25)` on the first 10 PCoA axes.
  - k-means needs coordinates in Euclidean space, which the PCoA axes provide.
  - 10 axes are used rather than 2 so that clustering uses more of the variation in
    the data than is visible in the 2D plot.
  - `nstart = 25` runs each clustering from 25 random starting points and keeps the
    best, reducing the chance of a poor solution.
  - `set.seed(123)` makes the random starts, and therefore the results, reproducible.
- **Choosing k:** for k = 2 to 10, the mean silhouette width is calculated with
  `cluster::silhouette()` using Euclidean distances in the PCoA space. Silhouette
  width measures how well each sequence fits its own cluster compared with the
  nearest other cluster (from -1 to 1; higher means better separated). The range of
  k is capped below the number of distinct sequences, because k-means cannot form
  more clusters than there are distinct points.
- **Result:** k = [k] had the highest mean silhouette width ([value]; full table in
  `results/silhouette_widths.csv`, plot in `figures/silhouette_widths.png`).
  [Explain whether this agrees with the clusters visible in Figure 1.]
- **Figure 2** (`figures/fig2_pcoa_by_group.png`): PCoA1 vs PCoA2, points coloured
  by inferred sequence group.
- The distribution of groups across outlets is saved in
  `results/groups_by_outlet.csv`, and each sequence's group in
  `results/sequence_groups.csv`. Group numbers assigned by k-means are arbitrary
  labels.

---

## Step 3: Representative sequences (R)

**Script:** `scripts/03_sequence_groups.R` (section 7)

- For each group, the representative is the **medoid**: the sequence with the
  smallest mean JC69 distance to all other sequences in the same group.
- **Rationale:** the medoid is a real, observed sequence (unlike an average or
  consensus) and is the most typical member of its group, so it is unlikely to be
  an outlier or to carry unusual sequencing errors.
- Selected representatives: [group 1: id], [group 2: id], [...]
- Alignment gaps are removed (`ape::del.gaps()`) and the sequences are saved, with
  the group added to the header, to `results/representatives.fasta` for ORF
  prediction.

---

## Limitations so far

- Clustal Omega does not detect sequences submitted in reverse-complement
  orientation; any such sequence would align poorly and appear as an outlier in
  the PCoA.
- k-means assumes clusters that are roughly spherical and similar in size, which may
  not hold for genetic data.
- Silhouette width is a guide to the best-supported number of groups, not proof that
  the groups are biologically distinct; the choice of k was therefore checked against
  the PCoA plot and later confirmed with BLAST and phylogenetic analysis.
- JC69 assumes equal base frequencies and equal substitution rates between all bases.

---

## Steps 4-7

[To be added: ORF prediction and translation (Python), BLASTP and reference
selection, protein alignment, phylogeny.]
