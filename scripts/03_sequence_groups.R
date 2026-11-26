# =============================================================================
# 03_sequence_groups.R
# Bristol Bacteria Meat Mystery
#
# Step 2: identify sequence groups
#   - pairwise JC69 distances
#   - PCoA (Figure 1: coloured by outlet)
#   - k-means clustering + silhouette analysis to choose k
#   - PCoA coloured by group (Figure 2)
#   - table of groups vs outlets
# Step 3 (start): pick one representative sequence per group (the medoid)
#
# Open Bristol-Bacteria-Meat-Mystery.Rproj first, then run this script.
# Input : results/all_sequences_aligned.fasta (Clustal Omega, scripts/02_align_clustalo.sh)
# Output: figures/ and results/
# =============================================================================


# ---- 1. Packages ------------------------------------------------------------

for (pkg in c("ape", "cluster", "ggplot2", "here")) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# here() builds paths from the project folder, so the script works whether
# you run it line by line or with File > Compile Report
dir.create(here("figures"), showWarnings = FALSE)
dir.create(here("results"), showWarnings = FALSE)


# ---- 2. Read the aligned sequences ------------------------------------------
# The raw sequences differ in length (3,126-3,615 bp), so they were aligned
# with Clustal Omega on the HPC before this script (see the README).
# Input: results/all_sequences_aligned.fasta

seqs <- read.FASTA(here("results", "all_sequences_aligned.fasta"))

# keep only the first word of each header, e.g. "outlet_03_seq014"
names(seqs) <- sub(" .*", "", names(seqs))

length(seqs)            # should be 100
table(lengths(seqs))    # should be a single length (the alignment length)

if (length(unique(lengths(seqs))) > 1) {
  stop("Sequences are not all the same length - check the alignment file.")
}

# outlet label taken from the header added in step 1
outlet <- sub("^(outlet_[0-9]+)_.*", "\\1", names(seqs))
table(outlet)


# ---- 3. Pairwise distances --------------------------------------------------

dna <- as.matrix(seqs)

# JC69 corrects raw p-distances for multiple substitutions at the same site
jc_dist <- dist.dna(dna, model = "JC69", pairwise.deletion = TRUE)

summary(as.vector(jc_dist))
# JC69 cannot be calculated if two sequences differ at >= 75% of sites
if (any(!is.finite(jc_dist))) {
  warning("Some JC69 distances are NaN/Inf - check for very divergent sequences.")
}


# ---- 4. Principal coordinates analysis (PCoA) -------------------------------

pcoa_result <- cmdscale(jc_dist, k = 10, eig = TRUE)

pcoa_df <- as.data.frame(pcoa_result$points)
colnames(pcoa_df) <- paste0("PCoA", seq_len(ncol(pcoa_df)))
pcoa_df$id     <- names(seqs)
pcoa_df$outlet <- outlet

# % of variation explained by each axis (from the positive eigenvalues)
eig <- pcoa_result$eig
variance <- 100 * eig / sum(eig[eig > 0])
round(variance[1:5], 1)

axis_x <- paste0("PCoA1 (", round(variance[1], 1), "%)")
axis_y <- paste0("PCoA2 (", round(variance[2], 1), "%)")

# Figure 1: PCoA coloured by food outlet
fig_outlet <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, colour = outlet)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(x = axis_x, y = axis_y, colour = "Food outlet") +
  theme_classic()

fig_outlet
ggsave(here("figures", "fig1_pcoa_by_outlet.png"), fig_outlet,
       width = 7, height = 5, dpi = 300)


# ---- 5. k-means + silhouette analysis to choose k ---------------------------

pcoa_coords <- pcoa_df[, grep("^PCoA", colnames(pcoa_df))]
pcoa_dist   <- dist(pcoa_coords)

# k-means can't make more groups than there are distinct points
# (identical sequences sit on top of each other)
n_distinct <- nrow(unique(round(pcoa_coords, 10)))
k_values   <- 2:min(10, n_distinct - 1)

set.seed(123)   # k-means starts from random centres
mean_silhouette <- sapply(k_values, function(k) {
  km <- kmeans(pcoa_coords, centers = k, nstart = 25)
  mean(silhouette(km$cluster, pcoa_dist)[, "sil_width"])
})

sil_table <- data.frame(k = k_values, mean_silhouette = round(mean_silhouette, 3))
sil_table
write.csv(sil_table, here("results", "silhouette_widths.csv"), row.names = FALSE)

fig_sil <- ggplot(sil_table, aes(x = k, y = mean_silhouette)) +
  geom_line() +
  geom_point(size = 3) +
  scale_x_continuous(breaks = k_values) +
  labs(x = "Number of clusters (k)", y = "Mean silhouette width") +
  theme_classic()

fig_sil
ggsave(here("figures", "silhouette_widths.png"), fig_sil,
       width = 6, height = 4, dpi = 300)

# k with the highest mean silhouette width
# (check this agrees with the clusters you can see in Figure 1)
optimal_k <- k_values[which.max(mean_silhouette)]
optimal_k


# ---- 6. Final clustering ----------------------------------------------------

set.seed(123)
km_final <- kmeans(pcoa_coords, centers = optimal_k, nstart = 25)
pcoa_df$group <- factor(km_final$cluster)   # group numbers are arbitrary labels

# Figure 2: PCoA coloured by inferred sequence group
fig_group <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, colour = group)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(x = axis_x, y = axis_y, colour = "Sequence group") +
  theme_classic()

fig_group
ggsave(here("figures", "fig2_pcoa_by_group.png"), fig_group,
       width = 7, height = 5, dpi = 300)

# Which groups are found at which outlets?
group_by_outlet <- table(Outlet = pcoa_df$outlet, Group = pcoa_df$group)
group_by_outlet
write.csv(as.data.frame.matrix(group_by_outlet),
          here("results", "groups_by_outlet.csv"))

write.csv(pcoa_df[, c("id", "outlet", "group")],
          here("results", "sequence_groups.csv"), row.names = FALSE)


# ---- 7. Representative sequence for each group (step 3) ---------------------
# The medoid: the sequence with the smallest average distance to all other
# sequences in its group, i.e. the most "typical" real sequence in the group.

dmat <- as.matrix(jc_dist)

representatives <- sapply(levels(pcoa_df$group), function(g) {
  members <- pcoa_df$id[pcoa_df$group == g]
  within  <- dmat[members, members, drop = FALSE]
  members[which.min(rowMeans(within))]
})
representatives

# Save them (gaps removed) for ORF prediction in Python
rep_seqs <- del.gaps(seqs[representatives])
names(rep_seqs) <- paste0("group", names(representatives), "_", representatives)
write.FASTA(rep_seqs, here("results", "representatives.fasta"))

# Summary of each sequence: group, ungapped length, ambiguous bases, centrality
ungapped <- del.gaps(seqs)
seq_summary <- data.frame(
  id          = names(seqs),
  group       = pcoa_df$group,
  length_bp   = lengths(ungapped),
  n_ambiguous = sapply(as.character(ungapped), function(s) sum(!s %in% c("a","c","g","t"))),
  mean_dist_to_group = sapply(seq_along(seqs), function(i) {
    members <- pcoa_df$id[pcoa_df$group == pcoa_df$group[i]]
    mean(dmat[names(seqs)[i], members])
  })
)


install.packages("dplyr")
library(dplyr)

# For each group: the 5 most central sequences, with their lengths
seq_summary |>
  dplyr::arrange(group, mean_dist_to_group) |>
  dplyr::group_by(group) |>
  dplyr::slice_head(n = 5)

# Range of lengths in each group, for comparison
tapply(seq_summary$length_bp, seq_summary$group, range)

write.csv(seq_summary, here("results", "representative_selection.csv"), row.names = FALSE)

# ---- 8. Record package versions for the README ------------------------------

writeLines(capture.output(sessionInfo()),
           here("results", "sessionInfo_R.txt"))
