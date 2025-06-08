# Identify non-empty drops with EmptyDrops R package

setwd(rprojroot::find_rstudio_root_file())
renv::load()
suppressPackageStartupMessages({
  source(file.path("utils", "env.R"))
  library(argparse)
  library(assertthat)
  library(tidyverse)
  library(qs)
  library(DelayedMatrixStats)
  library(BiocParallel)
  library(DropletUtils)
  library(Seurat)
})

# Parameters

parser <- ArgumentParser(description = "Identify empty drops with EmptyDrops")
parser$add_argument("--counts_dir", type = "character", required = TRUE,
                    help = "10x outs directory")
parser$add_argument("--seed_val", type = "integer", default = 827,
                    help = "Seed value")
parser$add_argument("--n_cores", type = "integer", default = 1,
                    help = "Number of cores")
parser$add_argument("--emptydrops_alpha", type = "double", default = 0.001,
                    help = "Alpha parameter")
parser$add_argument("--out_dir", type = "character", required = TRUE,
                    help = "Output directory")
args <- parser$parse_args()
for (i in seq_along(args)) {assign(names(args)[i], args[[i]])}

#----- MAIN -----

pdf(file.path(out_dir, "plots.pdf"), height = 9, width = 16)
par(mfrow = c(1,1))

# Convert read counts data to SingleCellExperiment

if (!dir.exists(counts_dir)) {
  counts_dir = paste0(counts_dir, ".h5")
  # To bypass problem that DropletUtils::read10xCounts() loads matrix as DelayedMatrix instead of dgCMatrix
  x <- Seurat::Read10X_h5(counts_dir, use.names = FALSE, unique.features = FALSE)
  x <- do.call("rbind", x)
  sce <- DropletUtils::read10xCounts(c(raw = counts_dir), col.names = TRUE)
  if (identical(dimnames(x), dimnames(counts(sce)))) {
    counts(sce) <- x
  }
} else {
  sce <- DropletUtils::read10xCounts(c(raw = counts_dir), col.names = TRUE)
}
is_hto <- grepl("HTO", rowData(sce)$ID, ignore.case = FALSE)
is_antibody_capture <- grepl("Antibody Capture", rowData(sce)$Type)
rowData(sce)$Type[is_hto & is_antibody_capture] <- "Multiplexing Capture"

# PLOT: Experiments present and show number of features per experiment
plot_df <- enframe(table(rowData(sce)$Type), name = "experiment", value = "n_features") %>% 
  mutate(n_features = as.numeric(n_features))
ggplot(plot_df, aes(experiment, n_features)) +
  geom_col() +
  geom_text(aes(label = n_features), vjust = -0.5)

# Split experiments
sce <- SingleCellExperiment::splitAltExps(sce, rowData(sce)$Type)

## Checks
checks <- vapply(altExpNames(sce), USE.NAMES = FALSE, FUN.VALUE = logical(1), FUN = function(x) {
  identical(colnames(sce), colnames(altExp(sce, x)))
})
assert_that(all(checks),
            msg = "Experiments have different droplet barcode names after splitting")

assert_that(mainExpName(sce) == "Gene Expression",
            msg = "mainExpName(sce) not Gene Expression")

# Remove barcodes with 0 counts
bctotalcount <- DelayedMatrixStats::colSums2(
  DelayedArray(assay(sce, "counts"))
  #assay(sce, "counts")
)
# This applies to all experiments
sce_nonzero <- sce[, bctotalcount > 0]

# DropletUtils::emptyDrops()

set.seed(seed_val)
emptydrops_out <- DropletUtils::emptyDrops(assay(sce_nonzero, "counts"),
                                           lower = 100,
                                           # "Setting test.ambient=TRUE will also modify the p-values prior to correction"
                                           # Set to true in a separate run of function to check the distribution of p-values for low-total barcodes
                                           test.ambient = FALSE,
                                           BPPARAM = MulticoreParam(n_cores, RNGseed = seed_val)
                                           )

# Check
# "emptyDrops() assumes that barcodes with low total UMI counts are empty droplets. 
# Thus, the null hypothesis should be true for all of these barcodes. 
# We can check whether the hypothesis testing procedure holds its size by examining the distribution of p-values for low-total barcodes with test.ambient=TRUE.
# Ideally, the distribution should be close to uniform. 
# Large peaks near zero indicate that barcodes with total counts below lower are not all ambient in origin. 
# This can be resolved by decreasing lower further to ensure that barcodes corresponding to droplets with very small cells are not used to estimate the ambient profile." (OSCA book)
set.seed(seed_val)
emptydrops_all_out <- DropletUtils::emptyDrops(assay(sce_nonzero, "counts"),
                                               lower = 100,
                                               test.ambient = TRUE,
                                               BPPARAM = MulticoreParam(n_cores, RNGseed = seed_val))
lower_limit = 100
hist(emptydrops_all_out$PValue[emptydrops_all_out$Total <= lower_limit & 
                                 emptydrops_all_out$Total > 0],
     xlab = "P-value",col = "grey80", cex.main = 0.5,
     main = paste0("DropletUtils::emptyDrops(test.ambient = TRUE), lower = ", lower_limit,
                   "\n'lower' value used to determine ambient cells needs to be adjusted?",
                   "\nIf not uniform and larger peaks in 0, some drops assumed to be ambient may not be and so decrease 'lower' further to ensure that barcodes corresponding to droplets with very small cells are not used to estimate the ambient profile")
)
rm(emptydrops_all_out)

# PLOT: Remaining non-empty drops at different FDR thresholds and whether some droplets Limited and need to rerun emptyDrops()
# "If any non-significant barcodes are TRUE for Limited, we may need to increase the number of iterations" (OSCA book)
plot_df <- emptydrops_out[, c("FDR", "Limited")]
for (alpha_val in c(0.001, 0.01, 0.05)) {
  n_notemptydrops <- sum(plot_df$FDR < alpha_val, na.rm = TRUE)
  print(
    ggplot(plot_df, aes(FDR < alpha_val)) +
      geom_bar() +
      geom_text(stat = 'count', aes(label = after_stat(count), vjust = -0.5)) +
      labs(title = paste0(n_notemptydrops, " non-empty drops at FDR < ", alpha_val, 
                          "; facet is whether Limited or not"),
           x = paste0("FDR < ", alpha_val)) +
      facet_grid(~ Limited)
  )
}
# PLOT: distribution of FDRs to help decide what FDR alpha to use 
# e.g. there's a peak between 0.05 and 0.01 then maybe don't be too stringent or investigate further what these are
ggplot(plot_df, aes(-log10(FDR))) + 
  geom_histogram() +
  scale_x_continuous(breaks = -log10(c(0, 0.001, 0.01, 0.05, 1))) +
  labs(title = "Count of non-empty drops based on FDR alpha")

# Save outputs

# Whitelist of non-empty barcodes
whitelist_inds <- which(emptydrops_out$FDR < emptydrops_alpha)
writeLines(
  rownames(emptydrops_out)[whitelist_inds], file.path(out_dir, "whitelist.txt")
)
writeLines(
  rownames(emptydrops_out)[setdiff(1:nrow(emptydrops_out), whitelist_inds)],
  file.path(out_dir, "blacklist.txt")
)
# Ambient profile here can be used for measuring ambient contamination with DecontX
emptydrops_out$Barcode <- rownames(emptydrops_out)
qsave(emptydrops_out, file.path(out_dir, "output.qs"))

dev.off()

# Session info

sessioninfo::session_info(to_file = file.path(out_dir, "session_info.txt"))
sessioninfo::session_info()

# rm(list = ls()); gc()
