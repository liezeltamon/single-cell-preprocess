# Identify doublets with scDblFinder
#Rscript doublet_scdblfinder.R --counts_dir='data/cellranger/SLE_1A_total/raw_feature_bc_matrix' --out_dir=test/doublet --whitelist_path=results/process_droplets/empty/emptydrops_alpha0.001_seed827/SLE_1A_total/whitelist.txt --dehash_path="test/barcode_metadata.csv" &

setwd(rprojroot::find_rstudio_root_file())
renv::load()
suppressPackageStartupMessages({
  source(file.path("utils", "env.R"))
  source(file.path("utils", "sc_helpers.R"))
  library(argparse)
  library(cowplot)
  library(qs)
  library(tidyverse)
  library(BiocParallel)
  library(Seurat)
  library(scDblFinder)
})

## Parameters

parser <- ArgumentParser(description = "Identify doublets with scDblFinder")
parser$add_argument("--counts_dir", type = "character", required = TRUE,
                    help = "10x counts directory")
parser$add_argument("--out_dir", type = "character", required = TRUE,
                    help = "Output directory")
parser$add_argument("--whitelist_path", type = "character", required = TRUE,
                    help = "Path to empty droplet barcode whitelist")
parser$add_argument("--seed_val", type = "integer", default = 827,
                    help = "Random seed value (default: 827)")
parser$add_argument("--n_cores", type = "integer", default = 1,
                    help = "Number of CPU cores to use (default: 1)")
parser$add_argument("--doublet_expected_rate", type = "numeric", default = NULL,
                    help = "Expected doublet rate (default: NULL, 1% per 1000 cells by default in scDblFinder)")
parser$add_argument("--dehash_path", type = "character", default = NULL,
                    help = "Path to dehashing results CSV file (optional, for plotting)")
args <- parser$parse_args()
for (i in seq_along(args)) {assign(names(args)[i], args[[i]])}
out_dir <- create_dir(out_dir)

#seq_run_id = "sle_batch2_pool2"
#counts_dir = file.path("data", "cellranger", seq_run_id, "raw_feature_bc_matrix")
#seed_val = 827
#n_cores = 2
## 0.8% for every 1,000 cells (For the standard 3' v3.1 assay, may be different depending on kit)
## https://kb.10xgenomics.com/hc/en-us/articles/360054599512-What-is-the-cell-multiplet-rate-when-using-the-3-CellPlex-Kit-for-Cell-Multiplexing
#doublet_expected_rate = NULL # 1% per 1000 cells by default in scDblFinder()
#whitelist_path = file.path(
#  paste0("results/process_droplets/empty/emptydrops_alpha0.001_seed827"),
#  seq_run_id, "whitelist.txt"
#)

# (Optional) Plot output with demultiplex results (if available)
#dehash_method_id = "hasheddrops"
#dehash_path = file.path(
#  paste0("results/process_droplets/dehash/", dehash_method_id, "_emptydrops_alpha0.001_seed827"),
#  seq_run_id, "barcode_metadata.csv"
#)
#out_dir = create_dir(file.path("results", "process_droplets", "multiplet"))

#for (param in names(configs)) {
# eval(parse(text = paste0(param, " <- configs$", param)))
#}
#cat("Parameters: \n")
#print(configs)

#out_dir <- create_dir(file.path(out_dir, paste0("scdblfinder_seed", seed_val,
#                                                "_doublet_expected_rate", 
#                                                doublet_expected_rate), seq_run_id))

# Fixed for now
dehash_assignment_var = "assignment"
dehash_hashing_var = "hashing_var"
dehash_is_multiplet_var = "is_multiplet"
dehash_is_confident_var = "is_confident"

# ----- MAIN -----

pdf(file.path(out_dir, "plots.pdf"), height = 9, width = 16)
par(mfrow = c(1,1))

# Convert read counts data to SingleCellExperiment

sce <- read_sc_data(counts_dir)
is_hto <- grepl("HTO", rowData(sce)$ID, ignore.case = FALSE)
is_antibody_capture <- grepl("Antibody Capture", rowData(sce)$Type)
rowData(sce)$Type[is_hto & is_antibody_capture] <- "Multiplexing Capture"
# Split experiments
sce <- SingleCellExperiment::splitAltExps(sce, rowData(sce)$Type)

## Checks
checks <- vapply(altExpNames(sce), USE.NAMES = FALSE, FUN.VALUE = logical(1), FUN = function(x) {
  identical(colnames(sce), colnames(altExp(sce, x)))
})
if (!all(checks)) {
  stop("Experiments have different droplet barcode names after splitting")
}
if (mainExpName(sce) != "Gene Expression") {
  stop("mainExpName(sce) not Gene Expression")
}

# Prepare sce needed downstream

not_empty_barcodes <- readLines(whitelist_path)
sce_nonzero_notempty <- sce[,not_empty_barcodes]

# Identify doublets

if (!is.null(doublet_expected_rate) && doublet_expected_rate == "use_demuxafy") {
  doublet_expected_rate <- ncol(sce_nonzero_notempty) / 1000 * 0.008
}
message("doublet_expected_rate: ", doublet_expected_rate)

set.seed(seed_val)
scdblfinder_sce <- scDblFinder(counts(sce_nonzero_notempty), 
                               dbr = doublet_expected_rate,
                               BPPARAM = MulticoreParam(n_cores, RNGseed = seed_val))
scdblfinder_out_df <- colData(scdblfinder_sce) %>% 
  as.data.frame() %>% 
  rownames_to_column("Barcode")

# Outputs to save

qsave(scdblfinder_out_df, file.path(out_dir, "output.qs"))
# Abridged scdblfinder_out_df that can be added as cell metadata
scdblfinder_mdta_df <- scdblfinder_out_df %>% 
  mutate(assignment = scDblFinder.class,
         score = scDblFinder.score,
         is_multiplet = scDblFinder.class == "doublet") %>% 
  select(Barcode, assignment, score, is_multiplet)
write.csv(scdblfinder_mdta_df, file.path(out_dir, "barcode_metadata.csv"), row.names = FALSE)
scdblfinder_whitelist <- scdblfinder_mdta_df %>%
  filter(!is_multiplet) %>%
  pull(Barcode)
writeLines(scdblfinder_whitelist, file.path(out_dir, "whitelist.txt"))

if (!is.null(dehash_path)) {
  
  dehash_mdta_df <- read.csv(dehash_path)
  dehash_mdta_df <- dehash_mdta_df[ ,!duplicated(as.list(dehash_mdta_df))]
  dehash_mdta_df <- dehash_mdta_df %>% 
    dplyr::rename(dehash_assignment = dehash_assignment_var,
                  dehash_hashing_var = dehash_hashing_var,
                  dehash_is_multiplet = dehash_is_multiplet_var,
                  dehash_is_confident = dehash_is_confident_var)
  plot_df <- cbind(dehash_mdta_df, scdblfinder_mdta_df)
  plot_df <- plot_df[ ,!duplicated(as.list(plot_df))]
  
} else {
  plot_df <- scdblfinder_mdta_df
}

# PLOT: Doublets from demultiplexing (if done) and scdblfinder
multiplet_perc_combined <- round(sum(plot_df$dehash_is_multiplet | 
                                     plot_df$is_multiplet) /
                                   ncol(sce_nonzero_notempty) * 100, digits = 4)
p1 <- ggplot(plot_df, aes(is_multiplet, fill = dehash_is_multiplet)) +
  geom_bar(position = position_dodge(1)) +
  geom_text(stat = 'count', aes(label = after_stat(count)), position = position_dodge(1), vjust = -0.5) +
  labs(title = paste0("scdblfinder doublet_expected_rate: ", doublet_expected_rate, 
                      "\nDoublet rate combined in %",
                      "\n(not including dehashing Singlet but not confident): ", 
                      multiplet_perc_combined)) +
  theme(plot.title = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
p2 <- ggplot(plot_df, aes(is_multiplet, fill = dehash_is_multiplet)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  theme(plot.title = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
print(cowplot::plot_grid(p1, p2, ncol = 2))

# PLOT: Final numbers for each HTO sample to take downstream
if (identical(dehash_mdta_df$Barcode, scdblfinder_mdta_df$Barcode)) {
  
  plot_df$`dehash_is_confident\ndehash_is_multiplet\nis_multiplet` <- interaction(
    plot_df$dehash_is_confident,
    plot_df$dehash_is_multiplet,
    plot_df$is_multiplet
  )
  
  print(
    ggplot(plot_df, aes(interaction(dehash_hashing_var, dehash_assignment), 
                        fill = `dehash_is_confident\ndehash_is_multiplet\nis_multiplet`)) +
      geom_bar(position = position_dodge(1)) +
      geom_text(stat = 'count', aes(label = after_stat(count)), position = position_dodge(1), vjust = -0.5) +
      labs(title = paste0("scdblfinder doublet_expected_rate: ", doublet_expected_rate, 
                          "\nDoublet rate combined in %",
                          "\n(not including hasheddrops Singlet but not confident): ", 
                          multiplet_perc_combined)) +
      facet_grid(. ~ interaction(dehash_hashing_var, dehash_assignment), scales = "free_x") +
      scale_fill_brewer(palette = "Set3", direction = -1) +
      theme(plot.title = element_text(size = 8),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
  )
  
}

dev.off()

# Session info

sessioninfo::session_info(to_file = file.path(out_dir, "session_info.txt"))
sessioninfo::session_info()

# rm(list = ls()); gc()
