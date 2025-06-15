# Filtering of cells using conventional unimodal method e.g. scuttle::isOutlier(nmads = 3)
# Rscript filter_conventional.R --counts_dir='data/cellranger/SLE_1A_total/raw_feature_bc_matrix' --whitelists "test/whitelist.txt" "test/doublet/whitelist.txt" --sample_id_src_path=test/barcode_metadata.csv --out_dir=test/filter &

setwd(rprojroot::find_rstudio_root_file())
renv::load()
suppressPackageStartupMessages({
  source(file.path("utils", "env.R"))
  source(file.path("utils", "sc_helpers.R"))
  source(file.path("utils", "wrapper.R"))
  library(argparse)
  library(tidyverse)
  library(qs)
  library(scuttle)
  library(Seurat)
})

# Parameters

parser <- ArgumentParser(description = "Filter cells using conventional unimodal method e.g. scuttle::isOutlier(nmads = 3)")
parser$add_argument("--counts_dir", type = "character", default = file.path("data", "cellranger", "sle_batch2_pool2", "raw_feature_bc_matrix"),
                    help = "Path to counts directory (default: data/cellranger/sle_batch2_pool2/raw_feature_bc_matrix)")
parser$add_argument("--whitelists", type = "character", nargs = "+", required = TRUE, help = "Paths to barcode whitelists (required)")
parser$add_argument("--sample_id_src_path", type = "character", required = TRUE, help = "Path to sample ID source CSV file (required)")
parser$add_argument("--out_dir", type = "character", required = TRUE, help = "Output directory (required)")
parser$add_argument("--sample_id_var", type = "character", default = "hashing_var",
                    help = "Variable name for sample ID in cell metadata (default: hashing_var)")
parser$add_argument("--is_outlier_nmads", type = "numeric", default = 3,
                    help = "Number of MADs for outlier detection (default: 3)")
parser$add_argument("--is_outlier_fields_id", type = "character", default = "default",
                    help = "ID for outlier fields (default: default, options: default, include_adt)")
args <- parser$parse_args()
for (i in seq_along(args)) {assign(names(args)[i], args[[i]])}
out_dir <- create_dir(out_dir)

# Parameters

#seq_run_id = seq_run_id #"sle_batch2_pool2"
#counts_dir = file.path("data", "cellranger", seq_run_id, "raw_feature_bc_matrix")
#whitelists = c(file.path("results/process_droplets/dehash/hasheddrops_emptydrops_alpha0.001_seed827/",
#                         seq_run_id, "whitelist.txt"),
#               file.path("results/process_droplets/multiplet/scdblfinder_seed827_doublet_expected_rate",
#                         seq_run_id, "whitelist.txt")
#               )
# sample_id is used to split into groups when getting outliers. See wrapper/get_outliers
#sample_id_src_path = file.path("results/process_droplets/dehash/hasheddrops_emptydrops_alpha0.001_seed827", 
#                               seq_run_id, "barcode_metadata.csv")
#sample_id_var = "hashing_var"
#is_outlier_nmads = 3
#is_outlier_fields_id = "default"

# Outlier definition
if (is_outlier_fields_id == "default") {
  is_outlier_fields = NULL
} else if (is_outlier_fields_id == "include_adt") {
  is_outlier_fields = c("low_lib_size",
                        "low_n_features",
                        "high_subsets_mito_genes_percent",
                        "altexps_Antibody Capture_sum.low_lib_size",
                        "altexps_Antibody Capture_detected.low_n_features")
  
}

#out_dir = create_dir(file.path("results", "process_droplets", "filter", 
#                               paste0("conventional_", is_outlier_fields_id),
#                               seq_run_id))

#------ MAIN -----

pdf(file.path(out_dir, "plots.pdf"), height = 9, width = 50)
par(mfrow = c(1,1))

# Get barcode whitelist (cell singlets)
barcode_singlets <- Reduce(intersect, lapply(whitelists, function(pth) readLines(pth)))

# Create sce with cell singlets only

sce <- read_sc_data(counts_dir)
is_hto <- grepl("HTO", rowData(sce)$ID, ignore.case = FALSE)
is_antibody_capture <- grepl("Antibody Capture", rowData(sce)$Type)
rowData(sce)$Type[is_hto & is_antibody_capture] <- "Multiplexing Capture"
sce <- SingleCellExperiment::splitAltExps(sce, rowData(sce)$Type)
# Split experiments
sce <- sce[, barcode_singlets]

# Add hashing_var to cell metadata
sample_id_barcode_mapping <- read_csv(sample_id_src_path) %>% 
  select(c("Barcode", sample_id_var)) %>% 
  tibble::deframe()
colData(sce) <- cbind(colData(sce),
                      sample_id = unname(sample_id_barcode_mapping[colData(sce)$Barcode]))

# Add qc metrics for all relevant experiments e.g. RNA + ADT

# Get uniqified gene symbols as in Seurat object
feature_names_uniqified <- rownames(make_uniqueAsInSeu(sce))
sce <- do_basicQC(counts_mat = NULL,
                  sce = sce,
                  feature_names = feature_names_uniqified,
                  use.altexps = TRUE)

plot_basicQC(colData(sce),
             group = "sample_id", 
             metrics = names(which(unlist(lapply(colData(sce), function(x) is.numeric(x))))), 
             transform_log10 = TRUE)

# Get outliers

outlier_mx <- get_outliers(x = colData(sce), 
                           sample_ids = sce$sample_id,
                           is_outlier_fields = is_outlier_fields,
                           is_outlier_nmads = is_outlier_nmads)
barcodes <- rownames(outlier_mx)
outlier_mx <- apply(outlier_mx, MARGIN = 2, as.numeric)
rownames(outlier_mx) <- barcodes
# table(outlier_mx[,"outlier"], sce[[sample_id_var]])

# Save output

writeLines(rownames(outlier_mx)[!outlier_mx[ ,"outlier"]], file.path(out_dir, "whitelist.txt"))
outlier_df <- outlier_mx %>% as.data.frame() %>% 
  rownames_to_column(var = "Barcode")
write.csv(outlier_df, file.path(out_dir, "barcode_metadata.csv"), row.names = FALSE, quote = FALSE)

dev.off()

# Session info

sessioninfo::session_info(to_file = file.path(out_dir, "session_info.txt"))
sessioninfo::session_info()

# rm(list = ls()); gc()
