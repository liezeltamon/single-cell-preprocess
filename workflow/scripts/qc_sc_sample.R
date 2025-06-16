# Single-cell qc

# Setup

setwd(rprojroot::find_rstudio_root_file())
renv::load()
suppressPackageStartupMessages({
  source(file.path("utils", "env.R"))
  source(file.path("utils", "sc_helpers.R"))
  source(file.path("utils", "wrapper.R"))
  library(argparse)
  library(dplyr)
  library(pheatmap)
})

# Parameters

parser <- ArgumentParser(description = "Single-cell quality control")
parser$add_argument("--counts_dir", type = "character", required = TRUE, help = "Path to 10x count directory")
parser$add_argument("--group_id_src_path", type = "character", required = TRUE, help = "Cell metadata CSV")
parser$add_argument("--group_id_var", type = "character", default = "hashing_var", help = "Variable name for group ID in cell metadata (default: hashing_var)")
parser$add_argument("--whitelist_path", type = "character", required = TRUE, help = "Path to barcode whitelist file")
parser$add_argument("--out_dir", type = "character", required = TRUE, help = "Output directory")
args <- parser$parse_args()
for (i in seq_along(args)) {assign(names(args)[i], args[[i]])}
out_dir <- create_dir(out_dir)

# ----- MAIN -----

sce <- read_sc_data(counts_dir)
# Split experiments
is_hto <- grepl("HTO", rowData(sce)$ID, ignore.case = FALSE)
is_antibody_capture <- grepl("Antibody Capture", rowData(sce)$Type)
rowData(sce)$Type[is_hto & is_antibody_capture] <- "Multiplexing Capture"
sce <- SingleCellExperiment::splitAltExps(sce, rowData(sce)$Type)

# Get barcode whitelist (cell singlets)
barcode_singlets <- readLines(whitelist_path)
sce <- sce[, barcode_singlets]

# Add hashing_var to cell metadata
group_id_barcode_mapping <- read_csv(group_id_src_path) %>% 
    select(c("Barcode", !!sym(group_id_var))) %>% 
    tibble::deframe()
colData(sce)[[group_id_var]] <- unname(group_id_barcode_mapping[colData(sce)$Barcode])
#sce$seq_run_id <- seq_run_id
sce$group_id <- sce[[group_id_var]]

# Add qc metrics for all relevant experiments e.g. RNA + ADT

# Get uniqified gene symbols as in Seurat object
feature_names_uniqified <- rownames(make_uniqueAsInSeu(sce))
sce <- do_basicQC(
    counts_mat = NULL,
    sce = sce,
    feature_names = feature_names_uniqified,
    use.altexps = TRUE
)

pdf(file.path(out_dir, "plots.pdf"), height = 9, width = 50)

plot_basicQC(
    colData(sce),
    group = "group_id", 
    metrics = names(which(unlist(lapply(colData(sce), function(x) is.numeric(x))))), 
    transform_log10 = TRUE
)

dev.off()

qc_per_group_df <- colData(sce) %>%
    as.data.frame() %>%
    group_by(group_id) %>%
    summarise(across(where(is.numeric),
                        list(
                            median = ~ median(.x, na.rm = FALSE),
                            mean = ~ mean(.x, na.rm = TRUE),
                            min = ~ min(.x, na.rm = TRUE),
                            max = ~ max(.x, na.rm = TRUE),
                            sd = ~ sd(.x, na.rm = TRUE)
                        ),
                        .names = "{.col}-{.fn}"),
                barcode_count = n()) %>%
    mutate(seq_run_id = unique(sce$Sample)) %>%
    as.data.frame()

write.csv(qc_per_group_df, file.path(out_dir, "metrics.csv"))

# rm(list = ls()); gc()
