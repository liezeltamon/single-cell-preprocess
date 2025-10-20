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
parser$add_argument("--grouping_var_src_path", type = "character", required = TRUE, help = "Cell metadata CSV")
parser$add_argument("--grouping_var", type = "character", default = "hashing_var", help = "Variable name for group ID in cell metadata (default: hashing_var)")
parser$add_argument("--whitelist_path", type = "character", required = TRUE, help = "Path to barcode whitelist file")
parser$add_argument("--out_dir", type = "character", required = TRUE, help = "Output directory")
args <- parser$parse_args()
for (i in seq_along(args)) {assign(names(args)[i], args[[i]])}
out_dir <- create_dir(out_dir)

# Functions

collapse_unique <- function(x) {
  unique_x <- unique(na.omit(x))
  if (length(unique_x) == 0) NA_character_ else paste(unique_x, collapse = ";")
}

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

# Create cell metadata including column for grouping_var and other covariates
barcode_metadata_df <- read_csv(grouping_var_src_path)
rownames(barcode_metadata_df) <- barcode_metadata_df$Barcode

# Match order of barcodes in sce
barcode_metadata_df <- barcode_metadata_df[rownames(colData(sce)), , drop = FALSE]
colData(sce) <- DataFrame(barcode_metadata_df)
sce$grouping_var <- sce[[grouping_var]]
assert_that(!any(duplicated(colnames(colData(sce)))))

# Add qc metrics for all relevant experiments e.g. RNA + ADT

# Get uniqified gene symbols as in Seurat object
feature_names_uniqified <- rownames(make_uniqueAsInSeu(sce))
sce <- do_basicQC(
    counts_mat = NULL,
    sce = sce,
    feature_names = feature_names_uniqified,
    use.altexps = TRUE
)

metadata_df <- as.data.frame(colData(sce))
write.csv(
    metadata_df, file.path(out_dir, "barcode_metadata.csv"), row.names = FALSE, quote = FALSE
)

pdf(file.path(out_dir, "plots.pdf"), height = 9, width = 50)

plot_basicQC(
    colData(sce),
    group = "grouping_var", 
    metrics = names(which(unlist(lapply(colData(sce), function(x) is.numeric(x))))), 
    transform_log10 = TRUE
)

dev.off()

qc_per_group_df <- colData(sce) %>%
    as.data.frame() %>%
    group_by(grouping_var) %>%
    summarise(
        across(
            where(~ !is.numeric(.x)),
            ~ collapse_unique(.x),
            .names = "{.col}"
        ),
        across(
            where(is.numeric),
            list(
                # na.rm = FALSE in of unexpected NAs
                median = ~ median(.x, na.rm = FALSE),
                mean = ~ mean(.x, na.rm = FALSE),
                min = ~ min(.x, na.rm = FALSE),
                max = ~ max(.x, na.rm = FALSE),
                sd = ~ sd(.x, na.rm = FALSE)
            ),
            .names = "{.col}-{.fn}"
        ),
        barcode_count = n()
    ) %>%
    as.data.frame()

write.csv(qc_per_group_df, file.path(out_dir, "metrics.csv"))

# rm(list = ls()); gc()
