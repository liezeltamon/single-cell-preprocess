# Single-cell qc
# Rscript qc_sc.R --parent_counts_dir='data/cellranger_dir1' --group_id_dir='test/dehash' --group_id_var='hashing_var' --whitelist_dir='test/filter' --out_dir='test/qc_sc' &

# Setup

renv::load()
setwd(rprojroot::find_rstudio_root_file())
suppressPackageStartupMessages({
  source(file.path("utils", "env.R"))
  source(file.path("utils", "sc_helpers.R"))
  source(file.path("utils", "wrapper.R"))
  library(argparse)
  library(dplyr)
  library(pheatmap)
})

calculate_mads <- function(x) {
  med <- median(x, na.rm = TRUE)
  mad <- mad(x, constant = 1, na.rm = TRUE)
  mads <- (x - med) / mad
  return(mads)
}

# Parameters

parser <- ArgumentParser(description = "Single-cell quality control")
parser$add_argument("--parent_counts_dir", type = "character", required = TRUE, help = "Path to directory with 10x counts directories")
parser$add_argument("--group_id_dir", type = "character", required = TRUE, help = "Path to directory with metadata CSV")
parser$add_argument("--group_id_var", type = "character", default = "hashing_var", help = "Variable name for group ID in cell metadata (default: hashing_var)")
parser$add_argument("--whitelist_dir", type = "character", required = TRUE, help = "Path to directory with barcode whitelist files")
parser$add_argument("--out_dir", type = "character", required = TRUE, help = "Output directory")
args <- parser$parse_args()
for (i in seq_along(args)) {assign(names(args)[i], args[[i]])}
out_dir <- create_dir(out_dir)

# ----- MAIN -----

seq_run_ids <- list.files(parent_counts_dir, full.names = FALSE, recursive = FALSE)
qc_agg_df_lst <- lapply(seq_run_ids, function(seq_run_id) {
  
  counts_dir = file.path(parent_counts_dir, seq_run_id, "raw_feature_bc_matrix")
  sce <- read_sc_data(counts_dir)
  # Split experiments
  is_hto <- grepl("HTO", rowData(sce)$ID, ignore.case = FALSE)
  is_antibody_capture <- grepl("Antibody Capture", rowData(sce)$Type)
  rowData(sce)$Type[is_hto & is_antibody_capture] <- "Multiplexing Capture"
  sce <- SingleCellExperiment::splitAltExps(sce, rowData(sce)$Type)
  
  # Get barcode whitelist (cell singlets)
  barcode_singlets <- readLines(file.path(whitelist_dir, seq_run_id, "whitelist.txt"))
  sce <- sce[, barcode_singlets]
  
  # Add hashing_var to cell metadata
  group_id_barcode_mapping <- read_csv(file.path(group_id_dir, seq_run_id, "barcode_metadata.csv")) %>% 
    select(c("Barcode", "hashing_var")) %>% 
    tibble::deframe()
  colData(sce) <- cbind(
    colData(sce),
    hashing_var = unname(group_id_barcode_mapping[colData(sce)$Barcode])
  )
  sce$seq_run_id <- seq_run_id
  sce$group_id <- sce[[group_id_var]]
  
  # Add qc metrics for all relevant experiments e.g. RNA + ADT
  
  # Get uniqified gene symbols as in Seurat object
  feature_names_uniqified <- rownames(make_uniqueAsInSeu(sce))
  sce <- do_basicQC(counts_mat = NULL,
                    sce = sce,
                    feature_names = feature_names_uniqified,
                    use.altexps = TRUE)
  
  out_run_dir = create_dir(file.path(out_dir, seq_run_id))
  pdf(file.path(out_run_dir, "plots.pdf"), height = 9, width = 50)
  
  plot_basicQC(colData(sce),
               group = "group_id", 
               metrics = names(which(unlist(lapply(colData(sce), function(x) is.numeric(x))))), 
               transform_log10 = TRUE)
  
  dev.off()
  
  write.csv(as.data.frame(colData(sce)), 
            file.path(create_dir(file.path(out_dir, seq_run_id)), "metrics.csv"),
            row.names = FALSE)
  
  qc_agg_df <- colData(sce) %>%
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
    mutate(seq_run_id = seq_run_id) %>% 
    as.data.frame()
  
  return(qc_agg_df)
  
})

qc_agg_df <- do.call("rbind", qc_agg_df_lst)
out_dir <- create_dir(file.path(out_dir, "aggregate"))
write.csv(qc_agg_df, file.path(out_dir, "metrics.csv"))

# **(Optional) Exclude some columns from the heatmap**
#qc_agg_df_orig = qc_agg_df
#qc_agg_df = qc_agg_df_orig
qc_agg_df = qc_agg_df[ ,grepl("-median|barcode_count|seq_run_id|group_id",
                               colnames(qc_agg_df))]
heatmap_fontsize = 10

# Define the breaks for the pheatmap colour palette
breaks <- seq(-3, 3, by = 0.5)
custom_colors <- colorRampPalette(c("blue", "white", "red"))(length(breaks))

pdf(file.path(out_dir, "heatmaps.pdf"), height = 4, width = 10)

for (group_id_lvl in unique(qc_agg_df$group_id)) {
  
  qc_agg_sub_df <- qc_agg_df %>% 
    filter(group_id == group_id_lvl)
  
  data_mx <- qc_agg_sub_df %>% 
    select(where(is.numeric)) %>% 
    as.matrix()
  rownames(data_mx) <- paste0(qc_agg_sub_df$seq_run_id, "_", qc_agg_sub_df$group_id)

  mads_mx <- apply(data_mx, 2, calculate_mads)
  #allfinite_colnames <- apply(mads_mx, 2, function(x) all(!is.na(x)))
  allfinite_colnames <- colnames(data_mx)
  mads_mx[!is.finite(mads_mx)] <- NA
  
  # Generate the heatmap using pheatmap with custom colors and breaks
  display_numbers_mx <- data_mx[ ,allfinite_colnames]
  display_numbers_mx <- round(display_numbers_mx, digits = 2)
  mx <- mads_mx[ , allfinite_colnames]
  print(
    pheatmap(mx, 
             color = custom_colors, 
             breaks = breaks,
             cluster_cols = FALSE,
             cluster_rows = FALSE,
             display_numbers = display_numbers_mx,
             number_format = "%.2f", 
             fontsize_number = heatmap_fontsize,
             na_col = "grey80"
    )
  )

}

dev.off()

# rm(list = ls()); gc()
