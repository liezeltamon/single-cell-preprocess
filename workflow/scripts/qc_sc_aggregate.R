# Aggregate single-cell qc per sample/sequencing run results

# Setup

setwd(rprojroot::find_rstudio_root_file())
renv::load()
suppressPackageStartupMessages({
  source(file.path("utils", "env.R"))
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

parser <- ArgumentParser(description = "Aggregate single-cell quality control results per sameple / sequencing run")
parser$add_argument("--src_dir", type = "character", required = TRUE, help = "Path to directory containing single-cell qc results")
parser$add_argument("--out_dir", type = "character", required = TRUE, help = "Output directory")
#parser$add_argument("--heatmap_breaks", type = "double", nargs = 3, default = c(-3, 3, 0.5), help = "Three numeric values: min max step for seq() (default: -3 3 0.5)")
parser$add_argument("--heatmap_breaks", type = "character", default = "-3,3,0.5", help = "Comma-separated numeric values: min max step for seq()")
parser$add_argument("--heatmap_colors", type = "character", nargs = "+", default = c("blue", "white", "red"), help = "List of colors, e.g. blue white red")
parser$add_argument("--heatmap_fontsize", type = "integer", default = 10, help = "Font size for heatmap numbers (default: 10)")
args <- parser$parse_args()
for (i in seq_along(args)) {assign(names(args)[i], args[[i]])}
out_dir <- create_dir(out_dir)
heatmap_breaks = as.numeric(strsplit(heatmap_breaks, split = ",")[[1]])

# ----- MAIN -----

seq_run_ids <- list.files(src_dir, full.names = FALSE, recursive = FALSE)
qc_agg_df_lst <- lapply(seq_run_ids, function(seq_run_id) {
  read.csv(file.path(src_dir, seq_run_id, "metrics.csv"))
})

qc_agg_df <- do.call("rbind", qc_agg_df_lst)
write.csv(qc_agg_df, file.path(out_dir, "metrics.csv"))

# **(Optional) Exclude some columns from the heatmap**
qc_agg_df = qc_agg_df[
     , grepl("-median|barcode_count|seq_run_id|group_id", colnames(qc_agg_df))
]

breaks <- seq(args$breaks[1], args$breaks[2], by = args$breaks[3])
custom_colors <- colorRampPalette(heatmap_colors)(length(breaks))
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
