# Aggregate single-cell qc per sample/sequencing run results

# Setup

setwd(rprojroot::find_rstudio_root_file())
renv::load()
suppressPackageStartupMessages({
  source(file.path("utils", "env.R"))
  library(assertthat)
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
parser$add_argument("--heatmap_breaks", type = "character", default = "-5,5,0.5", help = "Comma-separated numeric values: min max step for seq()")
parser$add_argument("--heatmap_colors", type = "character", nargs = "+", default = c("blue", "white", "red"), help = "List of colors, e.g. blue white red")
parser$add_argument("--heatmap_fontsize", type = "integer", default = 10, help = "Font size for heatmap numbers (default: 10)")
args <- parser$parse_args()
for (i in seq_along(args)) {assign(names(args)[i], args[[i]])}
out_dir <- create_dir(out_dir)
heatmap_breaks = as.numeric(strsplit(heatmap_breaks, split = ",")[[1]])

# src_dir = "results/process_droplets_pipeline/config/qc_sc_sample"
# out_dir = "."
# heatmap_breaks = "-5,5,0.5"
# heatmap_colors = c("blue", "white", "red")
# heatmap_fontsize = 10

# ----- MAIN -----

seq_run_ids <- list.files(src_dir, full.names = FALSE, recursive = FALSE)
qc_agg_df_lst <- lapply(seq_run_ids, function(seq_run_id) {
  df <- read.csv(file.path(src_dir, seq_run_id, "metrics.csv"))
  df$seq_run_id <- seq_run_id
  return(df)
})

qc_agg_df <- do.call("rbind", qc_agg_df_lst)
write.csv(qc_agg_df, file.path(out_dir, "metrics.csv"), row.names = FALSE, quote = FALSE)

# **(Optional) Exclude some columns from the heatmap**
qc_agg_df = qc_agg_df[
     , grepl(".median|barcode_count|seq_run_id|group_id", colnames(qc_agg_df))
]

heatmap_breaks <- seq(heatmap_breaks[1], heatmap_breaks[2], by = heatmap_breaks[3])
heatmap_colors <- colorRampPalette(heatmap_colors)(length(heatmap_breaks))

# Prepare data for plotting

#qc_agg_df$run <- "a" # In case you want all samples in one heatmap, set run to a constant value

qc_agg_sub_df <- qc_agg_df #%>% 
  #filter(!!sym(split_key) == split_key_lvl)

data_mx <- qc_agg_sub_df %>% 
  select(where(is.numeric)) %>% 
  as.matrix()
rownames(data_mx) <- paste0(qc_agg_sub_df$seq_run_id, "_", qc_agg_sub_df$group_id)

mads_mx <- apply(data_mx, 2, calculate_mads)
#allfinite_colnames <- apply(mads_mx, 2, function(x) all(!is.na(x)))
allfinite_colnames <- colnames(data_mx)
mads_mx[!is.finite(mads_mx)] <- NA
assert_that(identical(dim(mads_mx), dim(data_mx)))

# Plot

split_key = "group_id"
plot_num_rows <- max(table(qc_agg_df[[split_key]]))
pdf(file.path(out_dir, "heatmaps.pdf"), width = (plot_num_rows * (5/6)), height = 8)

for (split_key_lvl in unique(qc_agg_df[[split_key]])) {
  
  is_plot_rows <- grepl(paste0("_", split_key_lvl), rownames(mads_mx), fixed = TRUE)

  # Plot
  display_numbers_mx <- data_mx[is_plot_rows, allfinite_colnames]
  display_numbers_mx <- round(display_numbers_mx, digits = 2)
  mx <- mads_mx[is_plot_rows, allfinite_colnames]
  print(
    pheatmap(t(mx), 
             color = heatmap_colors,
             breaks = heatmap_breaks,
             cluster_cols = FALSE,
             cluster_rows = FALSE,
             display_numbers = t(display_numbers_mx),
             number_format = "%.2f", 
             fontsize_number = heatmap_fontsize,
             na_col = "grey80"
    )
  )

}

dev.off()

# rm(list = ls()); gc()
