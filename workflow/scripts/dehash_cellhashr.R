# Dehash with multiple algorithms using cellhashR
# "Rscript dehash_cellhashr.R --counts_dir="data/cellranger/SLE_1A_total/raw_feature_bc_matrix" --out_dir=test --whitelist_path=results/process_droplets/empty/emptydrops_alpha0.001_seed827/SLE_1A_total/whitelist.txt --htotosampletsv_path=results/process_droplets/hto_to_sample_mapping/SLE_1A_total/hto_to_sample_mapping.tsv"
#sbatch -J dehash_SLE_1A \
#  --output=%x_%j.log.out \
#  --error=%x_%j.log.err \
#  --wrap="Rscript dehash_cellhashr.R --counts_dir='data/cellranger/SLE_1A_total/raw_feature_bc_matrix' --out_dir=test --whitelist_path=results/process_droplets/empty/emptydrops_alpha0.001_seed827/SLE_1A_total/whitelist.txt --htotosampletsv_path=results/process_droplets/hto_to_sample_mapping/SLE_1A_total/hto_to_sample_mapping.tsv"

setwd(rprojroot::find_rstudio_root_file())
renv::load()
suppressPackageStartupMessages({
  source(file.path("utils", "env.R"))
  library(DropletUtils)
  library(Seurat)
  source(file.path("utils", "sc_helpers.R"))
  library(assertthat)
  library(argparse)
  library(tidyverse)
  library(cellhashR)
})

# Parameters

parser <- ArgumentParser(description = "Dehash multiplexed samples with cellhashR")
parser$add_argument("--counts_dir", type = "character", required = TRUE,
                    help = "10x counts directory")
parser$add_argument("--out_dir", type = "character", required = TRUE,
                    help = "Output directory")
parser$add_argument("--whitelist_path", type = "character", required = TRUE,
                    help = "Path to cell barcode whitelist")
parser$add_argument("--chemistry_10x", type = "character", default = "10xV3",
                    help = "10x chemistry type (default: 10xV3)")
parser$add_argument("--methods", type = "character", nargs = "+",
                    # demuxem - installed but not detected
                    # demuxmix - very slow
                    default = c("htodemux", "multiseq", "dropletutils", "gmm_demux", "bff_raw", "bff_cluster"),
                    help = "Methods to use for dehashing, see ?cellhashR::GenerateCellHashingCalls for options")
parser$add_argument("--htotosampletsv_path", type = "character", required = TRUE,
                    help = "Path to HTO to sample mapping TSV file")
args <- parser$parse_args()
for (i in seq_along(args)) {assign(names(args)[i], args[[i]])}
out_dir <- create_dir(out_dir)

#counts_dir = "data/cellranger/SLE_1A_total/raw_feature_bc_matrix"
#counts_dir = "code/process_droplets/tmp/cell_type_counts.csv"
#whitelist_path = "results/process_droplets_pipeline/emptydrops_alpha00010_n_cores2_seed_val827/SLE_1A_total/whitelist.txt"
#chemistry_10x = "10xV3"
#methods = c("htodemux", "demuxmix")
#htotosampletsv_path = "results/process_droplets/hto_to_sample_mapping/SLE_1A_total/hto_to_sample_mapping.tsv"

# ----- MAIN -----

# Load to sce

sce <- read_sc_data(counts_dir)
is_hto <- grepl("HTO", rowData(sce)$ID, ignore.case = FALSE)
is_antibody_capture <- grepl("Antibody Capture", rowData(sce)$Type)
rowData(sce)$Type[is_hto & is_antibody_capture] <- "Multiplexing Capture"
hto_names <- rowData(sce)$ID[is_hto]
# Split experiments
sce <- SingleCellExperiment::splitAltExps(sce, rowData(sce)$Type)

pdf("plots.pdf", width = 16, height = 9)

# QC

whitelist <- readLines(whitelist_path)
whitelist_stripped <- unlist(lapply(
  strsplit(whitelist, split = "-", fixed = TRUE),
  function(x) x[1]
))
# Run to get QC plots
barcode_data <- cellhashR::ProcessCountMatrix(
  rawCountData = counts_dir,
  # Because cell barcode whitelist is given containing called non-empty droplets by DropletUtils::emptyDrops()
  minCountPerCell = 0,
  barcodeWhitelist = hto_names,
  cellbarcodeWhitelist = whitelist_stripped,
  datatypeName = "Antibody Capture"
)
#saveRDS(barcode_data, file.path(out_dir, "barcode_data.rds"))
#PlotNormalizationQC(barcode_data)

# Hashing

mx <- as.matrix(counts(altExp(sce, "Multiplexing Capture")))
mx <- mx[hto_names, whitelist]
#mx_tmp <- mx; colnames(barcode_data) <- colnames(mx_tmp) <- NULL
#assert_that(identical(mx_tmp, barcode_data))
assert_that(identical(dim(barcode_data), dim(mx)))
assert_that(identical(rownames(barcode_data), rownames(mx)))

df <- GenerateCellHashingCalls(
  barcodeMatrix = mx,
  methods = methods,
  #majorityConsensusThreshold = (3 / 4), # Allow 4 out of 5 methods agreeing, or 3 out of 4 (in case 1 method returns negative) - might be too strict
  metricsFile = file.path(out_dir, paste0("metrics.csv")),
  doTSNE = FALSE,
  chemistry = chemistry_10x,
  rawFeatureMatrixH5 = paste0(counts_dir, ".h5")
)
write.csv(
  df, file = file.path(out_dir, "output.csv"), row.names = FALSE, quote = FALSE
)

dev.off()

# Generate barcode_metadata.csv and whitelist.txt outputs

htotosampletsv_df <- read.table(
  htotosampletsv_path, col.names = c("hashing_var", "hto_id"), row.names = "hto_id"
)

df <- df %>%
  rename(
    Barcode = cellbarcode,
    assignment = consensuscall
  )
df <- df %>%
  select(Barcode, assignment) %>%
  mutate(
    hashing_var = htotosampletsv_df[df$assignment, "hashing_var"],
    score = NA,
    is_multiplet = !(df$assignment %in% hto_names),
    # To add column, consider all singlets as confident calls
    is_confident = !is_multiplet
  )
write.csv(df, file.path(out_dir, "barcode_metadata.csv"), row.names = FALSE)

whitelist <- df %>% 
  filter(!is_multiplet) %>% 
  pull(Barcode)
writeLines(whitelist, file.path(out_dir, "whitelist.txt"))

sessioninfo::session_info(to_file = file.path(out_dir, "session_info.txt"))
sessioninfo::session_info()

# rm(list = ls()); gc()
