## Workflow overview

This workflow is a best-practice workflow for preprocessing multiplexed 10X single-cell RNA-seq data produced by the Cell Ranger `multi` pipeline.
The workflow is built using [Snakemake](https://snakemake.readthedocs.io/en/stable/) and consists of the following steps:

1. Remove empty droplets from the raw feature-barcode matrix (`DropletUtils::emptyDrops`)
2. Demultiplex samples using HTO/hashtag oligonucleotide information (`cellhashR`, consensus across 6 methods)
3. Detect and flag doublets (`scDblFinder`)
4. Filter cells by conventional QC metrics using MAD-based outlier detection (`scuttle::isOutlier`)
5. Compute per-sample QC summary statistics and plots (`qc_sc_sample`)
6. Aggregate QC metrics across all samples and render cross-sample heatmaps (`qc_sc_aggregate`)

## Running the workflow

### Input data

The workflow auto-detects samples as subdirectories of `data/cellranger/`.
Each sample must have a Cell Ranger `multi` output H5 matrix and an HTO-to-sample mapping file.

| Input | Path pattern | Format |
| --- | --- | --- |
| Raw feature-barcode matrix | `data/cellranger/{sample}/outs/multi/count/raw_feature_bc_matrix.h5` | 10X HDF5; must contain Gene Expression and Multiplexing Capture experiments |
| HTO-to-sample mapping | `results/process_droplets/hto_to_sample_mapping/{sample}/hto_to_sample_mapping.tsv` | Tab-separated; columns: `hto_id`, `sample_name` |

### Parameters

This table lists all parameters that can be set in `config/config.yml`.

| Section | Parameter | Type | Description | Default |
| --- | --- | --- | --- | --- |
| **empty_emptydrops** | | | | |
| | `seed_val` | int | Random seed for reproducibility | `827` |
| | `n_cores` | int | Number of cores for `emptyDrops` | `1` |
| | `emptydrops_alpha` | float | FDR threshold for retaining barcodes as non-empty | `0.0001` |
| **dehash_cellhashr** | | | | |
| | `chemistry_10x` | str | 10X Genomics chemistry version (e.g. `10xV3`) | `10xV3` |
| | `methods` | list | Demultiplexing algorithms to run; subset of `[htodemux, multiseq, dropletutils, gmm_demux, bff_raw, bff_cluster]` | all 6 methods |
| **doublet_scdblfinder** | | | | |
| | `seed_val` | int | Random seed for reproducibility | `827` |
| | `n_cores` | int | Number of cores for `scDblFinder` | `1` |
| **filter_conventional** | | | | |
| | `sample_id_var` | str | Column in barcode metadata used to group cells per sample for outlier detection | `hashing_var` |
| | `is_outlier_nmads` | int | Number of MADs beyond which a cell is flagged as an outlier | `3` |
| | `is_outlier_fields_id` | str | Set of QC metric fields to use for outlier detection; `default` applies the built-in field set | `default` |
| **qc_sc_sample** | | | | |
| | `grouping_var` | str | Metadata column used to group cells when computing per-sample QC summaries | `hashing_var` |
| **qc_sc_aggregate** | | | | |
| | `grouping_var` | str | Metadata column used to group cells when aggregating QC metrics across samples | `group_id` |
| | `heatmap_colors` | list | Three colours defining the low / mid / high gradient of the QC heatmap | `[blue, white, red]` |
| | `heatmap_fontsize` | int | Font size used in the QC heatmap | `10` |
