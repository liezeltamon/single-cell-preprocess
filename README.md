# Snakemake workflow: `single-cell-preprocess`

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/liezeltamon/single-cell-preprocess/actions/workflows/main.yml/badge.svg?branch=main)](https://github.com/liezeltamon/single-cell-preprocess/actions/workflows/main.yml)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![workflow catalog](https://img.shields.io/badge/Snakemake%20workflow%20catalog-darkgreen)](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/liezeltamon/single-cell-preprocess)

A Snakemake workflow for preprocessing multiplexed 10X single-cell RNA-seq data: empty droplet removal, HTO-based sample demultiplexing, doublet detection, QC filtering, and quality metric reporting.

- [Snakemake workflow: `single-cell-preprocess`](#snakemake-workflow-single-cell-preprocess)
  - [Overview](#overview)
  - [Input data](#input-data)
  - [Output](#output)
  - [Usage](#usage)
  - [Deployment options](#deployment-options)
  - [Authors](#authors)
  - [References](#references)

## Overview

The workflow is built using [Snakemake](https://snakemake.readthedocs.io/en/stable/) and consists of the following steps:

1. **Empty droplet removal** — `DropletUtils::emptyDrops()` distinguishes cell-containing droplets from empty droplets using the raw Cell Ranger count matrix, retaining barcodes that pass the configured FDR threshold.
2. **Sample demultiplexing** — `cellhashR` runs up to six HTO-based demultiplexing algorithms (HTODemux, MultiSeqDemux, DropletUtils, GMM-Demux, BFF-raw, BFF-cluster) and calls a consensus singlet assignment per barcode by majority vote.
3. **Doublet detection** — `scDblFinder` simulates artificial doublets in PCA space and classifies each droplet as a singlet or doublet using a random forest classifier.
4. **Conventional QC filtering** — Per-sample MAD-based outlier removal on library size, feature count, and mitochondrial fraction using `scuttle::isOutlier()`.
5. **Per-sample QC metrics** — Computes and plots summary statistics (median, mean, SD, min, max) per sample or grouping variable.
6. **Aggregate QC** — Aggregates QC metrics across all samples, normalises values to MADs, and renders cross-sample heatmaps for a global quality overview.

Detailed information about input data and workflow configuration can be found in the [`config/README.md`](config/README.md).

## Input data

The workflow expects **Cell Ranger `multi` pipeline** outputs. Samples are auto-detected as subdirectories of `data/cellranger/`.

| Input | Path | Notes |
| --- | --- | --- |
| Raw feature-barcode matrix | `data/cellranger/{sample}/outs/multi/count/raw_feature_bc_matrix.h5` | Must contain both Gene Expression and Multiplexing Capture (HTO) libraries |
| HTO-to-sample mapping | `results/process_droplets/hto_to_sample_mapping/{sample}/hto_to_sample_mapping.tsv` | Tab-separated; columns: `hto_id`, `sample_name` |

## Output

All outputs are written to `results/process_droplets_pipeline/config/`.

| Directory | Key output files |
| --- | --- |
| `empty/{sample}/` | `whitelist.txt`, `blacklist.txt`, `output.qs`, `plots.pdf`, `session_info.txt` |
| `dehash/{sample}/` | `whitelist.txt`, `barcode_metadata.csv`, `metrics.csv`, `output.csv`, `plots.pdf`, `session_info.txt` |
| `doublet/{sample}/` | `whitelist.txt`, `barcode_metadata.csv`, `output.qs`, `plots.pdf`, `session_info.txt` |
| `filter/{sample}/` | `whitelist.txt`, `barcode_metadata.csv`, `plots.pdf`, `session_info.txt` |
| `qc_sc_sample/{sample}/` | `metrics.csv`, `plots.pdf` |
| `qc_sc_aggregate/` | `metrics.csv`, `heatmaps.pdf` |

`whitelist.txt` files at each step carry the set of high-quality barcodes surviving that step; `barcode_metadata.csv` files carry per-cell annotations accumulated across steps.

## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/liezeltamon/single-cell-preprocess).

If you use this workflow in a paper, please cite the repository URL or its DOI and the tools listed in the [References](#references) section.

## Deployment options

Change to the workflow directory and adjust options in `config/config.yml`.

```bash
cd path/to/single-cell-preprocess
```

Perform a dry run to check the workflow before execution:

```bash
snakemake --dry-run
```

Run with test files using **conda**:

```bash
snakemake --cores 2 --sdm conda --directory .test
```

Run with **apptainer** / **singularity**:

```bash
snakemake --cores 2 --sdm conda apptainer --directory .test
```

Run on an HPC cluster via **SLURM** (recommended for production):

```bash
# Load required modules first
module load R/4.3.2-gfbf-2023a

sbatch -J process_droplets_pipeline -p short,long \
  --mem=80G --cpus-per-task=4 \
  --output=%x.log.out --error=%x.log.err \
  --wrap="snakemake -s Snakefile --cores 4 --rerun-incomplete"
```

## Authors

- Liezel Tamon
  - University of Oxford
  - [ORCID profile](https://orcid.org/0000-0003-3705-6019)
  - https://www.imm.ox.ac.uk/people/liezel-tamon

## References

> Köster, J., Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., & Nahnsen, S. _Sustainable data analysis with Snakemake_. F1000Research, 10:33, **2021**. https://doi.org/10.12688/f1000research.29032.2

> Lun, A. T. L., Riesenfeld, S., Andrews, T., Dao, T. P., Gomes, T., & Marioni, J. C. _EmptyDrops: distinguishing cells from empty droplets in droplet-based single-cell RNA sequencing data_. Genome Biology, 20:63, **2019**. https://doi.org/10.1186/s13059-019-1662-y

> Bimber, B. N., & Cole, B. L. _CellHashR: an R package for cell hashing demultiplexing_. bioRxiv, **2021**. https://doi.org/10.1101/2021.09.03.458947

> Germain, P.-L., Lun, A., Garcia Meixide, C., Macnair, W., & Robinson, M. D. _Doublet identification in single-cell sequencing data using scDblFinder_. F1000Research, 10:979, **2021**. https://doi.org/10.12688/f1000research.73600.2

> McCarthy, D. J., Campbell, K. R., Lun, A. T. L., & Willis, Q. F. _Scater: pre-processing, quality control, normalisation and visualisation of single-cell RNA-seq data in R_. Bioinformatics, 33(8), 1179–1186, **2017**. https://doi.org/10.1093/bioinformatics/btw777
