
# TCGA Data Collection Process

This document describes all the steps taken to collect the TCGA data from the GDC portal.


## Downloading Data from the GDC Portal

Manifests were downloaded on Feb 12th 2025, from the GDC portal (https://portal.gdc.cancer.gov/) as follows:

- In the cohort builder pane, select all TCGA samples
- In the repositories pane, select only `Open Access` data
- In the repositories pane, select each of these data types separately, and download the manifest for each, as well as the sample sheet file (sample sheet not necessary for the clinical data):
    - **Clinical Data**: `Data Category: clinical` + `Data Type: Clinical Supplement` + `Data Format: bcr xml`, for the XML files with the metadata (11,167 files, 11,167 cases, 510.73MB, stored in `data/manifests/gcd_manifest_tcga_open_clinical.txt`)

    - **Copy Number Variations**: `Data Category: copy number variation` + `Data Type: Gene Level Copy Number` + `Workflow Type: ASCAT3`, for the CNV data (10,632 files, 10,632 cases, 36.64GB, stored in `gdc_manifest_tcga_open_gene-level-cn_ascat3.txt`)
    
    - **Gene Expression**: `Data Category: transcriptome profiling` + `Data Type: Gene Expression Quantification`, for the gene expression data (11,499 files, 10,511 cases, 48.72GB, stored in `data/manifests/gcd_manifest_tcga_open_expression.txt`)
        - NOTE: We have more files than cases here, because some samples have both tumor and normal tissue expression.
    
    - **Micro RNA Expression**: `Data Category: transcriptome profiling` + `Data Type: miRNA Expression Quantification` + `Data Format: txt`, for the miRNA expression data (11,082 files, 10,250 cases, 557.23MB stored in `data/manifests/gcd_manifest_tcga_open_mirna.txt`)

    - **DNA Methylation**: `Data Category: DNA methylation` + `Data Type: Methylation Beta Value`, for the DNA methylation data (12,527 files, 11,000 cases, 131.93 GB, stored in `data/manifests/gcd_manifest_tcga_open_meth.txt`)

Next, the files were downloaded using the `gdc-client` tool, which was downloaded from https://gdc.cancer.gov/access-data/gdc-data-transfer-tool, and was then run in the terminal with the following command for each manifest:

```bash
gdc-client download -m data/manifests/{file}.txt -d data/raw/{data_type}/
```

## Data Processing

Using the `process_clinical.py` script, first the individual data types and the clinical metadata are parsed individually, and then they are aligned and combined (TODO).

The resulting individual output files are:
- `data/processed/clinical.csv`
- `data/processed/cnv.pkl`
- `data/processed/expression.pkl`
- `data/processed/mirna.pkl` (TODO)
- `data/processed/meth.pkl`

From these files, there are ... samples that have all the data types available (TODO).


## Things to check

- Some data is only available for a subset of the samples. How will we handle missing data? Do we drop samples that don't have all the data types available?