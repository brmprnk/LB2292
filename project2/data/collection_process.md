
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

Using the `process_tcga.py` script (make sure to correctly fill out the paths to the data folders and sample sheets), first the individual data types and the clinical metadata are parsed individually (as described in the subsections below). The resulting individual output files are:
- `data/processed/clinical_full.csv`
- `data/processed/cnv_full.pkl`
- `data/processed/expression_full.pkl`
- `data/processed/meth_full.pkl`
- `data/processed/mirna_full.pkl`

From these files, there are $9648$ samples that have all the data types available, so we also create subset files for each type of data that only contain these samples. These files are:

- `data/processed/clinical.csv`
- `data/processed/cnv.pkl`
- `data/processed/expression.pkl`
- `data/processed/meth.pkl`
- `data/processed/mirna.pkl`


### Gene expression processing

The gene expression data was processed as follows:
- We used the transcripts per million (TPM) values from the raw downloaded data files.
- For some gene names, several ENSEMBL IDs were found, so we took the summed value of all values for the same gene name, preserving the TPM totals.
- Next, we log-transformed the data using the formula $log(TPM+1)$.
- NOTE: **gene expression data was NOT normalized**, as for some tests we want to know the original mean of the distribution.

### DNA Methylation processing

The DNA methylation data was processed as follows:
- We used the beta values from the raw downloaded data files.
- To link the probe IDs to gene names, we used the EPICv2 annotation file.
- We then took the mean beta value for all probes linked to the same gene. (They were considered linked if the probe was within 2000bp of the gene's transcription start site).


### Micro RNA processing

Similar to gene expression, the miRNA data was normalized to 1 million counts per sample (CPM), and log-transformed, but not standardized.

### Copy Number Variation processing

We just used the copy_number column from the raw data. No additional processing was performed.


### Clinical Data Processing

The clinical data was processed by parsing each of the XML files separately, and finally combining them into a single pandas dataframe. The parsing of the XML files is implemented in `process_clinical.py`, and this file also includes extensive documentation for the decisions made for each field in the XML file.