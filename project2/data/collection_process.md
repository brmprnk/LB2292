
# TCGA Data Collection Process

This document describes all the steps taken to collect the TCGA data from the GDC portal.


## Downloading Data from the GDC Portal

Manifests were downloaded on Feb 12th 2025, from the GDC portal (https://portal.gdc.cancer.gov/) as follows:

- In the cohort builder pane, select all TCGA samples
- In the repositories pane, select only `Open Access` data
- In the repositories pane, select each of these data types separately, and download the manifest for each:
    - **Clinical Data**: `Data Category: clinical` + `Data Type: Clinical Supplement` + `Data Format: bcr xml`, for the XML files with the metadata (11,167 files, 11,167 cases, 510.73MB, stored in `data/manifests/gcd_manifest_tcga_open_clinical.txt`)

    - **Copy Number Variations**: `Data Category: copy number variation` + `Data Type: Gene Level Copy Number` + `Workflow Type: ASCAT3`, for the CNV data (10,632 files, 10,632 cases, 36.64GB, stored in `gdc_manifest_tcga_open_gene-level-cn_ascat3.txt`)
    
    - **Gene Expression**: `Data Category: transcriptome profiling` + `Data Type: Gene Expression Quantification`, for the gene expression data (11,499 files, 10,511 cases, 48.72GB, stored in `data/manifests/gcd_manifest_tcga_open_expression.txt`)
        - NOTE: Why do we have more files than cases here? Do some samples have both tumor and normal tissue expression?
    
    - **Micro RNA Expression**: `Data Category: transcriptome profiling` + `Data Type: miRNA Expression Quantification` + `Data Format: txt`, for the miRNA expression data (11,082 files, 10,250 cases, 557.23MB stored in `data/manifests/gcd_manifest_tcga_open_mirna.txt`)

    - **DNA Methylation**: `Data Category: DNA methylation` + `Data Type: Methylation Beta Value`, for the DNA methylation data (12,527 files, 11,000 cases, 131.93 GB, stored in `data/manifests/gcd_manifest_tcga_open_meth.txt`)

Next, the files were downloaded using the `gdc-client` tool, which was downloaded from https://gdc.cancer.gov/access-data/gdc-data-transfer-tool, and was then run in the terminal with the following command for each manifest:

```bash
gdc-client download -m data/manifests/{file}.txt -d data/raw/{data_type}/
```

After downloading all data to the `data/raw/` directory, the directory structure was organized as follows:

```
data/
├── raw/
│   ├── clinical/
│   │   ├── *.xml
│   ├── cnv/
│   │   ├── *.seg
│   ├── expression/
│   │   ├── *.htseq.counts
│   ├── mirna/
│   │   ├── *.mirna.quantification.txt
│   ├── meth/
│   │   ├── *.methylation.quantification.txt
```

## Processing Clinical Data

The downloaded clinical data is in XML format, with one file for each sample. Using the `process_clinical.py` script, the data was combined into a single CSV file, stored in `data/processed/clinical.csv`.
