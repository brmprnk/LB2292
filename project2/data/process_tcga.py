import os
import pickle
import gzip
import shutil
import requests

import pandas as pd
import numpy as np

from tqdm import tqdm

from project2.data.process_clinical import parse_clinical
from project2.data.process_utils import find_data_files, load_samplesheet, from_sample_sheet


def process_cnv(folder: str, sample_sheet: pd.DataFrame):
    """
    Process the CNV data.
    
    Args:
    - folder (str): The path to the CNV data.
    - sample_sheet (pd.DataFrame): The sample sheet to get the case_id, sample_id and sample_type from.

    Returns:
    - pd.DataFrame: A dataframe with the
    """
    rows = []
    for file in tqdm(find_data_files(folder, extension="tsv"), desc="Processing CNV data"):
        case_id, sample_id, sample_type = from_sample_sheet(sample_sheet, file)
        case_id = case_id.split(", ")[0]  # All files have 2 samples, one for tumor and one for healthy tissue.
        # The case ID is the same for both, so we only take the first one.

        df = pd.read_csv(file, sep="\t", comment="#").set_index("gene_id")
        df = df.drop_duplicates(subset=['gene_name']).set_index("gene_name")  # Drop all duplicate gene_names
        rows.append((case_id, sample_id, sample_type, df["copy_number"].values))

    x = np.array([rows[i][3] for i in range(len(rows))])

    df = pd.DataFrame(x, columns=df.index, index=[row[1] for row in rows])
    df.index.name = "sample_id"
    df.insert(0, "patient_id", [row[0] for row in rows])  # Equivalent to bcr_patient_barcode in clinical data
    df.insert(1, "sample_type", [row[2] for row in rows])    
    return df


def process_gene_expression(path: str, sample_sheet: pd.DataFrame):
    """
    Process the gene expression data.
    
    Args:
    - path (str): The path to the gene expression data.
    - sample_sheet (pd.DataFrame): The sample sheet to get the case_id, sample_id and sample_type from.

    Returns:
    - pd.DataFrame: A dataframe with the gene expression data.
    """
    rows = []
    for file in tqdm(find_data_files(path, extension="tsv"), desc="Processing gene expression data"):
        case_id, sample_id, sample_type  = from_sample_sheet(sample_sheet, file)
        df = pd.read_csv(file, sep="\t", comment="#").iloc[4:]  # Read, skip first 4 rows, they are not needed
        df = df.drop_duplicates(subset=['gene_name'])  # drop all rows the gene_name is not unique (~1300 genes)
        df.set_index("gene_name", inplace=True)
        rows.append((case_id, sample_id, sample_type, df["unstranded"].values))  # Raw unstranded counts

    # Combine the expression values into a single matrix, normalize to counts per million, and log-normalize
    x = np.array([row[3] for row in rows])
    x = x / x.sum(axis=1)[:, None] * 1e6
    x = np.log10(x + 1)

    df = pd.DataFrame(x, columns=df.index, index=[row[1] for row in rows])
    df.index.name = "sample_id"
    df.insert(0, "patient_id", [row[0] for row in rows])  # Equivalent to bcr_patient_barcode in clinical data
    df.insert(1, "sample_type", [row[2] for row in rows])
    return df


def get_methylation_annotations() -> pd.DataFrame:
    """
    We have to download the methylation annotations, to know which probes are close to a TSS.

    Returns:
    - pd.DataFrame: A dataframe with the columns:
        ["geneName", "distToTSS", "transcriptType", "CpG_chrm", "probeID"]
    """

    if not os.path.exists('project2/data/raw/meth_probe_annotations.tsv'):
        # download the file if it doesn't exist already
        URL = "https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPICv2/EPICv2.hg38.manifest.gencode.v41.tsv.gz"
        response = requests.get(URL, stream=True)
        with open('project2/data/raw/meth_probe_annotations.tsv.gz', 'wb') as f:
            f.write(response.content)
        # unzip the file
        with gzip.open('project2/data/raw/meth_probe_annotations.tsv.gz', 'rb') as f_in:
            with open('project2/data/raw/meth_probe_annotations.tsv', 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    df = pd.read_csv('project2/data/raw/meth_probe_annotations.tsv', sep='\t')

    # remove trailing characters from probeID
    df["probeID"] = df["probeID"].map(lambda x: x.split("_")[0])

    # drop probes that don't have any TSS' associated with them
    df = df[df["distToTSS"].notna()]

    # split the columns for multiple genes, and store them in a long-form dataframe
    long_rows = []
    for _, row in tqdm(df.iterrows(), total=len(df)):
        genes = row["geneNames"].split(";")
        distances = row["distToTSS"].split(";")
        types = row["transcriptTypes"].split(";")
        others = [row["CpG_chrm"], row["probeID"]]

        for gene, distance, type in zip(genes, distances, types):
            long_rows.append([gene, distance, type] + others)
    df_long = pd.DataFrame(long_rows, columns=["geneName", "distToTSS", "transcriptType", "CpG_chrm", "probeID"])
    df_long["distToTSS"] = df_long["distToTSS"].astype(int)

    return df_long


def process_methylation(path: str, sample_sheet: pd.DataFrame, df_map: pd.DataFrame, max_dist: int):
    """
    Process the methylation data. This is a bit more involved, as we have to aggregate the probes over genes.
    To aggregate, we add up the beta values of all probes that are within max_dist of a TSS.
    For this we use the df_map, which maps the probe IDs to gene names, and it's downloaded by the get_methylation_annotations function.
    
    Args:
    - path (str): The path to the methylation data.
    - sample_sheet (pd.DataFrame): The sample sheet to get the case_id, sample_id and sample_type from.

    Returns:
    - pd.DataFrame: A dataframe with the methylation data.

    """
    # first, we filter the df_map to only include the probes that are within max_dist of a TSS
    df_map["distToTSS"] = df_map["distToTSS"].astype(int)
    df_map = df_map[df_map["distToTSS"] <= max_dist]
    df_map = df_map[df_map["distToTSS"] >= 0]  # If it's after the TSS, remove it.

    map_probe_to_gene = df_map.set_index("probeID")["geneName"].to_dict()

    # then we make a list of all the unique genes, which we will aggregate the probes over
    genes = df_map["geneName"].unique()
    print(f"Found {len(genes)} unique genes within {max_dist} bp of a TSS.")

    rows = []
    for file in tqdm(find_data_files(path, extension="txt")):
        case_id, sample_id, sample_type = from_sample_sheet(sample_sheet, file)
        df = pd.read_csv(file, sep="\t", header=None, index_col=0)  # (no headers in these files)

        # first we get only the probes that are in range of a TSS
        df_filtered = df[df.index.isin(set(df_map["probeID"].values))].copy()
        # next, we add the gene names to the df_filtered
        df_filtered["geneName"] = df_filtered.index.to_series().map(map_probe_to_gene)

        # presort by geneName, so we can group by geneName (and fill in missing genes with 0s)
        df_filtered.sort_values("geneName", inplace=True)
        df_accum = df_filtered.groupby("geneName").sum()
        gene_accum_meths = df_accum.reindex(genes).fillna(0).values

        rows.append((case_id, sample_id, sample_type, gene_accum_meths))

    x = np.array([np.array(row[3]) for row in rows]).reshape(len(rows), len(genes))

    df = pd.DataFrame(x, columns=genes, index=[row[1] for row in rows])
    df.index.name = "sample_id"
    df.insert(0, "patient_id", [row[0] for row in rows])  # Equivalent to bcr_patient_barcode in clinical data
    df.insert(1, "sample_type", [row[2] for row in rows])
    return df


def process_mirna_expression(path: str, sample_sheet: pd.DataFrame) -> pd.DataFrame:
    """
    Process the miRNA expression data. We use the reads_per_million_miRNA_mapped column,
    and log-transform the data.

    Args:
    - path (str): The path to the miRNA expression data.
    - sample_sheet (pd.DataFrame): The sample sheet to get the case_id, sample_id and sample_type from.

    Returns:
    - pd.DataFrame: A dataframe with the miRNA expression data.
    """
    rows = []
    for file in tqdm(find_data_files(path, extension="txt")):

        case_id, sample_id, sample_type = \
            sample_sheet[sample_sheet["File Name"] == file][["Case ID", "Sample ID", "Sample Type"]].values[0]

        df = pd.read_csv(file, sep="\t", comment="#")

        rows.append((case_id, sample_id, sample_type, df["reads_per_million_miRNA_mapped"].values))

    x = np.array([row[3] for row in rows])
    # also log-transform the miRNA data
    x = np.log10(x + 1)

    gene_ids = df["miRNA_ID"]

    df = pd.DataFrame(x, columns=gene_ids, index=[row[1] for row in rows])
    df.index.name = "sample_id"
    df.insert(0, "patient_id", [row[0] for row in rows])
    df.insert(1, "sample_type", [row[2] for row in rows])
    return df


if __name__ == "__main__":
    DATA_DIR = "project2/data"
    PATH_CLINICAL = os.path.join(DATA_DIR, "processed/clinical.csv")
    PATH_EXPRESSION = os.path.join(DATA_DIR, "processed/expression.pkl")
    PATH_CNV = os.path.join(DATA_DIR, "processed/cnv.pkl")
    PATH_METH = os.path.join(DATA_DIR, "processed/meth.pkl")
    PATH_MIRNA = os.path.join(DATA_DIR, "processed/mirna.pkl")

    # TODO: download data programatically through the GDC API, 
    #   instead of manually collecting manifests and downloading those


    # Clinical data
    df_clinical = parse_clinical(os.path.join(DATA_DIR, "raw/clinical"))
    df_clinical.to_csv(PATH_CLINICAL.replace(".csv", "_full.csv"), index=False)
    print(f"Done. Loaded clinical data for {len(df_clinical)} patients.\n")


    # Expression
    df_ge = process_gene_expression(
        path=os.path.join(DATA_DIR, "raw/expression"), 
        sample_sheet=load_samplesheet(filename=f"{DATA_DIR}/manifests/gdc_sample_sheet_tcga_open_expression.tsv")
    )
    with open(PATH_EXPRESSION.replace(".pkl", "_full.pkl"), "wb") as f:
        pickle.dump(df_ge, f)  
    df_ge_samples = df_ge["patient_id"].unique()
    del df_ge  # Free up memory
    print(f"Done. Loaded expression data for {len(df_ge)} samples.")
    

    # CNV
    df_cnv = process_cnv(
        path=os.path.join(f"{DATA_DIR}/raw/cnv"),
        sample_sheet=load_samplesheet(filename=f"{DATA_DIR}/manifests/gdc_sample_sheet_tcga_open_gene-level-cn_ascat3.tsv")
    )
    with open(PATH_CNV.replace(".pkl", "_full.pkl"), "wb") as f:
        pickle.dump(df_cnv, f)
    df_cnv_samples = df_cnv["patient_id"].unique()
    del df_cnv  # Free up memory
    print(f"Done. Loaded CNV data for {len(df_cnv)} samples.")


    # Methylation
    ANNOTATIONS_PATH = f"{DATA_DIR}/meth_probe_to_TSS_map.tsv"
    if not os.path.exists(ANNOTATIONS_PATH):
        get_methylation_annotations().to_csv(ANNOTATIONS_PATH, sep="\t")
    df_meth = process_methylation(
        path=f"{DATA_DIR}/raw/meth", 
        sample_sheet=load_samplesheet(filename=f"{DATA_DIR}/manifests/gdc_sample_sheet_tcga_open_meth.tsv"), 
        df_map=pd.read_csv(ANNOTATIONS_PATH, sep="\t"),
        max_dist=2000,
    )
    with open(PATH_METH.replace(".pkl", "_full.pkl"), "wb") as f:
        pickle.dump(df_meth, f)
    df_meth_samples = df_meth["patient_id"].copy()
    del df_meth  # Free up memory


    # miRNA
    df_mirna = process_mirna_expression(
        path=f"{DATA_DIR}/raw/mirna", 
        sample_sheet=load_samplesheet(filename=f"{DATA_DIR}/manifests/gdc_sample_sheet_tcga_open_mirna.tsv"),
    )
    with open(PATH_MIRNA.replace(".pkl", "_full.pkl"), "wb") as f:
        pickle.dump(df_mirna, f)
    df_mirna_samples = df_mirna["patient_id"].copy()
    del df_mirna  # Free up memory


    # Now we get for each data type the samples that are present in all data types
    samples_list = [df_ge_samples, df_cnv_samples, df_meth_samples, df_mirna_samples, df_clinical.index]
    samples = set(samples_list[0])
    for s in samples_list[1:]:
        samples = samples.intersection(set(s))
    print(f"\nFound {len(samples)} samples that are present in all data types.")

    # Now we filter the dataframes to only include these samples
    df_clinical = df_clinical[df_clinical["bcr_patient_barcode"].isin(samples)]
    df_clinical.to_csv(PATH_CLINICAL, index=False)
    df_ge = pd.read_pickle(PATH_EXPRESSION.replace(".pkl", "_full.pkl"))
    df_ge = df_ge[df_ge["patient_id"].isin(samples)]
    df_ge.to_pickle(PATH_EXPRESSION)
    del df_ge
    df_cnv = pd.read_pickle(PATH_CNV.replace(".pkl", "_full.pkl"))
    df_cnv = df_cnv[df_cnv["patient_id"].isin(samples)]
    df_cnv.to_pickle(PATH_CNV)
    del df_cnv
    df_meth = pd.read_pickle(PATH_METH.replace(".pkl", "_full.pkl"))
    df_meth = df_meth[df_meth["patient_id"].isin(samples)]
    df_meth.to_pickle(PATH_METH)
    del df_meth
    df_mirna = pd.read_pickle(PATH_MIRNA.replace(".pkl", "_full.pkl"))
    df_mirna = df_mirna[df_mirna["patient_id"].isin(samples)]
    df_mirna.to_pickle(PATH_MIRNA)
    del df_mirna
