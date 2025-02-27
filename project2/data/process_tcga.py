import os
import pickle
import pandas as pd
import gzip
import shutil
import requests

import pandas as pd
import numpy as np
import xml.etree.ElementTree as ET

from collections import defaultdict
from tqdm import tqdm



def load_samplesheet(filename: str) -> pd.DataFrame:
    """
    Loads the sample sheet from the given filename, 
        which is necessary to link the data to the sample metadata.
    This linking is done through the "Case ID" column, 
        which corresponds to the bcr_patient_barcode in the clinical data.
    The file should be downloaded from the GDC portal 
        and correspond exactly to the manifest file.
        
    Parameters:
    - filename (str): The filename of the sample sheet.

    Returns:
    - pd.DataFrame: The sample sheet with the columns:
        ["File Name", "Case ID", "Sample ID", "Sample Type"].
    """

    sample_sheet = pd.read_csv(filename, sep="\t")
    return sample_sheet[["File Name", "Case ID", "Sample ID", "Sample Type"]]


def get_files_with_ext(folder: str, extension: str, exclude=[]) -> list[tuple]:
    """
    Gets all files from the given folder with the given extension,
        excluding the ones that contain any of the given substrings.

    Parameters:
    - folder (str): The folder to search in.
    - extension (str): The extension of the files to search for.
    - exclude (list[str]): The substrings that the files should not contain.

    Returns:
    - list[tuple]: A list of tuples with the uuid and filename of the files.
    """

    # First we collect all the files that are there
    uuids = os.listdir(f"project2/data/raw/{folder}")
    files = []

    for dir in uuids:
        dir_files = os.listdir(f"project2/data/raw/{folder}/" + dir)

        # Find the last file with the given extension, or print an error if no such file exists.
        file_with_extension = None
        for file in dir_files:
            if file.endswith(f".{extension}") and not any([substr in file for substr in exclude]):
                file_with_extension = file
        if file_with_extension is None:
            print(f"ERROR: no .{extension.upper()} file in data dir {dir}")

        files.append((dir, file_with_extension))
    return files


def combine_files_cnv(folder: str, extension: str, exclude=[], sample_sheet=None):

    files = get_files_with_ext(folder, extension, exclude)

    rows = []

    gene_id = None  # Needed as index

    for uuid, file in tqdm(files):
        # Get the sample information from the sample sheet
        case_id, sample_id, sample_type = \
            sample_sheet[sample_sheet["File Name"] == file][["Case ID", "Sample ID", "Sample Type"]].values[0]

        case_id = case_id.split(", ")[0]  # All files have 2 samples, one for tumor and one for healthy tissue.
        # The case ID is the same for both, so we only take the first one.

        # Load the data
        df = pd.read_csv(f"project2/data/raw/{folder}/{uuid}/{file}", sep="\t")
        df.set_index("gene_id", inplace=True)

        # Drop all duplicate gene_names
        df = df.drop_duplicates(subset=['gene_name'])

        # Explicitly check that gene_id is in the same order for all files
        if gene_id is None:
            gene_id = df.index
        else:
            assert gene_id.equals(df.index), f"Gene ID mismatch: {gene_id} vs {df.index}"

        rows.append((case_id, sample_id, sample_type, df["copy_number"].values))

    x = np.array([rows[i][3] for i in range(len(rows))])

    # make a df, with the expression values as columns and the rows as samples
    # the index is the sample_id, the columns are the gene_id
    df = pd.DataFrame(x, columns=df.index, index=[row[1] for row in rows])
    df.index.name = "sample_id"
    df.insert(0, "patient_id", [row[0] for row in rows])  # Equivalent to bcr_patient_barcode in clinical data
    df.insert(1, "sample_type", [row[2] for row in rows])    

    return df


def combine_files_expression(folder: str, files: list[tuple], sample_sheet: pd.DataFrame):

    rows = []

    for uuid, file in tqdm(files):

        case_id, sample_id, sample_type = \
            sample_sheet[sample_sheet["File Name"] == file][["Case ID", "Sample ID", "Sample Type"]].values[0]

        df = pd.read_csv(f"project2/data/raw/{folder}/{uuid}/{file}", sep="\t", comment="#")

        # skip first 4 rows, they are not needed
        df = df.iloc[4:]
        df.reset_index(inplace=True)

        # drop all rows the gene_name is not unique
        df = df.drop_duplicates(subset=['gene_name'])
            
        rows.append((case_id, sample_id, sample_type, df["tpm_unstranded"].values))  # (TPM = transcripts per million)

    # Combine the expression values into a single matrix, and log-normalize with log2(x + 1)
    x = np.array([
        np.log10(np.array(row[3]) + 1)
        for row in rows
    ])

    # make a df, with the expression values as columns and the rows as samples
    # the index is the sample_id, the columns are the gene_name's
    df = pd.DataFrame(x, columns=df.index, index=[row[1] for row in rows])
    df.index.name = "sample_id"
    df.insert(0, "patient_id", [row[0] for row in rows])  # Equivalent to bcr_patient_barcode in clinical data
    df.insert(1, "sample_type", [row[2] for row in rows])

    return df


def get_methylation_annotations():

    if not os.path.exists('project2/data/meth_probe_annotations.tsv'):
        # download the file if it doesn't exist already
        URL = "https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPICv2/EPICv2.hg38.manifest.gencode.v41.tsv.gz"
        response = requests.get(URL, stream=True)
        with open('project2/data/meth_probe_annotations.tsv.gz', 'wb') as f:
            f.write(response.content)
        # unzip the file
        with gzip.open('project2/data/meth_probe_annotations.tsv.gz', 'rb') as f_in:
            with open('project2/data/meth_probe_annotations.tsv', 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    df = pd.read_csv('project2/data/meth_probe_annotations.tsv', sep='\t')

    # remove trailing characters from probeID
    df["probeID"] = df["probeID"].map(lambda x: x.split("_")[0])

    # drop probes that don't have any TSS' associated with them
    df = df[df["distToTSS"].notna()]

    # split the columns for multiple genes, and store them in a long-form dataframe
    long_rows = []
    for i, row in tqdm(df.iterrows(), total=len(df)):
        genes = row["geneNames"].split(";")
        distances = row["distToTSS"].split(";")
        types = row["transcriptTypes"].split(";")
        others = [row["CpG_chrm"], row["probeID"]]

        for gene, distance, type in zip(genes, distances, types):
            long_rows.append([gene, distance, type] + others)
    df_long = pd.DataFrame(long_rows, columns=["geneName", "distToTSS", "transcriptType", "CpG_chrm", "probeID"])
    df_long["distToTSS"] = df_long["distToTSS"].astype(int)

    return df_long


def combine_files_methylation(folder: str, files: list[tuple], sample_sheet: pd.DataFrame, df_map: pd.DataFrame, max_dist: int = 2000):

    # first, we filter the df_map to only include the probes that are within max_dist of a TSS
    df_map["distToTSS"] = df_map["distToTSS"].astype(int)
    df_map = df_map[df_map["distToTSS"] <= max_dist]
    df_map = df_map[df_map["distToTSS"] >= 0]  # If it's after the TSS, remove it.

    map_probe_to_gene = df_map.set_index("probeID")["geneName"].to_dict()

    # then we make a list of all the unique genes, which we will aggregate the probes over
    genes = df_map["geneName"].unique()
    print(f"Found {len(genes)} unique genes within {max_dist} bp of a TSS.")

    rows = []
    for uuid, file in tqdm(files):

        case_id, sample_id, sample_type = \
            sample_sheet[sample_sheet["File Name"] == file][["Case ID", "Sample ID", "Sample Type"]].values[0]
        
        # there are no headers in the files
        df = pd.read_csv(f"project2/data/raw/{folder}/{uuid}/{file}", sep="\t", header=None, index_col=0)
        # df is a single column, with index the probe id, and values the beta values.
        # we use df_map to map the probe IDs to gene IDs

        # first we get only the probes that are in range of a TSS
        df_filtered = df[df.index.isin(set(df_map["probeID"].values))].copy()

        # next, we add the gene names to the df_filtered
        df_filtered["geneName"] = df_filtered.index.to_series().map(map_probe_to_gene)

        # presort by geneName, so we can use groupby
        df_filtered.sort_values("geneName", inplace=True)

        df_accum = df_filtered.groupby("geneName").sum()

        # and now we have to arrange the genes in the order given by the "genes" array,
        # and fill in the missing ones with 0
        gene_accum_meths = df_accum.reindex(genes).fillna(0).values

        rows.append((case_id, sample_id, sample_type, gene_accum_meths))

    x = np.array([
        np.array(row[3])
        for row in rows
    ]).reshape(len(rows), len(genes))

    # make a df, with the expression values as columns and the rows as samples
    # the index is the sample_id, the columns are the gene_id
    df = pd.DataFrame(x, columns=genes, index=[row[1] for row in rows])
    df.index.name = "sample_id"
    df.insert(0, "patient_id", [row[0] for row in rows])  # Equivalent to bcr_patient_barcode in clinical data
    df.insert(1, "sample_type", [row[2] for row in rows])

    return df


def parse_clinical_list(tag: str, column: ET.Element) -> tuple:
    """
    Parse function for if the clinical data field is a list-type structure.
    """
    values = [v.text for v in column]
    if values == [None]:
        return tag, pd.NA
    else:
        return tag, ", ".join(values)
    

def combine_files_clinical(clinical_files: list[tuple]):

    patients = {}
    missing_col_counts = defaultdict(int)

    for uuid, file in tqdm(clinical_files):
        tree = ET.parse("project2/data/raw/clinical/" + uuid + "/" + file)
        patient = tree.getroot()[-1]
        assert patient.tag.endswith("patient")

        patient_dict = {}
        for column in patient:

            # By default, we just use the tag and value as the key/value-pair.
            tag = column.tag.split("}")[-1]
            value = column.text

            # Some values are lists, so we need to handle them differently
            if tag in ['race_list',  # Note: apparently there are never multiple in a "list"
                    'metastatic_site_list', 
                    'relation_testicular_cancer_list', 
                    'postoperative_tx_list'
                    ]:
                tag, value = parse_clinical_list(tag, column)

            if tag == "blood_relative_cancer_history_list":
                relatives = []
                for relative in column:
                    d = {}
                    for c in column[0]:
                        subtag = c.tag.split("}")[-1]
                        subvalue = c.text
                        d[subtag] = subvalue

                    type_col = "cancer_diagnosis_cancer_type_icd9_text_name"
                    if not type_col in d.keys():
                        type_col = "relative_family_cancer_hx_text"
                    if not type_col in d.keys():
                        type_col = "family_history_cancer_type_other"
                    if not type_col in d.keys():
                        type_col = "family_cancer_type_txt"

                    try:
                        comb = f'{d["family_member_relationship_type"]}:{d[type_col]}'
                    except:
                        print(f"Unable to parse family history: {d}")
                    if comb is not None and comb != "None:None" and comb not in relatives:
                        relatives.append(comb)

                if len(relatives) > 1:
                    print(relatives)

                if relatives == []:
                    value = pd.NA
                else:
                    value = ", ".join(relatives)

            # TODO: other columns that have a weird format, but could be included as well are:
            missing_columns = ["stage event", "new_tumor_events", "drugs", "radiations",
                "follow_ups", "history_hepato_carcinoma_risk_factors",
                "loss_expression_of_mismatch_repair_proteins_by_ihc_results",
                "antireflux_treatment_types", "sites_of_primary_melanomas",
                "viral_hepatitis_serologies", "prior_systemic_therapy_types",
                "anatomic_neoplasm_subdivisions", "first_nonlymph_node_metastasis_anatomic_sites",
                "patient_history_immune_system_and_related_disorders_names",
                "lymph_node_location_positive_pathology_names", "fdg_or_ct_pet_performed_outcomes",
                "diagnostic_mri_result_outcomes", "diagnostic_ct_result_outcomes",
                "human_papillomavirus_types","treatment", 
                "relation_testicular_cancer_list", "postoperative_tx_list"]

            if tag in missing_columns:
                nonempty = False
                try:
                    for c in column:
                        nonempty = True
                except:
                    pass

                if nonempty:
                    missing_col_counts[tag] += 1
                value = pd.NA

            # Sometimes there's a random newline in the value
            if value is not None and isinstance(value, str):
                value = value.strip()

            patient_dict[tag] = value
        
        patients[patient_dict["patient_id"]] = patient_dict

    print("\nThese columns are not correctly parsed because they need manual intervention\n" + \
        "     (amount of samples that have data here between brackets):")
    for col, count in missing_col_counts.items():
        print(f"- {col}: {count}")

    # Make it into a dataframe
    df = pd.DataFrame(patients).T
    return df


def combine_files_mirna(folder: str, files: list[tuple], sample_sheet: pd.DataFrame):

    rows = []

    for uuid, file in tqdm(files):

        case_id, sample_id, sample_type = \
            sample_sheet[sample_sheet["File Name"] == file][["Case ID", "Sample ID", "Sample Type"]].values[0]

        df = pd.read_csv(f"project2/data/raw/{folder}/{uuid}/{file}", sep="\t")

        rows.append((case_id, sample_id, sample_type, df["reads_per_million_miRNA_mapped"].values))

    x = np.array([row[3] for row in rows])

    gene_ids = df["miRNA_ID"]

    df = pd.DataFrame(x, columns=gene_ids, index=[row[1] for row in rows])
    df.index.name = "sample_id"
    df.insert(0, "patient_id", [row[0] for row in rows])
    df.insert(1, "sample_type", [row[2] for row in rows])

    return df


DATA_DIR = "project2/data"


if __name__ == "__main__":

    # TODO: download data programatically through the GDC API, 
    #   instead of manually collecting manifests and downloading those


    # Clinical data
    out_file_clinical = f"{DATA_DIR}/processed/metadata.csv"
    # if os.path.exists(out_file_clinical):
    #     print("Clinical data already processed, loading...")
    #     df_clinical = pd.read_csv(out_file_clinical, low_memory=False)
    # else:
    print("Processing clinical data...")
    clinical_files = get_files_with_ext(folder="clinical", extension="xml")
    df_clinical = combine_files_clinical(clinical_files=clinical_files)
    print("Saving clinical data...")
    df_clinical.to_csv(f"{DATA_DIR}/processed/metadata.csv", index=False)
    print(f"Done. Loaded clinical data for {len(df_clinical)} patients.")
    df_clinical_index = df_clinical["bcr_patient_barcode"].copy().values


    # Expression
    out_file_expression = f"{DATA_DIR}/processed/expression.pkl"
    # if os.path.exists(out_file_expression):
    #     print("Expression data already processed, loading...")
    #     with open(out_file_expression, "rb") as f:
    #         df_ge = pickle.load(f)
    # else:
    print("\nProcessing expression data...")
    samplesheet_expression = load_samplesheet(filename=f"{DATA_DIR}/manifests/gdc_sample_sheet_tcga_open_expression.tsv")
    files_expression = get_files_with_ext(folder="expression", extension="tsv")
    df_ge = combine_files_expression(folder="expression", files=files_expression, sample_sheet=samplesheet_expression)
    print("Saving expression data...")
    with open(out_file_expression, "wb") as f:
        pickle.dump(df_ge, f)  
    print(f"Done. Loaded expression data for {len(df_ge)} samples.")
    df_ge_samples = df_ge["patient_id"].unique()
    del df_ge  # Free up memory
        
    exit(0)
    
    # CNV
    out_file_cnv = f"{DATA_DIR}/processed/cnv.pkl"
    if os.path.exists(out_file_cnv):
        print("CNV data already processed, loading...")
        with open(out_file_cnv, "rb") as f:
            df_cnv = pickle.load(f)
    else:
        print("\nProcessing CNV data...")
        samplesheet_cnv = load_samplesheet(filename=f"{DATA_DIR}/manifests/gdc_sample_sheet_tcga_open_gene-level-cn_ascat3.tsv")
        df_cnv = combine_files_cnv(folder="cnv", extension="tsv", exclude=["annotation"], sample_sheet=samplesheet_cnv)
        print("Saving CNV data...")
        with open(out_file_cnv, "wb") as f:
            pickle.dump(df_cnv, f)
    print(f"Done. Loaded CNV data for {len(df_cnv)} samples.")
    df_cnv_samples = df_cnv["patient_id"].unique()
    del df_cnv  # Free up memory

    # Methylation
    print("\nProcessing methylation data...")
    if os.path.exists(f"{DATA_DIR}/processed/meth.pkl"):
        print("Methylation data already processed, loading...")
        with open(f"{DATA_DIR}/processed/meth.pkl", "rb") as f:
            df_meth = pickle.load(f)

    else:
        # Load (or download if necessary) the methylation annotations
        if not os.path.exists(f"{DATA_DIR}/meth_probe_to_TSS_map.tsv"):
            df_probe_ann = get_methylation_annotations()
            df_probe_ann.to_csv(f"{DATA_DIR}/meth_probe_to_TSS_map.tsv", sep="\t")
        else:
            df_probe_ann = pd.read_csv(f"{DATA_DIR}/meth_probe_to_TSS_map.tsv", sep="\t")

        samplesheet_meth = load_samplesheet(filename=f"{DATA_DIR}/manifests/gdc_sample_sheet_tcga_open_meth.tsv")
        files_meth = get_files_with_ext(folder="meth", extension="txt", exclude=["annotations"])
        df_meth = combine_files_methylation(folder="meth", files=files_meth, sample_sheet=samplesheet_meth, df_map=df_probe_ann)

        # check how many columns sum to 0
        print("Number of genes that are always 0:")
        col_sums = [df_meth[c].sum() for c in df_meth.columns if c not in ["patient_id", "sample_type"]]
        print(sum([s == 0 for s in col_sums]))  # Output: 6120

        print("Saving methylation data...")
        with open(f"{DATA_DIR}/processed/meth.pkl", "wb") as f:
            pickle.dump(df_meth, f)

    df_meth_samples = df_meth["patient_id"].copy()
    del df_meth  # Free up memory

    # miRNA
    if os.path.exists(f"{DATA_DIR}/processed/mirna.pkl"):
        print("miRNA data already processed, loading...")
        with open(f"{DATA_DIR}/processed/mirna.pkl", "rb") as f:
            df_mirna = pickle.load(f)

    else:
        samplesheet_mirna = load_samplesheet(filename=f"{DATA_DIR}/manifests/gdc_sample_sheet_tcga_open_mirna.tsv")
        files_mirna = get_files_with_ext(folder="mirna", extension="txt", exclude=["annotations"])
        df_mirna = combine_files_mirna(folder="mirna", files=files_mirna, sample_sheet=samplesheet_mirna)

        print("Saving miRNA data...")
        with open(f"{DATA_DIR}/processed/mirna.pkl", "wb") as f:
            pickle.dump(df_mirna, f)
        df_mirna.to_csv(f"{DATA_DIR}/processed/mirna.csv")

    # get the samples
    df_mirna_samples = df_mirna["patient_id"].copy()
    del df_mirna  # Free up memory

    # Now we get for each data type the samples that are present in all data types
    samples_list = [df_ge_samples, df_cnv_samples, df_meth_samples, df_mirna_samples, df_clinical_index]
    samples = set(samples_list[0])
    for s in samples_list[1:]:
        samples = samples.intersection(set(s))
    print(f"\nFound {len(samples)} samples that are present in all data types.")

    # Now we filter the dataframes to only include these samples
    df_ge = pd.read_pickle(out_file_expression)
    df_ge = df_ge[df_ge["patient_id"].isin(samples)]
    df_ge.to_pickle(out_file_expression.replace(".pkl", "_overlap.pkl"))
    del df_ge

    df_cnv = pd.read_pickle(out_file_cnv)
    df_cnv = df_cnv[df_cnv["patient_id"].isin(samples)]
    df_cnv.to_pickle(out_file_cnv.replace(".pkl", "_overlap.pkl"))
    del df_cnv

    df_meth = pd.read_pickle(f"{DATA_DIR}/processed/meth.pkl")
    df_meth = df_meth[df_meth["patient_id"].isin(samples)]
    df_meth.to_pickle(f"{DATA_DIR}/processed/meth_overlap.pkl")
    del df_meth

    df_mirna = pd.read_pickle(f"{DATA_DIR}/processed/mirna.pkl")
    df_mirna = df_mirna[df_mirna["patient_id"].isin(samples)]
    df_mirna.to_pickle(f"{DATA_DIR}/processed/mirna_overlap.pkl")
    del df_mirna

    df_clinical = pd.read_csv(out_file_clinical, low_memory=False)
    df_clinical = df_clinical[df_clinical["bcr_patient_barcode"].isin(samples)]
    df_clinical.to_csv(out_file_clinical.replace(".csv", "_overlap.csv"), index=False)
    del df_clinical