import os
import pickle

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

        # Explicitly check that gene_id is in the same order for all files
        if gene_id is None:
            gene_id = df.index
        else:
            assert gene_id.equals(df.index), f"Gene ID mismatch: {gene_id} vs {df.index}"

        rows.append((case_id, sample_id, sample_type, df["copy_number"].values))

    x = np.array([rows[i][3] for i in range(len(rows))])

    # make a df, with the expression values as columns and the rows as samples
    # the index is the sample_id, the columns are the gene_id
    df = pd.DataFrame(x, columns=gene_id, index=[row[1] for row in rows])
    df.index.name = "sample_id"
    df.insert(0, "patient_id", [row[0] for row in rows])  # Equivalent to bcr_patient_barcode in clinical data
    df.insert(1, "sample_type", [row[2] for row in rows])    

    return df

def combine_files_expression(folder: str, files: list[tuple], sample_sheet: pd.DataFrame):

    rows = []

    for uuid, file in tqdm(files):

        # case_id = sample_sheet[sample_sheet["File Name"] == file]["Case ID"].values[0]
        case_id, sample_id, sample_type = \
            sample_sheet[sample_sheet["File Name"] == file][["Case ID", "Sample ID", "Sample Type"]].values[0]

        # skip rows with # in the beginning
        df = pd.read_csv(f"project2/data/raw/{folder}/{uuid}/{file}", sep="\t", comment="#")

        # skip first 4 rows, they are not needed
        df = df.iloc[4:]
        df.reset_index(inplace=True)
            
        rows.append((case_id, sample_id, sample_type, df["tpm_unstranded"].values))

    x = np.array([
        np.array(row[3])
        for row in rows
    ])

    gene_ids, gene_names = df["gene_id"], df["gene_name"]
    # TODO: We want to use the gene_name for columns, and not the gene_id. 
    # However, the gene_name is not unique. Options:
    # - add a suffix to the gene_name to make it unique
    # - aggregate the values of genes with the same gene_name

    # make a df, with the expression values as columns and the rows as samples
    # the index is the sample_id, the columns are the gene_id
    df = pd.DataFrame(x, columns=gene_ids, index=[row[1] for row in rows])
    df.index.name = "sample_id"
    df.insert(0, "patient_id", [row[0] for row in rows])  # Equivalent to bcr_patient_barcode in clinical data
    df.insert(1, "sample_type", [row[2] for row in rows])

    return df


def combine_files_methylation(folder: str, files: list[tuple], sample_sheet: pd.DataFrame):

    rows = []

    for uuid, file in tqdm(files):

        # case_id = sample_sheet[sample_sheet["File Name"] == file]["Case ID"].values[0]
        try:
            case_id, sample_id, sample_type = \
                sample_sheet[sample_sheet["File Name"] == file][["Case ID", "Sample ID", "Sample Type"]].values[0]
        except:
            print(f"Error with file {file}, uuid {uuid}")
            exit(-1)

        # there are no headers in the files
        df = pd.read_csv(f"project2/data/raw/{folder}/{uuid}/{file}", sep="\t", header=None, index_col=0)

        # unfortunately the indices are not the same for all files
        # so we store both the values and the indices

        rows.append((case_id, sample_id, sample_type, df[1]))

    # Now we collect all the unique indices, and sort them
    indices = set()
    for _, _, _, row in rows:
        indices.update(row.index)
    indices = sorted(list(indices))

    # Now we create the matrix
    x = np.zeros((len(rows), len(indices)))
    for i, (_, _, _, row) in enumerate(tqdm(rows)):
        x[i, :] = row.reindex(indices).values

    # make a df, with the expression values as columns and the rows as samples
    # the index is the sample_id, the columns are the gene_id
    df = pd.DataFrame(x, columns=indices, index=[row[1] for row in rows])
    df.index.name = "sample_id"
    df.insert(0, "patient_id", [row[0] for row in rows])  # Equivalent to bcr_patient_barcode in clinical data
    df.insert(1, "sample_type", [row[2] for row in rows])

    return df


def combine_files_clinical(clinical_files: list[tuple]):

    patients = {}
    missing_col_counts = defaultdict(int)

    for uuid, file in tqdm(clinical_files):
        tree = ET.parse("project2/data/raw/clinical/" + uuid + "/" + file)

        patient = tree.getroot()[-1]
        assert patient.tag.endswith("patient")

        patient_dict = {}

        for column in patient:
            tag = column.tag.split("}")[-1]
            value = column.text

            # Some values are lists, so we need to handle them differently
            if tag in ['race_list', 
                    'metastatic_site_list', 'relation_testicular_cancer_list', 
                    'postoperative_tx_list']:
                values = []
                for v in column:
                    values.append(v.text)

                if values == [None]:
                    value = pd.NA
                else:
                    value = ", ".join(values)

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


DATA_DIR = "project2/data"


if __name__ == "__main__":

    # TODO: download data programatically through the GDC API, 
    #   instead of manually collecting manifests and downloading those


    # Clinical data
    out_file_clinical = f"{DATA_DIR}/processed/metadata.csv"
    if os.path.exists(out_file_clinical):
        print("Clinical data already processed, loading...")
        df_clinical = pd.read_csv(out_file_clinical, low_memory=False)
    else:
        print("Processing clinical data...")
        clinical_files = get_files_with_ext(folder="clinical", extension="xml")
        df_clinical = combine_files_clinical(clinical_files=clinical_files)
        print("Saving clinical data...")
        df_clinical.to_csv(f"{DATA_DIR}/processed/metadata.csv", index=False)
    print(f"Done. Loaded clinical data for {len(df_clinical)} patients.")
    df_clinical_index = df_clinical["bcr_patient_barcode"].copy().values

    # Expression
    out_file_expression = f"{DATA_DIR}/processed/expression.pkl"
    if os.path.exists(out_file_expression):
        print("Expression data already processed, loading...")
        with open(out_file_expression, "rb") as f:
            df_ge = pickle.load(f)
    else:
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
    samplesheet_meth = load_samplesheet(filename=f"{DATA_DIR}/manifests/gdc_sample_sheet_tcga_open_meth.tsv")
    files_meth = get_files_with_ext(folder="meth", extension="txt", exclude=["annotations"])
    df_meth = combine_files_methylation(folder="meth", files=files_meth, sample_sheet=samplesheet_meth)
    print("Saving methylation data...")
    with open(f"{DATA_DIR}/processed/meth.pkl", "wb") as f:
        pickle.dump(df_meth, f)
    df_meth_samples = df_meth["patient_id"].copy()
    del df_meth  # Free up memory
    print(df_meth_samples)

    print(f"Done. Loaded methylation data for {len(df_meth_samples)} samples.")
    overlap = set(df_meth_samples).intersection(set(df_clinical_index))
    print(f"And metadata is available for {len(overlap)}/{len(df_meth_samples)} of these samples.")
    overlap_ge_meth = set(df_ge_samples).intersection(set(df_meth_samples))
    print(f"Also, the expression data has {len(overlap_ge_meth)} samples in common with the methylation data.")
    overlap_cnv_meth = set(df_cnv_samples).intersection(set(df_meth_samples))
    print(f"Also, the CNV data has {len(overlap_cnv_meth)} samples in common with the methylation data.")
    overlap_all = set(df_ge_samples).intersection(set(df_cnv_samples)).intersection(set(df_meth_samples)).intersection(set(df_clinical_index))
    print(f"Finally, there are {len(overlap_all)} samples that have data for all three data types, and metadata available.")


    # TODO: miRNA
    print("\nTODO: miRNA data...")

