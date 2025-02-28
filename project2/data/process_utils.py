import os

import pandas as pd


def find_data_files(path: str, ext: str) -> list[str]:
    """
    Find all files with the given extension `ext` in the given path.
    We also exclude the annotations.txt files here.

    Parameters:
    - path (str): The path to search in.
    - ext (str): The extension of the files to search for.

    Returns:
    - list[str]: A list of the full paths of the files found.
    """
    files = []
    for uuid in os.listdir(path):
        uuid_file = None
        for file in os.listdir(os.path.join(path, uuid)):
            if file.endswith(ext) and not file == "annotations.txt":
                uuid_file = f"{path}/{uuid}/{file}"
                break
        if uuid_file is None:
            print(f"No file with extension {ext} found in {uuid}.")
        else:
            files.append(uuid_file)
    return files


def load_samplesheet(path: str) -> pd.DataFrame:
    """
    Loads the sample sheet from the given path, 
        which is necessary to link the data to the sample metadata.
    This linking is done through the "Case ID" column, 
        which corresponds to the bcr_patient_barcode in the clinical data.
    The file should be downloaded from the GDC portal 
        and correspond exactly to the manifest file.
        
    Parameters:
    - path (str): The path to the sample sheet.

    Returns:
    - pd.DataFrame: The sample sheet with the columns:
        ["File Name", "Case ID", "Sample ID", "Sample Type"].
    """

    sample_sheet = pd.read_csv(path, sep="\t")
    return sample_sheet[["File Name", "Case ID", "Sample ID", "Sample Type"]]


def from_sample_sheet(sample_sheet: pd.DataFrame, filepath: str) -> tuple:
    """
    Get the case_id, sample_id and sample_type from the sample sheet, based on the filepath.

    Args:
    - sample_sheet (pd.DataFrame): The sample sheet to look in.
    - filepath (str): The filepath to the file to get the information for.

    Returns:
    - tuple: The case_id, sample_id and sample_type.
    """
    return sample_sheet[sample_sheet["File Name"] == filepath.split("/")[-1]][["Case ID", "Sample ID", "Sample Type"]].values[0]
