import pandas as pd
import gzip
import shutil
import requests
import os
from tqdm import tqdm


def get_methylation_annotations():

    URL = "https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPICv2/EPICv2.hg38.manifest.gencode.v41.tsv.gz"

    if not os.path.exists('project2/data/methylation_data.tsv.gz'):
        # download the file if it doesn't exist already
        response = requests.get(URL, stream=True)
        with open('project2/data/methylation_data.tsv.gz', 'wb') as f:
            f.write(response.content)

        # unzip the file
        with gzip.open('project2/data/methylation_data.tsv.gz', 'rb') as f_in:
            with open('project2/data/methylation_data.tsv', 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

    df = pd.read_csv('project2/data/methylation_data.tsv', sep='\t')

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

if __name__ == "__main__":

    df = get_methylation_annotations()

    # save it
    df.to_csv("project2/data/meth_probe_to_TSS_map.tsv", sep="\t")