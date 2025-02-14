import pandas as pd
import gzip
import shutil
import requests



url = "https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPICv2/EPICv2.hg38.manifest.gencode.v41.tsv.gz"

def download_methylation_data(url):

    response = requests.get(url, stream=True)
    with open('methylation_data.tsv.gz', 'wb') as f:
        f.write(response.content)

    with gzip.open('methylation_data.tsv.gz', 'rb') as f_in:
        with open('methylation_data.tsv', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    df = pd.read_csv('methylation_data.tsv', sep='\t')
    return df

df = download_methylation_data(url)

print(df)