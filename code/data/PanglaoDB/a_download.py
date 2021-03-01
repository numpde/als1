# RA, 2020-12-10

"""
Cell type gene expression markers from PanglaoDB
https://panglaodb.se/markers.html?cell_type=%27choose%27

Oscar Franzén, Li-Ming Gan, Johan L M Björkegren
PanglaoDB: a web server for exploration of mouse and human
single-cell RNA sequencing data,
Database, 2019
http://dx.doi.org/10.1093/database/baz046
"""

URL = "https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz"

import pandas as pd
from pathlib import Path
from tcga.utils import download, pretty

download = download.to(abs_path=Path(__file__).with_suffix(''))
data = download(URL).now

print(data.meta)

with data.open(mode='rb') as fd:
    df = pd.read_csv(fd, compression='gzip', sep='\t')

celltypes = set(df['cell type'])
print('Meningeal cells' in celltypes)

print(df[df['cell type'] == 'Meningeal cells'].to_markdown())

# print(df['cell type'].value_counts().to_markdown())
