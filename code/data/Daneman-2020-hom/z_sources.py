# RA, 2021-07-19

"""
Single cell transcriptional profile of mouse spinal cord in EAE
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135185

from

CNS fibroblasts form a fibrotic scar in response to immune cell infiltration.
Dorrier CE, Aran D, Haenelt EA et al.
Nat Neurosci 24, 234--244 (2021).
https://doi.org/10.1038/s41593-020-00770-9
"""

import pandas
import tarfile
from pathlib import Path
from tcga.utils import download, unlist1, mkdir

URLS = {
    'data': "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE135185&format=file",
}

src_dir = mkdir(Path(__file__).with_suffix(''))
download = download.to(abs_path=(src_dir / "download_cache"))


def get_expr_from_gse():
    with download(URLS['data']).now.open(mode='rb') as fd:
        with tarfile.TarFile(fileobj=fd, mode='r') as tf:
            with tf.extractfile(unlist1(tf.getmembers())) as tar:
                return pandas.read_table(tar, compression='gzip', index_col=0).astype(int)


if __name__ == '__main__':
    df_expr = get_expr_from_gse()
    df_expr.to_csv(src_dir / "expr.csv.gz", compression='gzip', sep='\t')

df_expr = pandas.read_table(src_dir / "expr.csv.gz", sep='\t', index_col=0).astype(int)
