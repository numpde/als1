# RA, 2021-05-02

"""
Prepare the list of up- and down-regulated genes from the paper:

Identification of fibroblast activation-related genes in two acute kidney injury models.
Deng W, Wei X, et al.
https://doi.org/10.7717/peerj.10926
"""

from bugs import *
from tcga.utils import download

URLS = {
    # Supp Info 7: DE genes for GSE121190
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7982076/
    'SI7': "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7982076/bin/peerj-09-10926-s007.tsv",

    # Supp Info 8: DE genes for GSE62732
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7982076/
    'SI8': "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7982076/bin/peerj-09-10926-s008.tsv",
}

src_dir = next(p for p in Path.cwd().parents for p in p.glob("**/Deng-Wei-2021/z_*")).with_suffix('')
download = download.to(abs_path=(src_dir / "download_cache"))

if __name__ == '__main__':
    for (k, url) in URLS.items():
        with download(URLS[k]).now.open() as fd:
            df = pd.read_table(fd, index_col=0)
            df.to_csv(src_dir / f"{k}.tsv.gz", compression='gzip', sep='\t')

            print(df)

