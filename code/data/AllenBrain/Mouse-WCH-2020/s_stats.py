# RA, 2021-03-16

"""
A look a the complete dataset.
https://celltypes.brain-map.org/rnaseq/mouse_ctx-hip_10x

This takes forever.
"""

from bugs import *
from tcga.utils import at_most_n
from z_sources import URLS, download

out_dir = mkdir(Path(__file__).with_suffix(''))


def save_df(df: pd.DataFrame, name: str):
    df.to_csv(out_dir / f"{name}.csv", sep='\t')


with download(URLS['expr']).now.open() as fd:
    hist1 = pd.Series(name="nonzeros", data=sum(
        (df != 0).sum(axis=0)
        for df in pd.read_csv(fd, sep=',', chunksize=128, index_col=0)
    ))

    save_df(pd.DataFrame(hist1), "per_gene_nonzeros")

# df.astype(pd.SparseDtype(int, 0)).sparse.density
