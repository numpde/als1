# RA, 2020-12-20

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.manifold import TSNE

from tcga.utils import mkdir, assert_exists, relpath
from plox import Plox
from twig import log

seed = 2


def tsne(df, out_csv: Path):
    log.info(relpath(out_csv))

    if out_csv.is_file():
        X = pd.read_csv(out_csv, sep='\t', compression="infer", index_col=0)
        assert (X.index.name == "sample_name")
        assert (list(X.columns) == ['x', 'y'])
    else:
        assert (df.index.name == "gene_name")
        X = pd.DataFrame(
            index=pd.Series(df.columns, name="sample_name"),
            columns=['x', 'y'],
            data=TSNE(random_state=seed).fit_transform(df.T),
        )
        X.to_csv(out_csv, sep='\t', compression="gzip")

    # https://matplotlib.org/tutorials/introductory/customizing.html
    style = {
        'legend.fontsize': "xx-small",
        'legend.framealpha': 0.5,
    }

    with Plox(style) as px:
        px.a.plot(X.x, X.y, '.', ms=(10 / np.log10(len(X))))
        # px.a.legend()
        px.a.axis('off')
        px.f.savefig(out_csv.with_suffix(".png"))


def main():
    out_dir = mkdir(Path(__file__).with_suffix(''))

    datapath = Path(__file__).parent / "b_reduced"
    read_csv = (lambda f: pd.read_csv(assert_exists(f), sep='\t', index_col=0))

    data = read_csv(datapath / "data.csv.gz")
    meta = read_csv(datapath / "meta.csv.gz")

    # Normalize
    data = data / meta.libsize

    # By cell type
    for (t, df) in data.groupby(meta.celltype[data.columns], axis=1):
        tsne(df, out_dir / F"tsne_{t}.csv.gz")

    # All
    tsne(data, out_dir / F"tsne.csv.gz")


if __name__ == '__main__':
    main()
