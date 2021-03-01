# RA, 2021-03-07

"""
Extract cells of certain types, e.g. VLMC.
"""

from pathlib import Path

import pandas as pd

from twig import log
from tcga.utils import at_most_n, mkdir, relpath

from a_download import download, URLS


class PARAM:
    subclass_of_interest = ["VLMC", "Astro", "Endo", "OPC", "Micro-PVM"]  # , "Oligo"
    remote_sep = ','
    local_sep = '\t'
    DUMMY_MODE = False


def main():
    from a_download import df_meta

    out_dir = mkdir(Path(__file__).with_suffix(''))
    data_file = out_dir / "data.csv.gz"
    meta_file = out_dir / "meta.csv.gz"

    if True:
        df_meta = df_meta[df_meta.subclass_label.isin(PARAM.subclass_of_interest)]
        log.info(f"New subset of cells: {dict(df_meta.subclass_label.value_counts())}")

        df_meta.to_csv(meta_file, sep=PARAM.local_sep, compression='gzip')

        log.info(f"Size of reduced dataset: {len(df_meta)}.")
        log.info(f"Finished {relpath(meta_file)}")

    with download(URLS['expr']).now.open(mode='r') as fd:
        log.info(f"Reducing the expression data.")

        df_data = pd.concat(axis=0, objs=(
            df[df.index.isin(df_meta.index)]
            for df in pd.read_csv(fd, sep=PARAM.remote_sep, index_col=0, chunksize=1024)
            if any(df.index.isin(df_meta.index))
        ))

        # genes x samples
        df_data = df_data.T

        df_data.to_csv(data_file, sep=PARAM.local_sep, compression='gzip')

        log.info(f"Data has {len(df_data.columns)} samples, expected {len(df_meta)}.")
        log.info(f"Finished {relpath(data_file)}")


if __name__ == '__main__':
    main()
