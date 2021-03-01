# RA, 2021-03-07

"""
Download, or open form local zip, the big expression table
and filter it for certain cell types (see `subclass_of_interest`).
"""

import urllib.request
from contextlib import closing

from progressbar import progressbar

from bugs import *
from twig import log

from z_sources import URLS, download


class PARAM:
    subclass_of_interest = ["VLMC", "Astro", "Endo"]
    remote_sep = ','
    local_sep = '\t'

    # In dummy mode, the local dummy_data and dummy_meta will be sourced
    DUMMY_MODE = False


def main():
    if PARAM.DUMMY_MODE:
        log.info("WE'RE IN DUMMY MODE.")

    out_dir = mkdir(Path(__file__).with_suffix(''))

    if PARAM.DUMMY_MODE:
        data_file = out_dir / "dummy_data.csv.gz"
        meta_file = out_dir / "dummy_meta.csv.gz"
    else:
        data_file = out_dir / "data.csv.gz"
        meta_file = out_dir / "meta.csv.gz"

    def peek(x):
        log.info(x)
        return x

    def meta_open_remote():
        if PARAM.DUMMY_MODE:
            return open(unlist1(download.local_folder.glob("dummy_meta.csv")), mode='r')
        else:
            return download(URLS['meta']).now.open()

    def data_open_remote():
        if PARAM.DUMMY_MODE:
            return open(unlist1(download.local_folder.glob("dummy_data.csv")), mode='r')
        else:
            return closing(urllib.request.urlopen(url=URLS['expr']))

    # Make a reduced metadata file
    with meta_open_remote() as fd:
        if meta_file.exists():
            log.warning(f"File will be overwritten when done: {relpath(meta_file)}")

        df_meta = pd.read_csv(fd, sep=PARAM.remote_sep, index_col=0)
        assert (df_meta.shape == (len(df_meta), 56))

        nsamples_total = len(df_meta)
        log.info(f"Based on the metadata, there are {nsamples_total} in total.")

        # Subset df_meta to samples of interest

        if PARAM.DUMMY_MODE:
            ix = df_meta.sample(12, random_state=5, replace=False).index
        else:
            ix = df_meta.index[df_meta.subclass_label.isin(PARAM.subclass_of_interest)]

        df_meta = df_meta[df_meta.index.isin(ix)]
        df_meta.to_csv(meta_file, sep=PARAM.local_sep, compression='gzip')

        log.info(f"Size of reduced dataset: {len(df_meta)}.")
        log.info(f"Finished {relpath(meta_file)}")

    # Make a reduced expression data file
    with data_open_remote() as rd:
        if data_file.exists():
            log.warning(f"File will be overwritten when done: {relpath(data_file)}")

        if PARAM.DUMMY_MODE:
            chunksize = 24
        else:
            chunksize = 1024

        nchunks_expected = (nsamples_total // chunksize) + bool((nsamples_total % chunksize))

        log.info(f"Chunksize is {chunksize} rows. Expect {nchunks_expected} chunks.")
        log.info(f"Downloading.")

        df_data = pd.concat(axis=0, objs=[
            chunk[chunk.index.isin(df_meta.index)]
            for chunk in progressbar(
                pd.read_csv(rd, sep=PARAM.remote_sep, index_col=0, chunksize=chunksize),
                max_value=nchunks_expected
            )
            if any(chunk.index.isin(df_meta.index))
        ])

        # genes x samples
        df_data = df_data.T

        df_data.to_csv(data_file, sep=PARAM.local_sep, compression='gzip')

        log.info(f"Data has {len(df_data.columns)} samples, expected {len(df_meta)}.")
        log.info(f"Finished {relpath(data_file)}")


if __name__ == '__main__':
    main()
