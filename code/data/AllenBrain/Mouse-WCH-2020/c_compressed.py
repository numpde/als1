# RA, 2021-03-14

"""
Download the big expression table and
compress it on the fly.
"""

raise NotImplemented("This still blows the RAM.")

import urllib.request
from contextlib import closing

from progressbar import progressbar

from twig import log
from bugs import *

from a_download import URLS, download

from scipy.sparse import csr_matrix as matrix, vstack


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

    def peek(x, text=None):
        if text is None:
            log.info(x)
        else:
            log.info(text)
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

    # Metadata
    with meta_open_remote() as fd:
        if meta_file.exists():
            log.warning(f"File will be overwritten when done: {relpath(meta_file)}")

        df_meta = pd.read_csv(fd, sep=PARAM.remote_sep, index_col=0)
        assert (df_meta.shape == (len(df_meta), 56))

        nsamples_total = len(df_meta)
        log.info(f"Based on metadata, there are {nsamples_total} samples.")

        df_meta.to_csv(meta_file, sep=PARAM.local_sep, compression='gzip')

        log.info(f"Size of reduced dataset: {len(df_meta)}.")
        log.info(f"Finished {relpath(meta_file)}")

        del df_meta

    # Collect expression
    with data_open_remote() as rd:
        if data_file.exists():
            log.warning(f"File will be overwritten when done: {relpath(data_file)}")

        if PARAM.DUMMY_MODE:
            chunksize = 24
        else:
            chunksize = 128

        nchunks_expected = (nsamples_total // chunksize) + bool((nsamples_total % chunksize))

        log.info(f"Chunksize is {chunksize} rows. Expect {nchunks_expected} chunks.")
        log.info(f"Downloading.")

        df_data = pd.concat(axis=0, objs=(
            chunk.astype(pd.SparseDtype('int', fill_value=0))
            for chunk in progressbar(
            pd.read_csv(
                rd, sep=PARAM.remote_sep, index_col=0, chunksize=chunksize
            ),
            max_value=nchunks_expected
        )
        ))

        log.info(f"Sparse density: {df_data.sparse.density}")

        # genes x samples
        df_data = df_data.T

        df_data.to_csv(data_file, sep=PARAM.local_sep, compression='gzip')

        log.info(f"Data has {len(df_data.columns)} samples.")
        log.info(f"Finished {relpath(data_file)}")


if __name__ == '__main__':
    main()
