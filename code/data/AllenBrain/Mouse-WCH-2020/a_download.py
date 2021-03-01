# RA, 2021-03-06

"""
Download to the expression table and metadata into local zipfiles.

See z_sources.py.
"""

import urllib.request
from contextlib import closing

from bugs import *
from twig import log
from z_sources import download, URLS

out_dir = Path(__file__).with_suffix('')


def download_meta():
    log.info("Downloading the meta data.")
    with download(URLS['meta']).now.open() as rd:
        df_meta = pd.read_csv(rd, sep=',', index_col=0)

    # print(json.dumps(Counter(df_meta.subclass_label), indent=2))
    summary = {
        "NaN": 4014,
        "L2/3 IT ENTl": 5764, "L2 IT RHP": 7599, "L2/3 IT PPP": 34084, "L2 IT ENTl": 4068, "L4/5 IT CTX": 253722,
        "L5 PT CTX": 16783, "L5 IT TPE-ENT": 5525, "L2/3 IT CTX-1": 117565, "L3 IT ENT": 13789,
        "L3 RSP-ACA": 4214, "L2/3 IT CTX-2": 7141, "L6 IT CTX": 79403, "L5 PPP": 1240, "L6 IT ENTl": 1169,
        "L5 IT CTX": 44889, "L5 NP CTX": 29378, "L6 CT CTX": 135241, "L6b CTX": 13114, "L6b/CT ENT": 20789,
        "NP SUB": 1949, "NP PPP": 2695,
        "V3d": 66,
        "Meis2": 1,
        "Lamp5": 38464,
        "Vip": 41626,
        "Sncg": 11573,
        "Sst": 42310,
        "Pvalb": 31088,
        "Sst Chodl": 1906,
        "DG": 58754,
        "CA1-ProS": 16141,
        "Car3": 21538,
        "SUB-ProS": 4406,
        "CT SUB": 6012,
        "CA2": 336,
        "CA3": 1899,
        "CR": 268,
        "Oligo": 7685,
        "Astro": 3119,
        "SMC-Peri": 198,
        "Endo": 746,
        "VLMC": 129,
        "Micro-PVM": 636
    }

    log.info("Making the dummy set.")
    with closing(urllib.request.urlopen(url=URLS['expr'])) as rd:
        df_data = pd.read_csv(rd, sep=',', index_col=0, nrows=50).iloc[:, 0:101]
    df_meta = df_meta.loc[df_data.index]
    df_meta.to_csv(download.local_folder / "dummy_meta.csv", sep=',')
    df_data.to_csv(download.local_folder / "dummy_data.csv", sep=',')
    log.info("Dummy set done.")

    log.info("Run b_reduced.py to download the reduced dataset.")


def download_expr():
    log.info("Downloading the expr data.")
    log.info(download(URLS['expr']).now.meta)


def main():
    download_expr()
    download_meta()


if __name__ == '__main__':
    main()
