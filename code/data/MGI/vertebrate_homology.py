# RA, 2021-05-30


"""
Prepare some tables from
http://www.informatics.jax.org/downloads/reports/index.html
"""

from bugs import *
from tcga.utils import download


URLS = {
    'Vertebrate homology -- 1. Homology classes': "http://www.informatics.jax.org/downloads/reports/HOM_AllOrganism.rpt",
}

download = download.to(abs_path=(Path(__file__).parent / "download_cache"))
out_dir = mkdir(Path(__file__).parent / "data")

for (k, url) in URLS.items():
    with download(url).now.open(mode='r') as fd:
        df = pd.read_csv(fd, sep='\t', dtype='str')

        # df = df.pivot_table(index=hg, columns='NCBI Taxon ID', values='Symbol', aggfunc=','.join)

        hg = 'HomoloGene ID'
        eg = 'EntrezGene ID'
        tx = 'NCBI Taxon ID'
        sb = 'Symbol'

        df = df[[tx, eg, hg, sb]]

        df = df.drop_duplicates()

        df = df.set_index(df[tx] + "-" + df[eg], verify_integrity=True)
        df.columns = [c.replace(' ', '_') for c in df.columns]

        df.to_csv(out_dir / (Path(__file__).with_suffix('.tsv.gz').name), sep='\t', compression='gzip')
