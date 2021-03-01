# RA, 2021-04-24

"""
Prepare the scRNA-seq dataset from:

Disease-specific oligodendrocyte lineage cells arise in multiple sclerosis.
Falcão AM, van Bruggen D, Marques S, Meijer M et al.
Nat Med 2018 Dec;24(12):1837-1844.
https://pubmed.ncbi.nlm.nih.gov/30420755/

The GSE page contains "raw" data:
Submission date: May 02, 2018
Last update date: Nov 13, 2019
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113973

The metadata with cell types are from
https://cells.ucsc.edu/?ds=oligo-lineage-ms

Usage:
    from z_sources import df_meta, df_expr

df_meta is Samples x Features
df_expr is Samples x Genes
"""

from tcga.utils import download
import tcga.utils
import pandas
import pathlib
import os
import re

src_dir = next(p for p in pathlib.Path.cwd().parents for p in p.glob("**/CasteloBranco-2018/z_*")).with_suffix('')

download = download.to(abs_path=(src_dir / "download_cache"))

URLS = {
    'GSE98816_miniml': "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113973/miniml/GSE113973_family.xml.tgz",

    # From https://cells.ucsc.edu/?ds=oligo-lineage-ms
    'expr': "https://cells.ucsc.edu/oligo-lineage-ms/exprMatrix.tsv.gz",
    'meta': "https://cells.ucsc.edu/oligo-lineage-ms/meta.tsv",

    # The key is used for cell ID
    'SS2_16_290': "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113973/suppl/GSE113973_SS2_16_290_counts.tab.gz",
    'SS2_16_291': "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113973/suppl/GSE113973_SS2_16_291_counts.tab.gz",
    'SS2_16_292': "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113973/suppl/GSE113973_SS2_16_292_counts.tab.gz",
    'SS2_16_293': "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113973/suppl/GSE113973_SS2_16_293_counts.tab.gz",
    'SS2_16_294': "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113973/suppl/GSE113973_SS2_16_294_counts.tab.gz",
    'SS2_16_295': "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113973/suppl/GSE113973_SS2_16_295_counts.tab.gz",
}

# if __name__ == '__main__':
#     print("Download to:", tcga.utils.relpath(download.local_folder))
#
#     for (_, url) in URLS.items():
#         print("Downloading:", url)
#         json.dumps(download(url).now.meta, indent=2)


# Expression from UCSC cellbrowser (not useful)
def make_df_expr_ucsc():
    import gzip
    with download(URLS['expr']).now.open(mode='rb') as fd:
        with gzip.open(fd, mode='rb') as gz:
            return pandas.read_table(gz, index_col=0)


# Metadata from GSE (not useful)
def make_df_meta_gse() -> pandas.DataFrame:
    import xml.etree.ElementTree as ET
    import tarfile
    import re
    import pandas as pd
    from tcga.utils import unlist1

    with download(URLS['GSE98816_miniml']).now.open(mode='rb') as tf:
        with tarfile.open(fileobj=tf, mode='r') as tar:
            et = ET.parse(source=tar.extractfile(unlist1(tar))).getroot()

            # Namespace a la '{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}'
            ns = unlist1(re.findall(r"({.*}).*", et.tag))

            c1: ET.Element
            df_meta = pd.DataFrame(
                data=(
                    {
                        'gsm': c1.attrib['iid'],
                        'sra': unlist1(c1.findall("./*/[@type='SRA']")).attrib["target"].strip(),
                        'taxid': unlist1(c1.findall("*/*/[@taxid]")).attrib["taxid"].strip(),
                        'source': unlist1(c1.findall(f"*/{ns}Source")).text.strip(),
                        'biosample': unlist1(c1.findall("./*/[@type='BioSample']")).attrib["target"].strip(),
                        'strain': unlist1(c1.findall("*/*/[@tag='strain']")).text.strip(),
                        'tissue': unlist1(c1.findall("*/*/[@tag='tissue']")).text.strip(),
                        'celltype': unlist1(c1.findall("*/*/[@tag='celltype']")).text.strip(),
                        'title': unlist1(c1.findall(ns + "Title")).text.strip(),
                        'accession': unlist1(c1.findall(ns + "Accession")).text.strip(),
                        'description': sorted(unlist1(c1.findall(ns + "Description")).text.strip().split('\n')),
                    }
                    for c1 in et.findall(ns + "Sample")
                )
            )

            # Remove common prefix from the description column
            df_meta = df_meta.assign(
                sample=df_meta.description.str.slice(len(os.path.commonprefix(list(df_meta.description))))
            )

            df_meta = df_meta.drop(columns='description')

        return df_meta


# Metadata from UCSC cellbrowser (useful)
def make_df_meta_ucsc():
    with download(URLS['meta']).now.open(mode='r') as fd:
        df = pandas.read_table(fd, index_col=0)
        df = df.assign(celltype=df.Renamed_clusternames)
        p = re.compile(r"_([0-9]+)[.]tab[.]([A-Z][0-9]+)")
        df.index = ["SS2_16_{}_{}".format(*tcga.utils.unlist1(p.findall(i))) for i in df.index]
        return df


# Expression data from GSE (useful)
def make_df_expr_gse():
    def url2df(k):
        with download(URLS[k]).now.open(mode='rb') as fd:
            df = pandas.read_table(fd, compression='gzip', index_col=0).astype(int).sort_index(axis=1)
            df.columns = [f"{k}_{c}" for c in df.columns]
            return df

    df = pandas.concat(axis=1, verify_integrity=True, objs=(
        url2df(k)
        for k in URLS
        if k.startswith("SS2")
    ))

    # Collapse duplicate genes
    df = df.groupby(df.index).median().round().astype(int)

    # Transpose to Samples x Genes
    df = df.T

    return df


if __name__ == '__main__':
    df_meta = make_df_meta_ucsc().sort_index()
    df_expr = make_df_expr_gse().reindex(df_meta.index)

    assert df_meta.index.equals(df_expr.index)

    df_meta.to_csv(src_dir / "meta.csv.gz", compression='gzip', sep='\t')
    df_expr.to_csv(src_dir / "expr.csv.gz", compression='gzip', sep='\t')

df_meta = pandas.read_table(src_dir / "meta.csv.gz", sep='\t', index_col=0)
df_expr = pandas.read_table(src_dir / "expr.csv.gz", sep='\t', index_col=0)

assert df_meta.index.equals(df_expr.index)
