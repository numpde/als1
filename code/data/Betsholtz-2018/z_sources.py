# RA, 2021-03-24

"""
Prepare the scRNA-seq dataset from:

Single cell RNA-seq of mouse brain vascular transcriptomes
Vanlandewijck M, He L, Mäe M, Andrae J, Betsholtz C

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98816

Usage:
    from z_sources import df_meta, df_expr, df_mrkr


"""

from tcga.utils import download
import tcga.utils
import pandas
import pathlib
import os
import json
import openpyxl

src_dir = tcga.utils.mkdir(pathlib.Path(__file__).with_suffix(''))
download = download.to(abs_path=(src_dir / "download_cache"))

URLS = {
    'GSE98816': "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE98nnn/GSE98816/suppl/GSE98816_Brain_samples_raw_read_counts_matrix.txt.gz",
    'GSE98816_miniml': "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE98nnn/GSE98816/miniml/GSE98816_family.xml.tgz",

    # https://figshare.com/collections/_/4077260
    # The file contains the description for 3436 single cells from mouse brain and 1504 single cells from mouse lung
    'sc_description': "https://ndownloader.figshare.com/files/11188505",
}

# https://www.nature.com/articles/nature25739/figures/1
# 1c
markers = {
    'PC': "Pdgfrb Cspg4 Anpep Rgs5 Cd248 Abcc9 Vtn S1pr3",
    'SMC': "Acta2 Tagln Myh11 Myl9 Mylk Sncg Cnn1 Pln",
    'MG': "Csf1r Emr1 Cd68 Cd53 Cd48 Cd84 C1qa Fcgr1",
    'FB': "Pdgfra Lum Dcn Col3a1 Col5a1 Col8a2 Col12a1 Mmp2",
    'OL': "Mobp Plp1 Mog Cldn11 Mag Gjc2 Mal Cnp",
    'EC': "Pecam1 Kdr Flt1 Tie1 Tek Icam2 Podxl Ptprb",
    'AC': "Aldh1l1 Fgfr3 Slc4a4 Slc6a11 Slc7a10 Mlc1 Slc1a3 Cldn10",
}


def make_df_desc() -> pandas.DataFrame:
    import warnings
    with download(URLS['sc_description']).now.open(mode='rb') as fd:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            values = openpyxl.load_workbook(fd).active.values
        # Note: order of arguments matters
        df = pandas.DataFrame(columns=next(values), data=list(values))
        df = df.rename(columns={'GSM ID': "gsm", 'annoated cell types': "celltype"})
        return df


def make_df_meta() -> pandas.DataFrame:
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
            # c1 = first(et.findall(ns + "Sample"))
            df_meta = pd.DataFrame(
                data=(
                    {
                        'gsm': c1.attrib['iid'],
                        'sra': unlist1(c1.findall("./*/[@type='SRA']")).attrib["target"].strip(),
                        'taxid': unlist1(c1.findall("*/*/[@taxid]")).attrib["taxid"].strip(),
                        'biosample': unlist1(c1.findall("./*/[@type='BioSample']")).attrib["target"].strip(),
                        'strain': unlist1(c1.findall("*/*/[@tag='strain']")).text.strip().lower(),
                        'tissue': unlist1(c1.findall("*/*/[@tag='tissue']")).text.strip().lower(),
                        'genotype': unlist1(c1.findall("*/*/[@tag='genotype']")).text.strip().lower(),
                        'age': unlist1(c1.findall("*/*/[@tag='age']")).text.strip().lower(),
                        'title': unlist1(c1.findall(ns + "Title")).text.strip(),
                        'accession': unlist1(c1.findall(ns + "Accession")).text.strip(),
                        'description': unlist1(c1.findall(ns + "Description")).text.strip(),
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


def make_df_expr() -> pandas.DataFrame:
    with download(URLS['GSE98816']).now.open(mode='rb') as fd:
        # samples x genes
        df_expr = pandas.read_table(fd, compression='gzip', quotechar='"', index_col=0).T

    # Sort by sample ID
    df_expr = df_expr.sort_index()
    df_expr.index.name = "sample"

    assert df_expr.index.is_unique

    # Also remove common prefix to match df_meta
    df_expr.index = df_expr.index.str.slice(len(os.path.commonprefix(list(df_expr.index))))
    assert df_expr.index.is_unique
    return df_expr


def make_df_markers() -> pandas.DataFrame:
    df_markers = pandas.DataFrame(
        columns=['celltype', 'gene'],
        data=(
            (celltype, gene)
            for (celltype, genes) in markers.items()
            for gene in genes.split()
        )
    )

    df_markers = df_markers.assign(v=1).pivot_table(
        index='gene', columns='celltype', values='v', fill_value=0,
    )

    df_markers = df_markers.astype(int)

    return df_markers


if __name__ == '__main__':
    from tcga.utils import mkdir

    for (_, url) in URLS.items():
        json.dumps(download(url).now.meta, indent=2)

    df_meta = make_df_meta()
    df_meta = df_meta.merge(make_df_desc(), how="inner", on="gsm", suffixes=("", " (desc)"))
    df_meta = df_meta.set_index('sample', verify_integrity=True).sort_index()

    df_expr = make_df_expr()
    df_mrkr = make_df_markers()

    assert df_meta.index.equals(df_expr.index)

    df_meta.to_csv(src_dir / "meta.csv.gz", compression='gzip', sep='\t')
    df_expr.to_csv(src_dir / "expr.csv.gz", compression='gzip', sep='\t')
    df_mrkr.to_csv(src_dir / "mrkr.csv", sep='\t')

#

df_meta = pandas.read_table(src_dir / "meta.csv.gz", sep='\t', index_col=0)
df_expr = pandas.read_table(src_dir / "expr.csv.gz", sep='\t', index_col=0)
df_mrkr = pandas.read_table(src_dir / "mrkr.csv", sep='\t', index_col=0)

assert df_meta.index.equals(df_expr.index)
