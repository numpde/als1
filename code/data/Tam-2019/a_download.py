# RA, 2020-11-04

"""
References:

[1]
Postmortem cortex samples identify molecular subtypes of ALS:
retrotransposon activation, oxidative stress, and activated glia.
Tam et al. incl. The NYGC ALS Consortium.
Cell Rep, 2019.
https://www.sciencedirect.com/science/article/pii/S221112471931263X

[GSE124439]
RNA-seq analysis of ALS patient samples and
individuals without neurological disorders
Contributor: NYGC ALS Consortium
https://www.nygenome.org/als-consortium/
"""

import os
import re

from tcga.utils import download, unlist1
from pathlib import Path
from contextlib import redirect_stderr, redirect_stdout

import pandas as pd
import tarfile
import xml.etree.ElementTree as ET

datapath = Path(__file__).with_suffix('')
download = download.to(abs_path=(datapath / "cache"))

URLS = {
    # GSE124439_RAW.tar
    'GSE124439': "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE124439&format=file",

    # 'GSE124439_series': "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124439/matrix/GSE124439_series_matrix.txt.gz",

    # 'GSE124439_soft': "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124439/soft/GSE124439_family.soft.gz",
    'GSE124439_miniml': "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124439/miniml/GSE124439_family.xml.tgz",

    # For downloading Fasta files from FTP
    'manifest': "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA512012&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,fastq_ftp,submitted_ftp,sra_ftp&format=tsv&download=true&limit=0",
}

for url in URLS.values():
    data = download(url).now
    print(data.meta)

with download(URLS['manifest']).now.open(mode='r') as fd:
    df = pd.read_table(fd)
    df.to_csv(datapath / "GSE124439_manifest.tsv", compression=None, sep='\t', index=False)

# # Extract HRA-28394, HRA-89876-02 but not HRA-23938-b38
# extract_hra = (lambda t: unlist1(re.findall(r"(HRA-[0-9]+(?:-[0-9]+)?)", t)))
# Extract HRA-28394 only, drop suffix
extract_hra = (lambda t: unlist1(re.findall(r"(HRA-[0-9]+)", t)))

with download(URLS["GSE124439_miniml"]).now.open(mode='rb') as tf:
    with tarfile.open(fileobj=tf, mode='r') as tar:
        et = ET.parse(source=tar.extractfile(unlist1(tar))).getroot()
        namespace = "{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}"
        c1: ET.Element
        df = pd.DataFrame(
            data=(
                {
                    'gsm': c1.attrib['iid'],
                    'hra': extract_hra(unlist1(c1.findall(namespace + "Title")).text.strip()),
                    'sra': unlist1(c1.findall("./*/[@type='SRA']")).attrib["target"].strip(),
                    'biosample': unlist1(c1.findall("./*/[@type='BioSample']")).attrib["target"].strip(),
                    'tissue': unlist1(c1.findall("*/*/[@tag='tissue type']")).text.strip().lower(),
                    'cns_subregion': unlist1(c1.findall("*/*/[@tag='cns subregion']")).text.strip().lower(),
                    'subject_id': unlist1(c1.findall("*/*/[@tag='subject id']")).text.strip(),
                    'sample_group': unlist1(c1.findall("*/*/[@tag='sample group']")).text.strip(),
                }
                for c1 in et.findall(namespace + "Sample")
            )
        )

        df.sample_group = df.sample_group.map(
            {
                'Non-Neurological Control': "non_neuro",
                'Other Neurological Disorders': "other_neuro",
                'ALS Spectrum MND': "ALS"
            }
        )

        df = df.set_index("hra", verify_integrity=True)
        df.to_csv(datapath / "GSE124439_meta.txt.gz", compression="gzip", sep='\t')

        df = df.sort_index(axis=0)

        hra = df.index

        # print(df.sample_group.value_counts())

with download(URLS['GSE124439']).now.open(mode='rb') as tf:
    with tarfile.TarFile(fileobj=tf, mode='r') as tar:
        df = pd.concat(
            axis=1,
            ignore_index=False,
            verify_integrity=True,
            join="outer",
            objs=(
                pd.read_csv(tar.extractfile(f), compression="gzip", sep='\t', quotechar='"', index_col=0).astype(int)
                for f in tar
            )
        ).sort_index()

        df.columns = pd.Index(name="hra", data=list(map(extract_hra, df.columns)))
        df = df.sort_index(axis=1)

        assert df.columns.is_unique
        assert df.columns.equals(hra)

        df.to_csv(datapath / "GSE124439_data.txt.gz", compression="gzip", sep='\t')

# x = pd.DataFrame(data={'df': df.columns, 'hra': hra})
# x['='] = (x.df == x.hra)
# print(x.to_markdown())
