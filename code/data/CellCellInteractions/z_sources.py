# RA, 2021-07-06

"""
Cell-Cell Interaction Database by Bader lab
Version 1.0 - Built April 25, 2017 and contains iRefIndex version 14, Pathway Commons version 8 and BioGRID version 3.4.147
http://baderlab.org/CellCellInteractions#Download_Data

Numbers on the website:
    Receptors - 1851 genes
    Ligands - 1593 genes
    ECM - 433 genes

Fields:
    AliasA - main Alias for molecule A (often the recognized gene symbol)
    AliasB- main Alias for molecule B (often the recognized gene symbol)
    uidA - unique identifier for molecule A (depending on the source database this can be one of the following types uniprot, refseq, entrez gene id, ensembl)
    uidB - unique identifier for molecule A (depending on the source database this can be one of the following types uniprot, refseq, entrez gene id, ensembl)
    altA - list of alternate identifiers for molecule A.
    altB - list of alternate identifiers for molecule B.
    aliasA - list of alternate aliases for molecule A.
    aliasB - list of alternate aliases for molecule B.
    method - list of psi-mi terms indicating experimental methods used to discover interaction.
    author - text listing authors
    pmids - list of pmids associated with the interaction.
    taxa - taxon id for molecule A.
    taxb - taxon id for molecule B.
    interactionType - list of psi-mi terms indicating the type of interactions it is.
    sourcedb - source database.
    interactionIdentifier - source database interaction identifier
    confidence - confidence of interaction as supplied by database source
"""

import gzip
import pandas
from pathlib import Path
from tcga.utils.download import download

url = "http://baderlab.org/CellCellInteractions?action=AttachFile&do=get&target=receptor_ligand_interactions_mitab_v1.0_April2017.txt.zip"

download = download.to(abs_path=Path(__file__).with_suffix('').resolve())

with download(url).now.open(mode='rb') as fd:
    with gzip.open(fd, mode='r') as zf:
        df = pandas.read_table(zf, low_memory=False)
