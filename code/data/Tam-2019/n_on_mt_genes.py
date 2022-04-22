from pathlib import Path

import pandas as pd

from tcga.utils import unlist1

df = pd.read_table(unlist1(Path(__file__).parent.glob("**/*data.txt.gz")))

# PART I

# Some genes from the MitoCarta3.0
# https://www.broadinstitute.org/mitocarta/mitocarta30-inventory-mammalian-mitochondrial-proteins-and-pathways

# https://www.genecards.org/cgi-bin/carddisp.pl?gene=CYC1
CYC1 = ["CYC1", "MC3DN6", "UQCR4", "P08574", "gene:1537"]

# https://www.genecards.org/cgi-bin/carddisp.pl?gene=SDHB
SDHB = ["SDHB", "CWS2", "IP", "PGL4", "SDH", "SDH1", "SDH2", "SDHIP", "P21912", "gene:6390"]

# https://www.genecards.org/cgi-bin/carddisp.pl?gene=COQ7
COQ7 = ["COQ7", "CAT5", "CLK-1", "CLK1", "COQ10D8", "Q99807", "gene:10229"]

assert any(gene in list(df['gene/TE']) for gene in CYC1)
assert any(gene in list(df['gene/TE']) for gene in SDHB)
assert any(gene in list(df['gene/TE']) for gene in COQ7)

# PART II

# However, MT- genes seem absent

assert not any(df['gene/TE'].str.startswith("MT-"))
