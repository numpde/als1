# RA, 2021-05-20

"""
A look at the cell-to-cluster assignment.
"""

import re
from bugs import *

assignment = pd.read_csv(unlist1(Path(__file__).with_suffix('').glob("*.csv")), index_col=0)

assignment = assignment.assign(tag=[unlist1(re.findall(r"[TCGA]{16}", i)) for i in assignment.index])
assignment = assignment.assign(cls=[unlist1(re.findall(r"[12345]$", i)) for i in assignment.index])

cls = {'1': "healthy1", '2': "healthy2", '3': "healthy3", '4': "EAE1", '5': "EAE2"}
assignment = assignment.assign(sample=assignment.apply(lambda s: f"{cls[s.cls]}_{s.tag}", axis=1))

assignment = assignment.set_index('sample', verify_integrity=True).astype(str)

if __name__ == '__main__':
    df_meta = pd.read_table(unlist1(Path(__file__).parent.glob("*/meta.csv.gz")), sep='\t', index_col=0)
    df_meta = df_meta.assign(cluster=assignment.x)
    print(df_meta)
