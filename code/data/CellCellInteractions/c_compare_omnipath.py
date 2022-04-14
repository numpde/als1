# RA, 2021-07-06

import pandas as pd
from pathlib import Path
from z_sources import df as df_baderlab

df_omnipath = pd.read_table(next(p for p in Path(__file__).parents for p in p.glob("**/omnipath_lr.csv.gz")))

print(df_omnipath)
print(df_baderlab)
