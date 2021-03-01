# RA, 2021-03-24

from bugs import *
from z_sources import df_expr, df_meta

sample = pd.Series(df_expr.index)

with (mkdir(Path(__file__).with_suffix('')) / f"cells.txt").open(mode='w') as fd:
    sample.str.split('_').apply(first).value_counts().to_csv(fd, sep='\t', header=False)

# print("Total:", len(sample))
