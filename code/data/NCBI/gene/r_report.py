# RA, 2021-05-26

from pathlib import Path
import contextlib
import pandas as pd
from tcga.utils import mkdir

for src in Path(__file__).parent.glob("*table/*.tsv.gz"):
    df = pd.read_table(src)

    with (mkdir(Path(__file__).with_suffix('')) / f"{Path(src.stem)}.txt").open(mode='w') as fd:
        with contextlib.redirect_stdout(fd):
            print("Number of records:", len(df))
            print("gene_is is unique:", df.gene_id.is_unique)

            print("Preview:")
            print(df[~df.summary.isna()].head(10).to_markdown(index=False))

