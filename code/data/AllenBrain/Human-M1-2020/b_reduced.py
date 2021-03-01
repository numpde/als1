# RA, 2020-12-14

"""
From the AllenBrain M1 data,

compile a metadata file `meta.csv.gz` with
 - cell types
 - library sizes
 - etc.

and

subset the expression file to marker genes to `data.csv.gz`.
"""

from pathlib import Path
from datetime import datetime, timezone
from contextlib import redirect_stdout

import pandas as pd

from a_download import markers, URLS, download
from tcga.utils import mkdir

out_dir = mkdir(Path(__file__).with_suffix(''))

def main():

    # Save marker genes to file
    pd.Series(sorted(markers), name="markers").to_csv(out_dir / "marker_genes.csv", sep='\t', index=False)

    # Subset sample expression to marker genes
    with download(URLS['expr']).now.open() as fd:
        df_expr = pd.read_csv(fd, sep=',', index_col=0, usecols=['sample_name', *sorted(markers)]).T
        df_expr.index.name = "gene_name"
        df_expr.to_csv(out_dir / "data.csv.gz", sep='\t', compression='gzip')

        with redirect_stdout((out_dir / "data_readme.txt").open(mode='w')):
            print("Source:    ", URLS['expr'])
            print("Local copy:", download(URLS['expr']).now.local_file.name)
            print("Datetime:  ", datetime.now(tz=timezone.utc).strftime("%Z-%Y%m%d-%H%M%S"))
            print("Script:    ", Path(__file__).name)

    # Compute library size for each sample
    with download(URLS['expr']).now.open() as fd:
        libsize = pd.Series(
            data=pd.concat(
                df.sum(axis=1)
                for df in pd.read_csv(fd, sep=',', index_col=0, chunksize=1024)
                # for df in [pd.read_csv(fd, sep=',', index_col=0, nrows=2)]  # DEBUG
            ),
            name="libsize",
            dtype=int,
        )

    # Compile a reduced meta table
    with download(URLS['meta']).now.open() as fd:
        df_meta = pd.read_csv(fd, sep=',', index_col=0)
        assert (df_meta.shape == (len(df_meta), 38))

        df_meta.cluster_label.str.split(" ", n=2, expand=True).rename(
            columns={0: 'celltype', 1: 'layer', 2: 'markers'}
        ).assign(
            gender=df_meta.donor_sex_label,
            region=df_meta.region_label,
            donor=df_meta.external_donor_name_label,
            libsize=libsize,
            # Add other useful fields here
        ).astype(
            {
                'libsize': 'Int64',
            }
        ).to_csv(
            out_dir / "meta.csv.gz", sep='\t', compression='gzip',
        )

        with redirect_stdout((out_dir / "meta_readme.txt").open(mode='w')):
            print("Source:    ", URLS['meta'])
            print("Local copy:", download(URLS['meta']).now.local_file.name)
            print("Datetime:  ", datetime.now(tz=timezone.utc).strftime("%Z-%Y%m%d-%H%M%S"))
            print("Script:    ", Path(__file__).name)


if __name__ == '__main__':
    main()

