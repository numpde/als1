The script [a_download.py](a_download.py)
downloads 
the "Gene expression matrix"
([csv](https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/matrix.csv), 7GB)
and
the "Table of cell metadata"
([csv](https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/metadata.csv), 23MB)
from
the [Allen Brain M1](https://portal.brain-map.org/atlases-and-data/rnaseq/human-m1-10x)
dataset
directly to zipped files
of about 560MB and 2MB. 

The script [b_reduced.py](b_reduced.py)
subsets the gene expression matrix
to [140 genes](b_reduced/marker_genes.csv)
that are mentioned 
in the metadata
and writes this as (genes x cells) table.
It also processes the metadata
into a slimmer table.
The output is in
the folder [b_reduced](b_reduced).

The script [c_visualize.py](c_visualize.py)
computes
(an instance of)
t-SNE coordinates
for the reduced expression dataset,
for all cells and for each cell type individually.
The results are in the folder
[c_visualize](c_visualize).
