#!/usr/bin/env bash

DIR="download_cache"

echo "*.feather" > $DIR/.gitignore && git add $DIR/.gitignore

url="https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-500bp-upstream-7species.mc9nr.feather"
wget -nc $url -P $DIR

url="https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-tss-centered-10kb-7species.mc9nr.feather"
wget -nc $url -P $DIR

