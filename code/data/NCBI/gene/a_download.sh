#!/usr/bin/env bash

# This requires the ncbi tool `gene2xml`
# sudo apt install ncbi-tools-bin

DIR="download_cache"
mkdir -p "$DIR"

echo "*.xml.gz" > $DIR/.gitignore && git add $DIR/.gitignore

URL="ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/ASN_BINARY/Mammalia/Mus_musculus.ags.gz"
wget -S -O - "$URL" | gzip -d | gene2xml -b T -l T | gzip -9 > "$DIR/mm.xml.gz"

URL="ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/ASN_BINARY/Mammalia/Homo_sapiens.ags.gz"
wget -S -O - "$URL" | gzip -d | gene2xml -b T -l T | gzip -9 > "$DIR/hs.xml.gz"
