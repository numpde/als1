# RA, 2021-05-26

"""
On-the-fly conversion of *.xml.gz to a slim tab-separated table.

With help from
https://github.com/hsiaoyi0504/gene_dictionary/blob/4d84c5c4e4acb6103c6d22f8cc0e101e9535ab85/extract.py
"""

import xml.etree.ElementTree as ET
import gzip

from pathlib import Path

from tcga.utils import at_most_n, mkdir, unlist1, relpath


def parse(src):
    def find_text(e, path):
        e = e.findall(path)
        if e:
            return unlist1(e).text
        else:
            return ""

    yield ('gene_id', 'description', 'summary')

    with gzip.open(src, mode='r') as fd:
        e: ET.Element
        for (x, e) in ET.iterparse(fd, events=['end']):
            if e.tag == "Entrezgene":
                gene_id = find_text(e, "**/Gene-track_geneid")
                descptn = find_text(e, "**/Gene-ref_desc")
                summary = find_text(e, "Entrezgene_summary")

                yield (gene_id, descptn, summary)

                e.clear()


def main():
    for src in Path(__file__).parent.glob("download_cache/*.xml.gz"):
        trg = mkdir(Path(__file__).with_suffix('')) / Path(src.stem).with_suffix('.tsv.gz')
        if trg.is_file():
            print("Skipping existing file:", relpath(trg))
        else:
            print("Writing:", relpath(trg))
            with gzip.open(trg, mode='wt') as fd:
                for t in parse(src):
                    print(*t, sep='\t', file=fd)


if __name__ == '__main__':
    main()
