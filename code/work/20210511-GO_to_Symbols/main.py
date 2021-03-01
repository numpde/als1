# RA, 2021-05-11

from bugs import *
from networkx.algorithms.traversal.depth_first_search import dfs_tree

import obonet
import networkx as nx

URLS = {
    'expr': "https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hip_10x/matrix.csv",
    'meta': "https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hip_10x/metadata.csv",
    'obo': "http://purl.obolibrary.org/obo/go.obo",
}

roots = {
    # http://www.informatics.jax.org/vocab/gene_ontology/GO:1903715 (regulation of aerobic respiration)
    'GO:1903715',

    # http://www.informatics.jax.org/vocab/gene_ontology/GO:0009060 (aerobic respiration)
    'GO:0009060',

    # http://www.informatics.jax.org/vocab/gene_ontology/GO:0010762 (regulation of fibroblast migration)
    'GO:0010762',

    # http://www.informatics.jax.org/vocab/gene_ontology/GO:0006096 (glycolytic process / anaerobic glycolysis)
    'GO:0006096',

    # GO:1990349 (gap junction-mediated intercellular transport)
    # GO:0015139 (alpha-ketoglutarate transmembrane transporter activity)
    # Respirasome
    # Glycolysis
    # Motility / migration / locomotion / regulation of cell migration
    # Cell adhesion
}

glob = Path(__file__).parent.glob

sym2cat = pd.read_csv(unlist1(glob("../../data/*/*/goa_mouse_sym2cat.txt.gz")), sep='\t')

from tcga.utils import download

download = download.to(abs_path=unlist1(glob("../../data/*/Mouse-WCH*/*/download_cache")))

with download(URLS['obo']).now.open() as fd:
    G = obonet.read_obo(fd)
    assert isinstance(G, nx.MultiDiGraph)

for root in roots:
    # Take the subtree rooted at the `root` term
    # Keep only the genes implicated in those GO terms
    genes = sorted(set(sym2cat[sym2cat.goid.isin(G.subgraph(dfs_tree(G.reverse(), root)).nodes)].symbol))
    print(root, ":", genes)
