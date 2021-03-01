# RA, 2021-03-16

"""
Prototype for relating and visualizing GO categories.

TODO:
https://stackoverflow.com/questions/29586520/can-one-get-hierarchical-graphs-from-networkx-with-python-3
"""

# central nervous system development
"http://www.informatics.jax.org/vocab/gene_ontology/GO:0007417"
goid = "GO:0007417"

import obonet
import networkx as nx
from bugs import *
from twig import log
from plox import rcParam
from tcga.utils import from_iterable
from z_sources import URLS, download, url2file

out_dir = mkdir(Path(__file__).with_suffix(''))

sym_cat = pd.read_csv(url2file(URLS.goa.mouse), sep='\t')

"""
sym_cat is the dataframe:

               symbol        goid
0       0610009B22Rik  GO:0003674
1       0610009B22Rik  GO:0006888
...
296655            rp9  GO:0016607
"""

with download(URLS.obo).now.open() as fd:
    G = obonet.read_obo(fd)
    assert isinstance(G, nx.MultiDiGraph)

log.info(f"GO term {goid}")
for (k, v) in G[goid].items():
    log.info(f"{unlist1(v)}: {k}")

nlevels = 2

hierarchy = pd.DataFrame({'goid': [goid], 'level': 0}).set_index('goid').level
for level in range(nlevels):
    hierarchy = pd.concat(objs=(
        hierarchy,
        pd.DataFrame(
            from_iterable(
                # out_edges: hierarchy up; in_edges: hierarchy down
                [{'goid': a, 'level': level + 1} for (a, b, d) in G.in_edges(goid, data=True)]
                for goid in hierarchy[hierarchy == level].index
            )
        ).set_index('goid').level
    ))

# print(list(hierarchy.index))
# exit()

cat_name = nx.get_node_attributes(G, name='name')

g = G.subgraph(hierarchy.index)

# Genes x GOID incidence matrix
I = sym_cat[sym_cat.goid.isin(g.nodes)].assign(v=1).pivot(index="symbol", columns="goid", values='v')
I = I.T.reindex(hierarchy.index).T
I = I.fillna(0).astype(int)
log.info(I)

# Some g.nodes are not represented in sym_cat
# g = g.subgraph()

style = {
    rcParam.Xtick.labelsize: 1,
    rcParam.Ytick.labelsize: 1,
    rcParam.Figure.figsize: (12, 6),  # no effect?
    rcParam.Axes.linewidth: 0.1,
    rcParam.Xtick.Major.size: 0,
    rcParam.Ytick.Major.size: 0,
}

with Plox(style) as px:
    i = I.T  # goid x gene
    i = i.apply(lambda s: s * (1 + hierarchy), axis=0)
    px.a.imshow(i)
    px.a.set_xticks(range(len(i.columns)))
    px.a.set_xticklabels(list(i.columns), rotation=90)
    px.a.set_yticks(range(len(i.index)))
    px.a.set_yticklabels(list(i.index))
    px.f.savefig(out_dir / f"incidence_{goid}.png", dpi=600)

with Plox() as px:
    pos = nx.drawing.nx_agraph.graphviz_layout(g, prog="dot", root=goid)
    nx.draw(g, ax=px.a, pos=pos, with_labels=True, arrows=False, font_size=3, node_size=10)
    px.a.invert_yaxis()
    px.f.savefig(out_dir / f"graph_{goid}.png", dpi=300)
