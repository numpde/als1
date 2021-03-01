# RA, 2021-03-15

# Based on
# https://github.com/numpde/as/tree/master/p/single-cell/20180118-GO

from bugs import *

import networkx as nx
import obonet

from tcga.utils import whatsmyname, from_iterable, pretty

from z_sources import URLS, download

out_dir = Path(__file__).with_suffix('')


def read_goa(data: pd.DataFrame):
    """
    Read the GO annotations table
    Perform some consistency checks
    Returns the tuple
      (S0, G, H, S01)
    where
      S0 is the set of primary symbols
    G, H are bipartite graphs:
      H of primary symbols <--> synonyms
      G of primary symbols <--> GO terms
    and
      S01 is the set of synonyms that appear as primary symbols
      (they are omitted from the graphs)
    """

    from itertools import chain, groupby
    from operator import itemgetter

    # Check that each ID has at exactly one (primary) symbol
    # (but some symbols have multiple IDs)
    g = nx.Graph()
    g.add_edges_from(('ID:' + i, 'Sy:' + s) for (i, s) in zip(data['ID'], data['Symbol']))
    assert (all((d == 1) for (n, d) in g.degree() if n.startswith('ID:')))
    del g

    # Obtain the symbol synonyms (excluding itself)
    SYN = [
        (sym, set(str(syn).split('|')) - set(sym))
        for (sym, syn) in zip(data['Symbol'], data['Synonyms'])
    ]

    # Check that each symbol has a unique list of synonyms
    assert (all(
        (1 == len(set(tuple(sorted(s)) for (_, s) in S)))
        for (_, S) in groupby(SYN, itemgetter(0))
    ))

    # Now can collapse the list to a dict
    SYN = dict(SYN)

    # Primary symbols
    S0 = set(SYN.keys())

    # Primary symbols that appear as synonyms
    S01 = S0 & set(chain.from_iterable(SYN.values()))

    # Remove the primary symbols from the synonyms
    SYN = {s0: (syn - S01) for (s0, syn) in SYN.items()}

    # All synonyms (excluding the primary symbols)
    S1 = set(chain.from_iterable(SYN.values()))

    # They should be disjoint by now:
    assert (not (S0 & S1))

    # For bipartite graphs, see
    # https://networkx.github.io/documentation/stable/reference/algorithms/bipartite.html

    # Make a bipartite graph "ID" <--> "Synonyms"
    H = nx.Graph()
    H.add_nodes_from(S0, is_symbol=1, is_synonym=0)
    H.add_nodes_from(S1, is_symbol=0, is_synonym=1)
    H.add_edges_from((i, s) for (i, syn) in SYN.items() for s in syn)
    # There should be no self-loops
    assert (all((a != b) for (a, b) in H.edges()))
    # Each synonym belongs to at least one primary symbol
    assert (all(H.degree(s) for s in S1))

    # Make a bipartite graph "Primary gene symbols" <--> "GO terms"
    G = nx.Graph()
    assert (S0 == set(data['Symbol']))
    G.add_nodes_from(data['Symbol'], is_symbol=1, is_go=0)
    G.add_nodes_from(data['GO'], is_symbol=0, is_go=1)
    G.add_edges_from(zip(data['Symbol'], data['GO']))

    return (S0, G, H, S01)


def goa(url):
    import re
    organism = unlist1(re.findall(pattern=r"go/goa/([a-z]+)/", string=str.lower(url)))

    with download(url).now.open(mode='rb') as gz:
        # Use pandas to load the GO annotations table
        # http://www.geneontology.org/page/go-annotation-file-format-20
        columns = [
            'DB', 'ID', 'Symbol', 'Q', 'GO', 'DB Ref', 'Evidence', 'With/From', 'Aspect',
            'Name', 'Synonyms', 'Type', 'Taxon', 'Date', 'Assigned by', 'Extension', 'Gene product ID',
        ]

        df = pd.read_table(gz, compression='gzip', index_col=False, comment='!', sep='\t', header=None, names=columns)

        (S, G, H, _) = read_goa(df)

    # Symbols and synonyms

    filename = out_dir / f"{whatsmyname()}_{organism}_sym2syn.txt.gz"

    pd.DataFrame(
        columns=["symbol", "synonym", "is_ambiguous"],
        data=from_iterable([
            [
                # Ambiguous synonyms
                (symbol, syn, True)
                for syn in sorted(n for n in H.neighbors(symbol) if (H.degree(n) >= 2))
            ] + [
                # Non-ambiguous synonyms
                (symbol, syn, False)
                for syn in sorted(n for n in H.neighbors(symbol) if (H.degree(n) == 1))
            ]
            for symbol in sorted(S)
        ])
    ).to_csv(
        filename, sep='\t', index=False, compression='gzip',
    )

    with filename.with_suffix(".readme.txt").open(mode='w') as fd:
        print(pretty(download(url).now.meta), file=fd)

    # GO terms

    filename = out_dir / f"{whatsmyname()}_{organism}_sym2cat.txt.gz"

    pd.DataFrame(
        columns=["symbol", "goid"],
        data=((symbol, cat) for symbol in sorted(S) for cat in G.neighbors(symbol)),
    ).to_csv(
        filename, sep='\t', index=False, compression='gzip',
    )

    with filename.with_suffix(".readme.txt").open(mode='w') as fd:
        print(pretty(download(url).now.meta), file=fd)


def obo(url):
    with download(url).now.open(mode='r') as fd:
        # GO graph
        G = obonet.read_obo(fd)
        assert isinstance(G, nx.MultiDiGraph)

    # pickle.dump(G, open(OFILE['graph'], 'wb'))

    # GO names

    filename = out_dir / f"{whatsmyname()}_{G.name}.txt.gz"

    pd.DataFrame(
        columns=["goid", "term", "namespace"],
        data=(
            (n, G.nodes[n]['name'], G.nodes[n]['namespace'])
            for n in sorted(G.nodes)
        )
    ).to_csv(
        filename, sep='\t', index=False, compression='gzip',
    )

    with filename.with_suffix(".readme.txt").open(mode='w') as fd:
        print(pretty(download(url).now.meta), file=fd)


if __name__ == '__main__':
    goa(URLS.goa.mouse)
    obo(URLS.obo)
