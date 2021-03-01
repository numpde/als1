# RA, 2021-03-16

"""
Prototype functions for the GO x Symbol incidence matrix.
"""

from bugs import *
from twig import log
from tcga.utils import from_iterable
from scipy.sparse import csr_matrix

log.warning("This script doesn't do anything.")

organism = "mouse"


def get_links(organism):
    file = unlist1(Path(__file__).parent.glob(f"a_*/goa_{organism}_sym2cat.txt.gz"))
    log.info(f"Loading file {relpath(file)}.")
    df = pd.read_csv(file, sep='\t')
    return df


def get_sparse_df(organism):
    df = get_links(organism)

    categories = df.symbol.groupby(df.cat).agg(list).sort_index()

    genes = pd.Series({g: i for (i, g) in enumerate(set(from_iterable(categories)))}).sort_index()

    log.info("Collecting the sparsity pattern.")

    (ii, jj) = np.array([
        (i, j)
        for (i, gg) in enumerate(categories)
        for j in genes[gg]
    ]).T

    log.info("Converting to a sparse dataframe.")

    df = pd.DataFrame.sparse.from_spmatrix(
        csr_matrix((np.ones_like(ii), (ii, jj)), shape=(len(categories), len(genes))),
        index=pd.Series(categories.index, name="goid"),
        columns=pd.Series(genes.index, name="gene"),
    )

    # 800MB CSV compresses to 1.5MB
    # Takes a while to write
    # df.sparse.to_dense().to_csv('test.csv.gz', compression='gzip')

    return df
