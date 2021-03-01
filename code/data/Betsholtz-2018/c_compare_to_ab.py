# RA, 2021-04-21

import sys

from bugs import *
from tcga.utils import whatsmyname
from twig import log

from scipy.stats.kde import gaussian_kde
from sklearn.metrics.pairwise import cosine_similarity

import matplotlib.pyplot as plt

out_dir = mkdir(Path(__file__).with_suffix(''))

try:
    bh_dir = str(next((p.resolve() for p in Path(__file__).parents for p in p.glob("**/Betsholtz-2018")), None))
    (bh_dir in sys.path) or sys.path.append(bh_dir)

    from z_sources import df_expr as df_expr_bh, df_meta as df_meta_bh, df_mrkr as mrkr, markers
except ImportError:
    log.warning("Import from z_sources failed.")
    raise

try:
    ab_dir = next((p.resolve() for p in Path(__file__).parents for p in p.glob("**/Mouse-WCH-2020")), None)
    df_expr_ab = pd.read_table(unlist1(ab_dir.glob("*fewer_cells/*data*")), index_col=0).T
    df_meta_ab = pd.read_table(unlist1(ab_dir.glob("*fewer_cells/*meta*")), index_col=0)
    df_meta_ab = df_meta_ab.reindex(df_expr_ab.index)
except Exception:
    raise

# Only keep the genes in common and sort consistently
(df_expr_ab, df_expr_bh) = df_expr_ab.align(df_expr_bh, join="inner", axis=1)


def histograms():
    genes = [
        "Trpm3", "mt-Co1", "mt-Co3", "Nnat", "Ptgds", "Adam12", "Alcam",  # Up early
        "Itih5", "Malat1", "Zbtb20", "Spp1", "Col15a1", "Ece1", "Cemip",  # Up late
    ]

    ab = df_expr_ab.reindex(genes, axis=1).dropna(axis=1)
    bh = df_expr_bh.reindex(genes, axis=1).dropna(axis=1)
    assert list(ab.columns) == list(bh.columns)

    # Norm1 normalization
    # ab = ab.div(ab.sum(axis=1), axis=0)
    # bh = bh.div(bh.sum(axis=1), axis=0)

    # log1p trafo
    ab = ab.transform(lambda x: np.log(x + 1))
    bh = bh.transform(lambda x: np.log(x + 1))

    ab = ab[df_meta_ab.subclass_label == "VLMC"]
    bh = bh[df_meta_bh.celltype.str.startswith("FB")]
    assert list(ab.columns) == list(bh.columns)

    genes = sorted(ab.columns)

    # g = first(genes)

    colors = {
        'FB1': 'purple',
        'FB2': 'violet',
        '374_VLMC': "C0",
        '375_VLMC': "C1",
        '376_VLMC': "C2",
    }

    for g in genes:
        expr: pd.Series
        with Plox() as px:
            grps = [
                ab[g].groupby(df_meta_ab.reindex(ab.index).cell_type_alias_label),
                bh[g].groupby(df_meta_bh.reindex(bh.index).celltype),
            ]

            for grp in grps:
                for (label, expr) in grp:
                    if any(expr):
                        f = gaussian_kde(expr)
                        xx = np.linspace(0, max(expr) * 1.5, 100)
                        px.a.plot(xx, f(xx), label=f"{label} ({len(expr)})", color=colors[label])

            px.a.set_xlabel("log1p(count)")
            px.a.legend()
            px.a.set_yticks([])
            px.f.savefig(out_dir / f"hist_{g}.png")


def allvlmc_vs_allfb__mrkr():
    ab = df_expr_ab.reindex(mrkr.index, axis=1).dropna(axis=1)
    bh = df_expr_bh.reindex(mrkr.index, axis=1).dropna(axis=1)
    assert set(ab.columns) == set(bh.columns)

    ab = df_expr_ab[df_meta_ab.subclass_label == "VLMC"]
    bh = df_expr_bh[df_meta_bh.celltype.str.startswith("FB")]

    ab = ab.reindex(df_meta_ab.reindex(ab.index).cell_type_alias_label.sort_values().index)
    bh = bh.reindex(df_meta_bh.reindex(bh.index).celltype.sort_values().index)

    S = pd.DataFrame(index=bh.index, columns=ab.index, data=cosine_similarity(bh, ab))

    njj = ab.index.groupby(df_meta_ab.reindex(ab.index).cell_type_alias_label.to_list())
    mii = bh.index.groupby(df_meta_bh.reindex(bh.index).celltype.to_list())

    fig: plt.Figure
    (fig, AX) = plt.subplots(len(mii), len(njj))

    for (i, (m, ii)) in enumerate(mii.items()):
        for (j, (n, jj)) in enumerate(njj.items()):
            ax = AX[i, j]
            ax.imshow(S.loc[ii, jj], aspect='auto', cmap='Greys')
            ax.axis('off')

    # No effect?
    plt.xlabel("Mouse-WCH-2020")
    plt.ylabel("Betsholtz-2018")

    fig.savefig(out_dir / f"{whatsmyname()}.png")


def avg_cos_sim__mrkr():
    ab = df_expr_ab.reindex(mrkr.index, axis=1).dropna(axis=1)
    bh = df_expr_bh.reindex(mrkr.index, axis=1).dropna(axis=1)
    assert set(ab.columns) == set(bh.columns)

    # # DEBUG
    # df_expr_bh = df_expr_bh.sample(n=9, random_state=43)
    # df_expr_ab = df_expr_ab.sample(n=8, random_state=43)

    # Cell label
    anno_bh = df_meta_bh.reindex(bh.index).celltype
    anno_ab = df_meta_ab.reindex(ab.index).cell_type_alias_label

    bh = bh.groupby(by=list(anno_bh)).mean()
    ab = ab.groupby(by=list(anno_ab)).mean()

    S = pd.DataFrame(index=bh.index, columns=ab.index, data=cosine_similarity(bh, ab))

    with Plox() as px:
        px.a.imshow(S, cmap="Greys")
        px.a.set_xticks(range(len(S.columns)))
        px.a.set_xticklabels(S.columns, rotation=90)
        px.a.set_xlabel("Mouse-WCH-2020")
        px.a.set_yticks(range(len(S.index)))
        px.a.set_yticklabels(S.index)
        px.a.set_ylabel("Betsholtz-2018")
        px.f.savefig(out_dir / f"{whatsmyname()}.png")


histograms()
allvlmc_vs_allfb__mrkr()
avg_cos_sim__mrkr()
