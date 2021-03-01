# RA, 2021-04-21

import sys

from bugs import *
from tcga.utils import whatsmyname
from twig import log

from sklearn.metrics.pairwise import cosine_similarity

out_dir = mkdir(Path(__file__).with_suffix(''))

try:
    dataset_name = "CasteloBranco-2018"

    cb_dir = str(next((p.resolve() for p in Path.cwd().parents for p in p.glob(f"**/{dataset_name}")), None))
    (cb_dir in sys.path) or sys.path.append(cb_dir)

    from z_sources import df_expr as df_expr_cb, df_meta as df_meta_cb
except ImportError:
    log.warning("Import from z_sources failed.")
    raise

try:
    ab_dir = next((p.resolve() for p in Path.cwd().parents for p in p.glob("**/Mouse-WCH-2020")), None)
    df_expr_ab = pd.read_table(unlist1(ab_dir.glob("*fewer_cells/*data*")), index_col=0)
    df_meta_ab = pd.read_table(unlist1(ab_dir.glob("*fewer_cells/*meta*")), index_col=0)
    df_meta_ab = df_meta_ab.reindex(df_expr_ab.columns)
except Exception:
    raise


# Only keep the genes in common and sort consistently
(df_expr_ab, df_expr_cb) = df_expr_ab.align(df_expr_cb.T, join="inner", axis=0)


def avg_cos_sim__mrkr():
    # This is from Betsholtz-2018
    # https://www.nature.com/articles/nature25739/figures/1
    # 1c
    markers = {
        'PC': "Pdgfrb Cspg4 Anpep Rgs5 Cd248 Abcc9 Vtn S1pr3",
        'SMC': "Acta2 Tagln Myh11 Myl9 Mylk Sncg Cnn1 Pln",
        'MG': "Csf1r Cd68 Cd53 Cd48 Cd84 C1qa Fcgr1",  # removed Emr1 (missing in AB)
        'FB': "Pdgfra Lum Dcn Col3a1 Col5a1 Col8a2 Col12a1 Mmp2",
        'OL': "Mobp Plp1 Mog Cldn11 Mag Gjc2 Mal Cnp",
        'EC': "Pecam1 Kdr Flt1 Tie1 Tek Icam2 Podxl Ptprb",
        'AC': "Aldh1l1 Fgfr3 Slc4a4 Slc6a11 Slc7a10 Mlc1 Slc1a3 Cldn10",
    }

    markers = sorted(set(' '.join(markers.values()).split(' ')))

    ab = df_expr_ab.loc[markers].T
    cb = df_expr_cb.loc[markers].T
    assert set(ab.columns) == set(cb.columns)

    # # DEBUG
    # df_expr_cb = df_expr_cb.sample(n=9, random_state=43)
    # df_expr_ab = df_expr_ab.sample(n=8, random_state=43)

    # Cell label
    anno_cb = df_meta_cb.reindex(cb.index).celltype
    anno_ab = df_meta_ab.reindex(ab.index).cell_type_alias_label

    cb = cb.groupby(by=list(anno_cb)).mean()
    ab = ab.groupby(by=list(anno_ab)).mean()

    S = pd.DataFrame(index=cb.index, columns=ab.index, data=cosine_similarity(cb, ab))

    with Plox() as px:
        px.a.imshow(S, cmap="Greys")
        px.a.set_xticks(range(len(S.columns)))
        px.a.set_xticklabels(S.columns, rotation=90)
        px.a.set_xlabel("Mouse-WCH-2020")
        px.a.set_yticks(range(len(S.index)))
        px.a.set_yticklabels(S.index)
        px.a.set_ylabel(dataset_name)
        px.f.savefig(out_dir / f"{whatsmyname()}.png")


avg_cos_sim__mrkr()
