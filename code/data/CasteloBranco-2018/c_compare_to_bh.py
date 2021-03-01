# RA, 2021-04-21

from bugs import *
from tcga.utils import whatsmyname

from sklearn.metrics.pairwise import cosine_similarity

out_dir = mkdir(Path(__file__).with_suffix(''))

dir = next((p.resolve() for p in Path.cwd().parents for p in p.glob("**/Betsholtz-2018")), None)
df_expr_bh = pd.read_table(unlist1(dir.glob("*/*expr*")), index_col=0).T
df_meta_bh = pd.read_table(unlist1(dir.glob("*/*meta*")), index_col=0).reindex(df_expr_bh.columns)

dir = next((p.resolve() for p in Path.cwd().parents for p in p.glob("**/CasteloBranco-2018")), None)
df_expr_cb = pd.read_table(unlist1(dir.glob("*/*expr*")), index_col=0).T
df_meta_cb = pd.read_table(unlist1(dir.glob("*/*meta*")), index_col=0).reindex(df_expr_cb.columns)


# Only keep the genes in common and sort consistently
(df_expr_bh, df_expr_cb) = df_expr_bh.align(df_expr_cb, join="inner", axis=0)


def avg_cos_sim__mrkr():
    # This is from Betsholtz-2018
    # https://www.nature.com/articles/nature25739/figures/1
    # 1c
    markers = {
        'PC': "Pdgfrb Cspg4 Anpep Rgs5 Cd248 Abcc9 Vtn S1pr3",
        'SMC': "Acta2 Tagln Myh11 Myl9 Mylk Sncg Cnn1 Pln",
        'MG': "Csf1r Cd68 Cd53 Cd48 Cd84 C1qa Fcgr1 Emr1",
        'FB': "Pdgfra Lum Dcn Col3a1 Col5a1 Col8a2 Col12a1 Mmp2",
        'OL': "Mobp Plp1 Mog Cldn11 Mag Gjc2 Mal Cnp",
        'EC': "Pecam1 Kdr Flt1 Tie1 Tek Icam2 Podxl Ptprb",
        'AC': "Aldh1l1 Fgfr3 Slc4a4 Slc6a11 Slc7a10 Mlc1 Slc1a3 Cldn10",
    }

    markers = sorted(set(' '.join(markers.values()).split(' ')))

    # Convert to Samples x Marker genes
    bh = df_expr_bh.loc[markers].T
    cb = df_expr_cb.loc[markers].T
    assert set(bh.columns) == set(cb.columns)

    # # DEBUG
    # df_expr_cb = df_expr_cb.sample(n=9, random_state=43)
    # df_expr_ab = df_expr_ab.sample(n=8, random_state=43)

    # Cell label
    anno_cb = df_meta_cb.reindex(cb.index).celltype
    anno_bh = df_meta_bh.reindex(bh.index).celltype

    cb = cb.groupby(by=list(anno_cb)).mean()
    bh = bh.groupby(by=list(anno_bh)).mean()

    S = pd.DataFrame(index=cb.index, columns=bh.index, data=cosine_similarity(cb, bh))

    with Plox() as px:
        px.a.imshow(S, cmap="Greys")
        px.a.set_xticks(range(len(S.columns)))
        px.a.set_xticklabels(S.columns, rotation=90)
        px.a.set_xlabel("Betsholtz, 2018")
        px.a.set_yticks(range(len(S.index)))
        px.a.set_yticklabels(S.index)
        px.a.set_ylabel("Castelo-Branco, 2018")
        px.f.savefig(out_dir / f"{whatsmyname()}.png")


avg_cos_sim__mrkr()
