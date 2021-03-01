# RA, 2021-03-18

from bugs import *
from twig import log
from tcga.utils import First
from plox import rcParam
from b_poc2 import Setup
from sklearn.manifold import TSNE

src_file = max(Path(__file__).parent.glob("b_poc2/run*/*/latent*"))
log.info(f"Reading from {relpath(src_file)}.")

df_latent = pd.read_table(src_file, index_col=0)
log.info(f"df_latent: {len(df_latent)} samples x {len(df_latent.columns)} features.")

# go_latent = df_latent.T.apply(lambda s: pd.Series(s).nlargest(3).index).T[0]

read = First(Setup.path_to_reduced.glob).then(unlist1).then(lambda f: pd.read_table(f, index_col=0, low_memory=False))

df_meta = read("meta*").reindex(df_latent.index)

## Plot

# GROUP_SAMPLES_BY = 'cell_type_alias_label'
# GROUP_SAMPLES_BY = 'class_label'
GROUP_SAMPLES_BY = 'subclass_label'

# df: pd.DataFrame
# df = pd.DataFrame({
#     'celltype': df_meta.subclass_label,
#     'golatent': go_latent,
# }).pivot_table(
#     index='golatent', columns='celltype', aggfunc=len, fill_value=1
# )

df: pd.DataFrame
df = df_latent.groupby(df_meta[GROUP_SAMPLES_BY]).sum().T
df = df.apply((lambda s: s / (s.sum() or 1)), axis=0)

ii = df.sum(axis=1).nlargest(55).index
df = df[df.index.isin(ii)]  # .T.assign(other=df[~df.index.isin(ii)].sum(axis=0)).T
df = df.sort_index(axis=1)

go_cat_size: pd.Series
go_cat_size = Setup.sym2cat.groupby('goid').symbol.count()

df = df.reindex(go_cat_size.reindex(df.index).sort_values().index)

style = {
    rcParam.Xtick.labelsize: 3,
    rcParam.Ytick.labelsize: 2,
    rcParam.Axes.linewidth: 0.1,
    rcParam.Figure.figsize: (df.shape[1] / 15, df.shape[0] / 15),
}

with Plox(style) as px:
    px.a.imshow(df, cmap='Reds', zorder=1)

    px.a.set_ylabel("Latent GO category")
    px.a.set_yticks(range(len(df.index)))
    px.a.set_yticklabels(labels=[
        f"{i} ({n})"
        for (i, n) in go_cat_size.reindex(df.index).iteritems()
    ])
    # px.a.axis("left").major_ticklabels.set_ha("left")
    # [t.set_ha("left") for t in px.a.yaxis.get_majorticklabels()]
    for (n, name) in enumerate(Setup.get_go_terms().reindex(df.index)):
        px.a.text(x=0, y=n, s=name[0:60], fontdict={'fontsize': 3}, color='b', va='center', alpha=0.4, zorder=1000)

    px.a.set_xlabel("Cell subtype")
    px.a.set_xticks(range(len(df.columns)))
    px.a.set_xticklabels(rotation=90, labels=[
        f"{n} x {i}"
        for (i, n) in df_meta[GROUP_SAMPLES_BY].value_counts().reindex(df.columns).iteritems()
    ])

    # Grid
    for x in np.linspace(min(px.a.get_xlim()), max(px.a.get_xlim()), 1 + len(px.a.get_xticks())):
        px.a.plot([x, x], px.a.get_ylim(), 'k', lw=0.01, zorder=10)
    for y in np.linspace(min(px.a.get_ylim()), max(px.a.get_ylim()), 1 + len(px.a.get_yticks())):
        px.a.plot(px.a.get_xlim(), [y, y], 'k', lw=0.01, zorder=10)

    px.f.savefig(src_file.parent / f"incidence.png", dpi=300)

exit()

for (reduced, df) in enumerate([df_expr, df_latent]):
    X2 = pd.DataFrame(
        data=TSNE(n_components=2, random_state=43).fit_transform(df),
        index=df.index, columns=["x", "y"],
    )

    import matplotlib.colors as mcolors

    # https://matplotlib.org/stable/tutorials/intermediate/color_cycle.html
    from cycler import cycler

    style = {
        rcParam.Legend.fontsize: 3,
        rcParam.Axes.prop_cycle: cycler(color=mcolors.XKCD_COLORS),
    }

    with Plox(style) as px:
        for (label, x2) in X2.groupby(df_meta[GROUP_SAMPLES_BY]):
            px.a.scatter(x2.x, x2.y, label=label, s=1)
        px.a.set_axis_off()
        px.a.legend()
        px.f.savefig(src_file.parent / f"tsne__celltype__reduced={reduced}.png")

    with Plox(style) as px:
        for (goid, x2) in X2.groupby(go_latent):
            px.a.scatter(x2.x, x2.y, label=f"{goid} ({len(x2)})", s=1)
        px.a.set_axis_off()
        px.a.legend()
        px.f.savefig(src_file.parent / f"tsne__latentgo__reduced={reduced}.png")
