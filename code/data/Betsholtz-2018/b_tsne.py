# RA, 2021-04-14

from bugs import *
from twig import log

from sklearn.manifold import TSNE

out_dir = mkdir(Path(__file__).with_suffix(''))

# df_expr = pd.DataFrame(np.random.RandomState(0).random(size=(30, 5)))

from z_sources import df_expr, df_meta, df_mrkr

# # DEBUG -- subset samples
# df_expr = df_expr.sample(n=100, random_state=43)

# # Subset to marker genes
# df_expr = df_expr[df_mrkr.index]

# Order samples
df_meta = df_meta.reindex(df_expr.index)

styles = pd.DataFrame(index=['marker', 'c'], data={
    'PC': ('s', 'red'),
    'vSMC': ('s', 'green'),
    'aaSMC': ('o', 'green'),
    'aSMC': ('^', 'green'),
    'MG': ('s', 'gray'),
    'FB1': ('^', 'purple'),
    'FB2': ('s', 'violet'),
    'OL': ('s', 'brown'),
    'EC1': ('s', 'cyan'),
    'EC2': ('o', 'cyan'),
    'EC3': ('^', 'cyan'),
    'vEC': ('s', 'blue'),
    'capilEC': ('o', 'blue'),
    'aEC': ('^', 'blue'),
    'AC': ('s', 'orange'),
})


def cluster(df_expr):
    df = pd.DataFrame(index=df_expr.index, columns=["x", "y"], data=TSNE(random_state=43).fit_transform(df_expr))
    df = df.assign(celltype=list(df_meta.reindex(df.index).celltype))
    return df


# Full set of samples

with Plox() as px:
    for (label, X) in cluster(df_expr).groupby(by='celltype'):
        px.a.scatter(X.x, X.y, s=2, alpha=0.7, label=label, **styles[label].to_dict())
    px.a.legend(loc='lower right', fontsize=4)
    px.a.axis('off')
    px.f.savefig(out_dir / "tsne.png")

# Subset to fibroblasts

with Plox() as px:
    for (label, X) in cluster(df_expr.loc[df_meta[df_meta.celltype.isin(['FB1', 'FB2'])].index]).groupby(by='celltype'):
        px.a.scatter(X.x, X.y, s=2, alpha=0.7, label=label, **styles[label].to_dict())
    px.a.legend(loc='lower right', fontsize=4)
    px.a.axis('off')
    px.f.savefig(out_dir / "tsne_fb_only.png")
