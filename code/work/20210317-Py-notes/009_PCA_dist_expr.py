# RA, 2021-04-14

"""
PCA + multivariate normalization + TSNE
"""

from bugs import *

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import RobustScaler

src_dir = next((p for p in Path(__file__).parents for p in p.glob("**/Mouse-WCH-2020")), None)
out_dir = mkdir(Path(__file__).with_suffix(''))

df_meta = pd.read_table(unlist1(src_dir.glob("*fewer_cells/meta.csv.gz")), index_col=0)
df_expr = pd.read_table(unlist1(src_dir.glob("*fewer_cells/data.csv.gz")), index_col=0).T

# Note: np.cov(prepro.fit_transform(data).T) is identity
prepro = make_pipeline(PCA(n_components=30, random_state=43), RobustScaler(quantile_range=(0.05, 0.95)))
prepro.fit(df_expr)

X = pd.DataFrame(index=df_expr.index, data=prepro.transform(df_expr))

with Plox() as px:
    px.a.imshow(np.cov(X.T))
    px.f.savefig(out_dir / "cov.png")

with Plox() as px:
    px.a.plot(X[0], X[1], '.', ms=1, alpha=0.3)
    px.f.savefig(out_dir / "scatter.png")

X2 = pd.DataFrame(index=X.index, columns=["x", "y"], data=TSNE(random_state=43).fit_transform(X))

with Plox() as px:
    for (label, df) in X2.groupby(df_meta.subclass_label):
        px.a.plot(df.x, df.y, '.', ms=2, label=label)
    px.a.axis('off')
    px.a.legend()
    px.f.savefig(out_dir / "tsne.png")
