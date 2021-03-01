# RA, 2021-04-11

from bugs import *
from sklearn.decomposition import PCA

from numpy.random import default_rng
rng = default_rng(0)

n_cells = 33
n_genes = 20

X = rng.integers(low=0, high=20, size=(n_cells, n_genes))
X = X * np.atleast_2d(np.linspace(1, 1000, X.shape[1]))  # gene mean
# X = X * np.atleast_2d(np.linspace(1, 3, X.shape[0])).T  # library size

pca = PCA(n_components=5, random_state=0)
pca.fit(X)
x = pca.transform(X)

pca.explained_variance_

pca.inverse_transform(pca.fit_transform(X))

# Covariance-distance

C = (pca.get_covariance())
D = (lambda x, y: np.dot(x - y, C @ (x - y)))

M = np.array([[D(a, b) for a in X] for b in X])

# Covariance-distance in reduced space

c = np.cov(pca.transform(X).T)
assert c.shape == (pca.n_components, pca.n_components)
d = (lambda x, y: np.dot(x - y, c @ (x - y)))

m = np.array([[d(a, b) for a in x] for b in x])

#

with Plox() as px:
    px.a.imshow(M)
    px.a.axis('off')
    px.f.savefig(mkdir(Path(__file__).with_suffix('')) / "M.png")


with Plox() as px:
    px.a.imshow(m)
    px.a.axis('off')
    px.f.savefig(mkdir(Path(__file__).with_suffix('')) / "m.png")


#

# # loss prototype:
# p = np.linalg.inv(np.cov(pca.transform(X).T))
# d = pca.transform(X) - pca.transform(X * 1.01)
# np.sum((d @ p) * d, axis=1)

#

pca2 = PCA(n_components=2, random_state=0)
pca2.fit(X)
X2 = pca2.transform(X).astype(float)

assert pca2.get_covariance().shape == (n_genes, n_genes)

c2 = PCA(n_components=pca2.n_components, random_state=0).fit(X2).get_covariance()
[u, s, v] = np.linalg.svd(c2)
# u @ np.diag(s) @ v == c2
z2 = X2.mean(axis=0)

# tr = (lambda a, b: )

with Plox() as px:
    px.a.scatter(*X2.T, s=4)
    px.a.axis('off')
    px.a.plot(*np.array([z2, z2 + u[:, 0] * np.sqrt(s[0])]).T)
    px.a.plot(*np.array([z2, z2 + u[:, 1] * np.sqrt(s[1])]).T)
    # px.a.quiver(z2, z2 + u[:, 1])
    px.f.savefig(mkdir(Path(__file__).with_suffix('')) / "X2.png")

#

# Z = rng.random(size=(20, 30000))
# pca = PCA(n_components=10).fit(Z)
# print(np.cov(pca.transform(Z).T))
#
