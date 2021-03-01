# RA, 2021-04-11

from bugs import *
from sklearn.decomposition import PCA

from numpy.random import default_rng

rng = default_rng(0)

C = np.array([[4, 0], [0, 1]])

[u, s, v] = np.linalg.svd(C)

N = 100
X = rng.multivariate_normal(mean=[0, 0], cov=C, size=N)

assert X.shape == (N, 2)

C1 = np.cov(X.T)
C2 = np.cov(PCA(n_components=2).fit_transform(X).T)

