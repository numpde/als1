# RA, 2021-03-21

import numpy as np
from scipy.stats import nbinom
from numpy.random import default_rng

n = 10  # aka r
p = 0.3

# Note: p in scipy is (1 - p) on wiki
mu = n * (1 - p) / p
var = n * (1 - p) / (p ** 2)

D = nbinom(n=n, p=p)
# D = rng.negative_binomial(n=n, p=p)

x = np.arange(1, 100)
f = D.pmf(x)

print("Mean:", [np.dot(x, f), mu, D.mean()])
print("Var :", [var, D.var(), mu / p])

# Mean-variance relation

# Dispersion parameterization
# var = mu + mu ** 2 / phi
phi = mu / (var / mu - 1)
print("Dispersion: ", [phi, n])
# n
print("n = Mean^2 / (Var - Mean):", [mu / (var / mu - 1), n])
# This is the scipy p
print("scipy p = Mean / variance:", [mu / var, p, (1 / (1 + mu / phi))])
print("wiki p = 1 - Mean / variance:", [(var / mu - 1) * mu / var, 1 - p])

# Numpy negative_binomial
rng = default_rng(seed=43)
x = rng.negative_binomial(n=n, p=p, size=10000)
print("numpy mean, var:", [x.mean(), x.var()])


def gen(phi, ngenes=11, nsamples=7, seed=43):
    import pandas
    from numpy.random import default_rng
    from scipy.stats import nbinom
    rng = default_rng(seed=seed)
    means = rng.lognormal(mean=2, sigma=1, size=ngenes)
    gg = pandas.Index(range(ngenes), name="gene")
    ii = pandas.Index(range(nsamples), name="sample")
    expr = [
        rng.negative_binomial(n=phi, p=(1 / (1 + mu / phi)), size=nsamples)
        for mu in means
    ]
    return pandas.DataFrame(index=gg, columns=ii, data=expr)



df = gen(phi=phi, ngenes=1000, nsamples=1000)
M = df.mean(axis=1)
V = df.var(axis=1)

import statsmodels.api as sm
phi_hat1 = np.exp(np.mean(2 * np.log(M)) - np.mean(np.log(np.abs(V - M))))
phi_hat2 = sm.OLS(M ** 2, V - M).fit().params.x1
print("phi:", phi, "and phi estimate:", [phi_hat1, phi_hat2])

# from statsmodels.genmod.families.family import NegativeBinomial
# import statsmodels.formula.api as smf
# import pandas as pd
# df_meanvar = pd.DataFrame({'M': M, 'V': V}).sort_values(by='M')
# # https://campus.datacamp.com/courses/generalized-linear-models-in-python/modeling-count-data?ex=12
# nb = smf.glm(formula='V ~ M + I(M*M) - 1', data=df_meanvar, family=sm.families.Gaussian()).fit()
# print(nb.summary()) print(nb.params) print(nb.bse) print(nb.pvalues)
# np.exp(nb.params.Intercept)

from plox import Plox
from pathlib import Path
with Plox() as px:
    # Fraction of zeros
    px.a.plot(M, (df == 0).mean(axis=1), '.', label="observed")
    # Expected fraction of zeros
    px.a.plot(M.sort_values(), (1 / (1 + M.sort_values() / phi_hat1)) ** phi_hat1, '--', label="expected")
    px.a.set_xscale('log')
    px.a.set_ylabel("Fraction of zeros")
    px.a.set_xlabel("Mean expression")
    px.a.legend()
    px.f.savefig(Path(__file__).with_suffix('.png'))
