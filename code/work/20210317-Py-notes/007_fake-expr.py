# RA, 2021-03-21

import numpy as np
import pandas as pd
from plox import Plox
from numpy.random import default_rng

rng = default_rng(seed=0)

ngenes = 2000
nsamples = 1000

phi = 0.5

# target means from f_fewer_cells.stats.py
Emu = 0.08927231700286073
Vmu = 2.2484780351522575
# Gamma distro:
# https://en.wikipedia.org/wiki/Gamma_distribution
# mean = k * theta, var = k * theta^2
# In numpy: k = shape, theta = scale
# => theta = var / mean, k = mean / theta = mean^2 / var

# theta = Vmu / Emu
# k = Emu / theta

# target_means = rng.gamma(shape=((Emu ** 2) / Vmu), scale=(Vmu / Emu), size=ngenes)
target_means = rng.gamma(shape=0.1, scale=0.1, size=ngenes)

# with Plox() as px:
#     [hh, ee] = np.histogram(np.log(target_means), bins=20)  # 'scott')
#     ee = np.exp(ee)
#     px.a.bar(ee[0:-1], hh, width=np.diff(ee), edgecolor='w', align='edge')
#     px.a.set_yscale('log')
#     px.a.set_xscale('log')
#     px.show()
#     exit()

df_expr = pd.DataFrame(
    index=pd.Index(range(ngenes), name="gene"),
    columns=pd.Index(range(nsamples), name="sample"),
    data=[
        rng.negative_binomial(n=phi, p=(1 / (1 + mu / phi)), size=nsamples)
        for mu in target_means
    ],
    dtype=int,
)


M = df_expr.mean(axis=1)
V = df_expr.var(axis=1)

with Plox() as px:
    px.a.hist(df_expr.sum(axis=0))
    px.a.set_xlabel("Library size")
    px.a.set_ylabel("Number of samples")
    px.show()

exit()

# Mean-variance
with Plox() as px:
    px.a.plot(M, V, '.')
    px.a.set_xlabel("Mean")
    px.a.set_xlabel("Variance")
    px.a.set_xscale('log')
    px.a.set_yscale('log')
    px.show()

with Plox() as px:
    [hh, ee] = np.histogram(np.log(M[M > 0]), bins=20)  # 'scott')
    ee = np.exp(ee)
    px.a.bar(ee[0:-1], hh, width=np.diff(ee), edgecolor='w', align='edge')
    px.a.set_yscale('log')
    px.a.set_xscale('log')
    # px.a.set_ylabel(rf"Out of {sum(df_meanvar.M.dropna() > 0)} nonzero genes, how many ...")
    # px.a.set_xlabel(rf"... have this mean expression across samples")
    px.show()
