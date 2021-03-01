# RA, 2021-03-14

import torch
import torch.nn as nn

from bugs import *
from numpy.random import default_rng

n = 20  # number of genes
k = 5  # number of categories
rng = default_rng(seed=0)
cat_mask = pd.DataFrame(
    columns=pd.Series(name="category", data=range(k)),
    index=pd.Series(name="gene", data=range(n)),
    data=rng.choice([True, False], size=[n, k], p=(lambda p: (p, 1 - p))(2 * k / n)),
).astype(int)

print(cat_mask)

# Dream up per-category expression profiles
profiles = cat_mask * rng.integers(low=0, high=4, size=cat_mask.shape)
profiles = profiles / profiles.sum(axis=0)

print(profiles)

# Reconstruction of category mixture with NNLS
# provided the `profiles` are known
z = pd.Series(index=profiles.columns, data=rng.random(size=k))
x = profiles @ z
from scipy.optimize import nnls

print(pd.DataFrame({'[!]': z, '[?]': first(nnls(profiles, x))}))

# Now suppose the profiles are not known but we have N samples
N = 19
U = pd.DataFrame(index=profiles.columns, columns=pd.Series(range(N), name="samples"), data=rng.random(size=[k, N]))
X = profiles @ U
print(U)

mask = (lambda t: t * torch.tensor(cat_mask.values, requires_grad=False))


# The following is based on:
# https://pmelchior.net/blog/proximal-matrix-factorization-in-pytorch.html


class MF(nn.Module):
    def __init__(self):
        super().__init__()
        self.p = nn.Parameter(torch.rand(size=(n, k), requires_grad=True))
        self.u = nn.Parameter(torch.rand(size=(k, N), requires_grad=True))

    def forward(self):
        return self.p @ self.u


mf = MF()

for i in range(1000):
    loss = (mf() - torch.tensor(X.values, requires_grad=False)).pow(2).sum()
    print("Itertion", i, "loss:", loss.item())
    mf.zero_grad()
    loss.backward()

    with torch.no_grad():
        lr = 0.02
        for w in mf.parameters():
            w.data = (w - lr * w.grad)

            # Proximal operator: nonnegativity
            w.data = torch.clamp(w, min=0)

        mf.p.data = mask(mf.p)
        mf.p.data = mf.p / mf.p.sum(axis=0)

p = pd.DataFrame(mf.p.detach().numpy()).reindex_like(profiles)
print("Max abs error:", (p - profiles).abs().max().max())
