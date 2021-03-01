# RA, 2021-03-21

from bugs import *
from tcga.utils import First
from plox import Plox, rcParam
from numpy.random import default_rng

out_dir = mkdir(Path(__file__).with_suffix(''))

read = First(Path("f_fewer_cells").glob).then(unlist1).then(lambda f: pd.read_table(f, index_col=0))

# # Fake expression table
# rng = default_rng(seed=0)
# (ngenes, nsamples) = (2000, 1000)
# phi = 0.5  # dispersion
# df_expr = pd.DataFrame(
#     index=pd.Index(range(ngenes), name="gene"),
#     columns=pd.Index(range(nsamples), name="sample"),
#     data=[
#         rng.negative_binomial(n=phi, p=(1 / (1 + mu / phi)), size=nsamples)
#         for mu in rng.gamma(shape=0.1, scale=0.1, size=ngenes)
#     ],
#     dtype=int,
# )
# del phi

df_expr = read("data*")
# df_meta = read("meta*")

# Expression table should be genes x samples

df_genewise = pd.DataFrame({'M': df_expr.mean(axis=1), 'V': df_expr.var(axis=1)}).sort_values(by='M')
df_smplwise = pd.DataFrame({'S': df_expr.sum(axis=0), 'V': df_expr.var(axis=0)}).sort_values(by='S')

# print(f"For mu_g = mean expression of g across all samples:")
# print(f"    Mean of mu: ", df_genewise.M.mean())
# print(f"    Var  of mu: ", df_genewise.M.var())

from matplotlib.style import context as PlotStyle

common_style = {
    rcParam.Text.usetex: True,
    rcParam.Font.size: 14,
}

with PlotStyle(common_style):
    # Library size
    with Plox() as px:
        # [hh, ee] = np.histogram(np.log(df_smplwise.S[df_smplwise.S > 0]), bins=20)
        # ee = np.exp(ee)
        [hh, ee] = np.histogram(df_smplwise.S[df_smplwise.S > 0], bins=20)
        px.a.bar(ee[0:-1], hh, width=np.diff(ee), edgecolor='w', align='edge')
        px.a.set_yscale('log')
        # px.a.set_xscale('log')
        px.a.set_ylabel("Number of samples")
        px.a.set_xlabel("Library size (sum of raw counts)")
        # np.floor(np.log10(min(px.a.get_xlim())), np.ceil(np.log10(min(px.a.get_xlim()))
        px.a.plot(0, 0, '.', c='w', label=f"\# total: {len(df_smplwise.S)}")
        px.a.plot(0, 0, '.', c='w', label=f"\# zeros: {sum(df_smplwise.S.dropna() == 0)}")
        px.a.plot(0, 0, '.', c='w', label=f"\# n/a: {sum(df_smplwise.S.isna())}")
        px.a.legend(frameon=False, loc='upper right')
        px.f.savefig(out_dir / f"smplwise_libsize-frequency.png")

    # Distribution of genewise mean expression
    with Plox() as px:
        [hh, ee] = np.histogram(np.log(df_genewise.M[df_genewise.M > 0]), bins=20)
        ee = np.exp(ee)
        px.a.bar(ee[0:-1], hh, width=np.diff(ee), edgecolor='w', align='edge')
        px.a.set_yscale('log')
        px.a.set_xscale('log')
        px.a.set_ylabel(rf"Out of {sum(df_genewise.M.dropna() > 0)} nonzero genes, \textbf{{how many}} ...")
        px.a.set_xlabel(rf"... have this \textbf{{mean expression}} across {len(df_expr.columns)} samples")
        px.a.set_xlim(*px.a.get_xlim())
        px.a.set_ylim(*px.a.get_ylim())
        px.a.plot(0, 0, '.', c='w', label=f"\# total: {len(df_genewise.M)}")
        px.a.plot(0, 0, '.', c='w', label=f"\# zeros: {sum(df_genewise.M.dropna() == 0)}")
        px.a.plot(0, 0, '.', c='w', label=f"\# n/a: {sum(df_genewise.M.isna())}")
        px.a.legend(frameon=False, loc='upper right')
        px.f.savefig(out_dir / f"genewise_mean-frequency.png")

    # Mean-variance-zeros relations

    import statsmodels.api as sm

    V = df_genewise.V[~df_genewise.M.isna()]
    M = df_genewise.M[~df_genewise.M.isna()]
    # Remove "outliers" for nicer fit
    # ii = ((V - M) < pd.Series(V - M).quantile(q=1))
    phi_hat = sm.OLS(M ** 2, (V - M)).fit().params.x1
    # print("Estimated dispersion:", phi_hat)

    plot_kw = dict(marker='.', ls='none', alpha=0.3, ms=5, mfc='b', mec='none')

    # Variance vs Mean
    with Plox() as px:
        M = df_genewise.M[df_genewise.M > 0]
        V = df_genewise.V[df_genewise.M > 0]
        px.a.plot(M, V, **plot_kw)
        px.a.plot(0, 0, **{**plot_kw, **dict(ms=10)}, label="Nonzero genes")
        px.a.set_xscale('log')
        px.a.set_yscale('log')
        px.a.set_ylabel(rf"Genewise \textbf{{variance}} and ...")
        px.a.set_xlabel(rf"... \textbf{{mean}} of expression across {len(df_expr.columns)} samples")
        px.a.set_xlim(*px.a.get_xlim())
        px.a.set_ylim(*px.a.get_ylim())
        mm = np.logspace(np.log10(min(px.a.get_xlim())), np.log10(max(px.a.get_xlim())), 100)
        vv = mm + (mm ** 2) / phi_hat
        px.a.plot(mm, vv, '-', alpha=0.7, color='r', lw=2, label=rf"NB with $\phi = {phi_hat:.03g}$")
        px.a.legend(frameon=True, loc='upper left')
        px.f.savefig(out_dir / f"genewise_mean-variance.png")

    # Fraction-of-zeros vs Mean
    with Plox() as px:
        M = df_genewise.M[df_genewise.M > 0].sort_values()
        F = (1 / (1 + M / phi_hat)) ** phi_hat
        s = pd.Series(index=df_expr.mean(axis=1), data=(df_expr == 0).mean(axis=1).values).sort_index()
        px.a.plot(s.index, s.values, **plot_kw)
        h = px.a.plot(1, 1, **{**plot_kw, **dict(ms=10)}, label="Nonzero genes").pop()
        px.a.plot(M, F, '-', alpha=0.7, color='r', lw=2, label=rf"NB with $\phi = {phi_hat:.03g}$")
        px.a.set_xscale('log')
        px.a.set_ylabel(rf"Genewise \textbf{{fraction of zeros}} and ...")
        px.a.set_xlabel(rf"... \textbf{{mean}} of expression across {len(df_expr.columns)} samples")
        px.a.legend(frameon=True, loc='lower left')
        h.remove()
        px.f.savefig(out_dir / f"genewise_mean-zerofrac.png")
