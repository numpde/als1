# RA, 2021-05-02

"""
Prepare the scRNA-seq dataset from:

CNS fibroblasts form a fibrotic scar in response to immune cell infiltration.
Dorrier CE, Aran D, Haenelt EA et al.
Nat Neurosci 24, 234--244 (2021).
https://doi.org/10.1038/s41593-020-00770-9

GSE:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135186
"""

from bugs import *

from tcga.utils import download
import tcga.utils
import pandas
import pathlib

src_dir = next(p for p in pathlib.Path.cwd().parents for p in p.glob("**/Daneman-2020/z_*")).with_suffix('')

download = download.to(abs_path=(src_dir / "download_cache"))

URLS = {
    # 'all': "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE135186&format=file",
    'healthy1': "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3993nnn/GSM3993136/suppl/GSM3993136_Healthy1.txt.gz",
    'healthy2': "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3993nnn/GSM3993137/suppl/GSM3993137_Healthy2.txt.gz",
    'healthy3': "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4751nnn/GSM4751220/suppl/GSM4751220_health3.txt.gz",
    'EAE1': "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3993nnn/GSM3993138/suppl/GSM3993138_EAE.txt.gz",
    'EAE2': "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4751nnn/GSM4751221/suppl/GSM4751221_EAE2.txt.gz",
}

if __name__ == '__main__':
    print("Download to:", tcga.utils.relpath(download.local_folder))

    for (_, url) in URLS.items():
        print("Downloading:", url)
        print("OK:", download(url).now.meta)


def make_df_expr_gse():
    def url2df(k):
        with download(URLS[k]).now.open(mode='rb') as fd:
            df = pandas.read_table(fd, compression='gzip', index_col=0).astype(int).sort_index(axis=1)
            df.columns = [f"{k}_{c}" for c in df.columns]
            return df

    df = pandas.concat(axis=1, verify_integrity=True, objs=(
        url2df(k)
        for k in URLS
    ))

    df = df.dropna(axis=0).astype(int)

    # Expected number of samples
    assert (7205 == len(df.columns))

    assert df.index.is_unique

    # Collapse duplicate genes
    # df = df.groupby(df.index).median().round().astype(int)

    # Samples x Genes
    df = df.T

    return df


def cookup_meta(df_expr):
    s = pd.Series(index=df_expr.index, data=df_expr.index).str.split('_').apply(first)
    df_meta = pd.DataFrame(index=df_expr.index).assign(group=s.apply(lambda x: x[:-1])).assign(subgroup=s)
    return df_meta


if __name__ == '__main__':
    df_expr = make_df_expr_gse()
    df_meta = cookup_meta(df_expr)
    assert df_meta.index.equals(df_expr.index)

    from c_clusters import assignment
    df_meta = df_meta.assign(cluster=assignment.x)

    df_meta.to_csv(src_dir / "meta.csv.gz", compression='gzip', sep='\t')
    df_expr.to_csv(src_dir / "expr.csv.gz", compression='gzip', sep='\t')

df_meta = pandas.read_table(src_dir / "meta.csv.gz", sep='\t', index_col=0)
df_expr = pandas.read_table(src_dir / "expr.csv.gz", sep='\t', index_col=0)

assert df_meta.index.equals(df_expr.index)
