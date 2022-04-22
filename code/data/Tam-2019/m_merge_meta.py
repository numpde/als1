
from pathlib import Path

import pandas as pd

from tcga.utils import unlist1, mkdir

df_meta = pd.read_table(unlist1(Path(__file__).parent.glob("**/*meta.txt.gz")))

df_manifest = pd.read_table(unlist1(Path(__file__).parent.glob("**/*manifest.tsv")))

df_subtype_nygc = pd.read_table(unlist1(Path(__file__).parent.glob("**/S1B_NYGC_subtypes.csv")), comment="#")
df_subtype_uscd = pd.read_table(unlist1(Path(__file__).parent.glob("**/S1C_UCSD.csv")), comment="#")

assert not df_meta['subject_id'].is_unique
assert (95 == len(set(df_meta['subject_id'])))
assert df_subtype_nygc['patient_id'].is_unique
assert df_subtype_uscd['patient_id'].is_unique

# print(df_meta.to_markdown())
# print(df_manifest.to_markdown())

df_meta['experiment_accession'] = df_meta['sra'].str[-10:]
assert df_meta['experiment_accession'].is_unique

df_merged = df_meta.merge(df_manifest, how='inner', on='experiment_accession')

assert df_merged['run_accession'].is_unique  # e.g., SRR8375274
assert df_merged['sample_accession'].is_unique  # e.g., SAMN10656686
assert df_merged['experiment_accession'].is_unique  # e.g., SRX5185320

# print(df_subtype_nygc.to_markdown())
# print(df_subtype_uscd.to_markdown())

common_fields = ['patient_id', 'ALS_subtype', 'age_of_death', 'gender', 'site_of_onset']

df_subtype = pd.concat(
    objs=[df_subtype_nygc[common_fields], df_subtype_uscd[common_fields]],
    axis=0,
)

df_merged = df_merged.merge(df_subtype, how='left', left_on='subject_id', right_on='patient_id')

df_merged.to_csv(mkdir(Path(__file__).with_suffix('')) / "meta_merged.tsv", sep='\t', index=False)
