# RA, 2021-03-06

"""
This script downloads the human MCA dataset [2].

From [1]:

Individual layers of cortex were dissected
from tissues covering the middle temporal gyrus (MTG),
anterior cingulate gyrus (CgGr), primary visual cortex (V1C),
primary motor cortex (M1C), primary somatosensory cortex (S1C) and
primary auditory cortex (A1C) derived from human brain,
and nuclei were dissociated and sorted using the neuronal marker NeuN.
Nuclei were sampled from postmortem and neurosurgical (MTG only) donor brains,
and expression was profiled with SMART-Seq v4 or 10x v3 RNA-sequencing.

MCA: SMART-Seq v4

[1]
http://portal.brain-map.org/atlases-and-data/rnaseq/protocols-human-cortex

[2]
http://portal.brain-map.org/atlases-and-data/rnaseq/human-multiple-cortical-areas-smart-seq
"""

from pathlib import Path
from itertools import chain
from collections import Counter

import pandas as pd

from tcga.utils import download

URLS = {
    'expr': "https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/matrix.csv",
    'meta': "https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/metadata.csv",
}

out_dir = Path(__file__).with_suffix('')
download = download.to(abs_path=out_dir)

for (k, url) in URLS.items():
    (download(url).now.meta)

# with download(URLS['expr']).now.open() as fd:
#     df_expr_index = pd.read_csv(fd, sep=',', usecols=[0], index_col=0).index
#     assert (76533 == len(df_expr_index))
#
# with download(URLS['meta']).now.open() as fd:
#     df_meta_index = pd.read_csv(fd, sep=',', index_col=0).index
#     assert (df_expr_index.equals(df_meta_index[0:len(df_expr_index)]))

with download(URLS['expr']).now.open() as fd:
    df_expr = pd.read_csv(fd, sep=',', nrows=10, index_col=0).astype(int)
    assert (df_expr.shape == (len(df_expr), 50281))

with download(URLS['meta']).now.open() as fd:
    df_meta = pd.read_csv(fd, sep=',', index_col=0)
    assert (df_meta.shape == (len(df_meta), 40))


class meta_cols:
    exp_component_name = "exp_component_name"  # {..., 'SM-GE4QU_S191_E1-50': 1, 'SM-GE4QU_S192_E1-50': 1}
    cluster_label = "cluster_label"  # (!) Cell type and markers, e.g. "VLMC L1-3 CYP1B1"
    cluster_color = "cluster_color"
    cluster_order = "cluster_order"
    class_label = "class_label"  # {'Glutamatergic': 31229, 'GABAergic': 11450, 'Non-neuronal': 4753, nan: 1985}
    class_color = "class_color"
    class_order = "class_order"
    subclass_label = "subclass_label"  # {'VLMC': 11, 'IT': 21841, 'L4 IT': 3725, 'VIP': 3533, 'PVALB': 2800, 'L6 CT': 2556, 'LAMP5': 2434, 'SST': 2358, nan: 1985, 'Oligodendrocyte': 1930, 'Astrocyte': 1187, 'L6b': 1080, 'L5/6 IT Car3': 1053, 'L5/6 NP': 816, 'OPC': 773, 'Microglia': 750, 'PAX6': 325, 'L5 ET': 158, 'Endothelial': 70, 'Pericyte': 32}
    subclass_color = "subclass_color"
    subclass_order = "subclass_order"
    donor_sex_label = "donor_sex_label"  # {'M': 30899, 'F': 18518}
    donor_sex_color = "donor_sex_color"
    donor_sex_order = "donor_sex_order"
    region_label = "region_label"  # {'MTG': 16155, 'V1C': 8052, 'A1C': 6703, 'CgG': 6279, 'M1lm': 3371, 'S1ul': 3188, 'M1ul': 2864, 'S1lm': 2805}
    region_color = "region_color"
    region_order = "region_order"
    cortical_layer_label = "cortical_layer_label"  # {'L3': 9076, 'L5': 7682, 'L2': 7305, 'L6': 6620, 'L4': 5856, 'L1': 3812, 'L6a': 2443, 'L4c': 1977, 'L6b': 1830, 'L4ab': 954, 'L5b': 884, 'L5a': 668, 'WM': 310}
    cortical_layer_color = "cortical_layer_color"
    cortical_layer_order = "cortical_layer_order"
    cell_type_accession_label = "cell_type_accession_label"  # {'CS1910121062': 6502, 'CS1910121075': 2852, 'CS1910121079': 2254, 'CS1910121084': 2016, nan: 1985, 'CS1910121116': 1676, 'CS1910121055': 1556, 'CS1910121087': 1262, 'CS1910121057': 1171, 'CS1910121094': 1086, 'CS1910121040': 1015, 'CS1910121111': 966, 'CS1910121068': 893, 'CS1910121074': 885, 'CS1910121049': 880, 'CS1910121003': 859, 'CS1910121114': 773, 'CS1910121092': 768, 'CS1910121120': 750, 'CS1910121041': 734, 'CS1910121001': 704, 'CS1910121098': 682, 'CS1910121077': 626, 'CS1910121107': 584, 'CS1910121048': 570, 'CS1910121090': 495, 'CS1910121034': 443, 'CS1910121067': 433, 'CS1910121015': 428, 'CS1910121083': 428, 'CS1910121059': 428, 'CS1910121006': 423, 'CS1910121086': 410, 'CS1910121088': 390, 'CS1910121052': 379, 'CS1910121066': 362, 'CS1910121063': 360, 'CS1910121065': 357, 'CS1910121076': 349, 'CS1910121016': 329, 'CS1910121047': 325, 'CS1910121005': 315, 'CS1910121082': 310, 'CS1910121071': 281, 'CS1910121093': 273, 'CS1910121081': 268, 'CS1910121096': 259, 'CS1910121032': 258, 'CS1910121058': 256, 'CS1910121115': 254, 'CS1910121023': 251, 'CS1910121080': 240, 'CS1910121069': 211, 'CS1910121051': 202, 'CS1910121009': 197, 'CS1910121038': 191, 'CS1910121036': 189, 'CS1910121045': 177, 'CS1910121020': 171, 'CS1910121095': 170, 'CS1910121025': 162, 'CS1910121035': 148, 'CS1910121101': 147, 'CS1910121112': 146, 'CS1910121060': 138, 'CS1910121046': 137, 'CS1910121013': 130, 'CS1910121085': 130, 'CS1910121100': 129, 'CS1910121026': 120, 'CS1910121030': 113, 'CS1910121108': 102, 'CS1910121078': 101, 'CS1910121089': 99, 'CS1910121033': 99, 'CS1910121056': 98, 'CS1910121028': 96, 'CS1910121070': 94, 'CS1910121053': 93, 'CS1910121021': 89, 'CS1910121027': 89, 'CS1910121073': 88, 'CS1910121024': 87, 'CS1910121029': 85, 'CS1910121042': 83, 'CS1910121043': 82, 'CS1910121097': 82, 'CS1910121044': 82, 'CS1910121039': 81, 'CS1910121061': 78, 'CS1910121113': 75, 'CS1910121072': 71, 'CS1910121117': 70, 'CS1910121004': 69, 'CS1910121091': 69, 'CS1910121109': 68, 'CS1910121054': 66, 'CS1910121105': 66, 'CS1910121064': 58, 'CS1910121019': 57, 'CS1910121002': 54, 'CS1910121011': 52, 'CS1910121012': 49, 'CS1910121104': 44, 'CS1910121010': 43, 'CS1910121099': 40, 'CS1910121014': 37, 'CS1910121037': 35, 'CS1910121106': 35, 'CS1910121022': 34, 'CS1910121007': 33, 'CS1910121119': 32, 'CS1910121103': 28, 'CS1910121110': 27, 'CS1910121050': 26, 'CS1910121018': 25, 'CS1910121031': 25, 'CS1910121102': 20, 'CS1910121017': 19, 'CS1910121118': 11, 'CS1910121008': 10}
    cell_type_accession_color = "cell_type_accession_color"
    cell_type_accession_order = "cell_type_accession_order"
    cell_type_alias_label = "cell_type_alias_label"  # {'VLMC L1-3 CYP1B1': 11, 'Exc L2-3 LINC00507 RPL9P17': 6502, 'Exc L6 THEMIS LINC00343': 2852, 'Exc L4-5 RORB LINC01474': 2254, 'Exc L4-5 RORB RPL31P31': 2016, nan: 1985, 'Oligo L4-6 OPALIN': 1676, 'Exc L2-4 RORB GRIK1': 1556, 'Exc L4-5 RORB LCN15': 1262, 'Exc L4 RORB BHLHE22': 1171, 'Exc L6 FEZF2 FAM95C': 1086, 'Inh L3-5 SST MAFB': 1015, 'Astro L1-6 FGFR3 ETNPPL': 966, 'Exc L3-4 RORB PRSS12': 893, 'Exc L6 THEMIS EGR3': 885, 'Inh L2-4 PVALB C8orf4': 880, 'Inh L1-4 LAMP5 DUSP4': 859, 'OPC L1-6 MYT1': 773, 'Exc L6 FEZF2 VWA2': 768, 'Micro L1-6 C1QC': 750, 'Inh L4-6 SST MTHFD2P6': 734, 'Inh L1 LAMP5 NDNF': 704, 'Exc L6 FEZF2 KRT17': 682, 'Exc L3-5 RORB CMAHP': 626, 'Exc L5-6 FEZF2 MYBPHL': 584, 'Inh L5 PVALB CNTNAP3P2': 570, 'Exc L5-6 THEMIS GPR21': 495, 'Inh L1-3 VIP ZNF322P1': 443, 'Exc L3-4 RORB SEMA6D': 433, 'Inh L1 ADARB2 ADAM33': 428, 'Exc L4-5 RORB HNRNPA1P46': 428, 'Exc L4 RORB CACNG5': 428, 'Inh L5-6 LAMP5 SFTA3': 423, 'Exc L5 RORB SNHG7': 410, 'Exc L5-6 THEMIS IL7R': 390, 'Inh L3-4 PVALB HOMER3': 379, 'Exc L3-4 RORB FOLH1B': 362, 'Exc L3 LINC00507 PSRC1': 360, 'Exc L3 RORB CARTPT': 357, 'Exc L5-6 THEMIS TMEM233': 349, 'Inh L1 SST CXCL14': 329, 'Inh L5-6 PVALB STON2': 325, 'Inh L1-6 LAMP5 CA13': 315, 'Exc L5 RORB LINC01202': 310, 'Exc L3 LINC00507 CTXN3': 281, 'Exc L5-6 FEZF2 ANKRD20A1': 273, 'Exc L5-6 RORB LINC00320': 268, 'Exc L6 FEZF2 ETV4': 259, 'Inh L2-5 VIP TOX2': 258, 'Exc L4 RORB CCDC168': 256, 'Oligo L4-6 MOBP COL18A1': 254, 'Inh L2-6 VIP VIP': 251, 'Exc L4-6 RORB HPCA': 240, 'Exc L3-5 RORB HSPB3': 211, 'Inh L1-3 PVALB WFDC2': 202, 'Inh L1 PAX6 CA4': 197, 'Inh L5-6 SST ISOC1': 191, 'Inh L1-3 VIP GGH': 189, 'Inh L5-6 PVALB FAM150B': 177, 'Inh L1-5 VIP KCNJ2': 171, 'Exc L6 FEZF2 CPZ': 170, 'Inh L1-6 VIP RGS16': 162, 'Inh L3 VIP CBLN1': 148, 'Exc L6 FEZF2 TBCC': 147, 'Astro L1 FGFR3 MT1G': 146, 'Exc L4-5 RORB ASCL1': 138, 'Inh L2-4 SST AHR': 137, 'Inh L1 VIP PRSS8': 130, 'Exc L3-5 THEMIS UBE2F': 130, 'Exc L6 FEZF2 SLITRK6': 129, 'Inh L1-3 VIP SSTR1': 120, 'Inh L1-4 VIP CHRNA2': 113, 'Exc L5-6 FEZF2 CYP26B1': 102, 'Exc L3-5 RORB CD24': 101, 'Exc L6 THEMIS C6orf48': 99, 'Inh L2-4 VIP DSEL': 99, 'Exc L3-4 RORB RPS3P6': 98, 'Inh L1-3 VIP ACHE': 96, 'Exc L3-5 THEMIS ELOF1': 94, 'Inh L1-6 PVALB SCUBE3': 93, 'Inh L1 VIP PCDH20': 89, 'Inh L1-2 VIP RPL41P3': 89, 'Exc L5-6 THEMIS OR1J1': 88, 'Inh L3-6 VIP KCTD13': 87, 'Inh L2-4 VIP LGI2': 85, 'Inh L4-5 PVALB TRIM67': 83, 'Inh L5-6 SST TH': 82, 'Exc L6 FEZF2 TBC1D26': 82, 'Inh L6 LHX6 GLP1R': 82, 'Inh L5-6 SST KLHL14': 81, 'Exc L4-5 RORB AIM2': 78, 'Astro L1 FGFR3 FOS': 75, 'Exc L3 THEMIS PLA2G7': 71, 'Endo L2-5 CLDN5': 70, 'Inh L6 LAMP5 C1QL2': 69, 'Exc L5-6 THEMIS THTPA': 69, 'Exc L5-6 FEZF2 RSAD2': 68, 'Inh L3-6 PVALB MFI2': 66, 'Exc L5 FEZF2 MORN2': 66, 'Exc L3-5 LINC00507 SLN': 58, 'Inh L1-6 VIP PENK': 57, 'Inh L1 LAMP5 GGT8P': 54, 'Inh L1-3 PAX6 NABP1': 52, 'Inh L1-6 VIP RCN1': 49, 'Exc L5 FEZF2 SCN7A': 44, 'Inh L1 PAX6 GRIP2': 43, 'Exc L6 FEZF2 P4HA3': 40, 'Inh L1 VIP TNFAIP8L3': 37, 'Inh L6 SST NPY': 35, 'Exc L5 FEZF2 DYRK2': 35, 'Inh L1-2 VIP PPAPDC1A': 34, 'Inh L1-2 PAX6 SCGN': 33, 'Peri L1-6 MUSTN1': 32, 'Exc L3-5 FEZF2 DCN': 28, 'Exc L5-6 FEZF2 CABP7': 27, 'Inh L1-2 PVALB TAC1': 26, 'Inh L1 VIP SOX11': 25, 'Inh L1-3 VIP CCDC184': 25, 'Exc L3-5 FEZF2 ONECUT1': 20, 'Inh L1 ADARB2 DISP2': 19, 'Inh L6 LAMP5 ANKRD20A11P': 10}
    cell_type_alias_color = "cell_type_alias_color"
    # cell_type_alias_order = "cell_type_alias_order"
    cell_type_order = "cell_type_order"  # {62.0: 6502, 75.0: 2852, 79.0: 2254, ...}
    cell_type_alt_alias_label = "cell_type_alt_alias_label"  # {nan: 47626, 'Lamp5 Rosehip': 859, 'Lamp5 Lhx6 2': 423, 'Lamp5 Lhx6 1': 315, 'Chandelier 1': 93, 'Chandelier 2': 66, 'Sst Chodl': 35}
    cell_type_alt_alias_color = "cell_type_alt_alias_color"
    cell_type_alt_alias_order = "cell_type_alt_alias_order"
    cell_type_designation_label = "cell_type_designation_label"  # {'Neuron 062': 6502, 'Neuron 075': 2852, 'Neuron 079': 2254, 'Neuron 084': 2016, nan: 1985, 'Non-neuron 006': 1676, 'Neuron 055': 1556, 'Neuron 087': 1262, 'Neuron 057': 1171, 'Neuron 094': 1086, 'Neuron 040': 1015, 'Non-neuron 001': 966, 'Neuron 068': 893, 'Neuron 074': 885, 'Neuron 049': 880, 'Neuron 003': 859, 'Non-neuron 004': 773, 'Neuron 092': 768, 'Non-neuron 010': 750, 'Neuron 041': 734, 'Neuron 001': 704, 'Neuron 098': 682, 'Neuron 077': 626, 'Neuron 107': 584, 'Neuron 048': 570, 'Neuron 090': 495, 'Neuron 034': 443, 'Neuron 067': 433, 'Neuron 015': 428, 'Neuron 083': 428, 'Neuron 059': 428, 'Neuron 006': 423, 'Neuron 086': 410, 'Neuron 088': 390, 'Neuron 052': 379, 'Neuron 066': 362, 'Neuron 063': 360, 'Neuron 065': 357, 'Neuron 076': 349, 'Neuron 016': 329, 'Neuron 047': 325, 'Neuron 005': 315, 'Neuron 082': 310, 'Neuron 071': 281, 'Neuron 093': 273, 'Neuron 081': 268, 'Neuron 096': 259, 'Neuron 032': 258, 'Neuron 058': 256, 'Non-neuron 005': 254, 'Neuron 023': 251, 'Neuron 080': 240, 'Neuron 069': 211, 'Neuron 051': 202, 'Neuron 009': 197, 'Neuron 038': 191, 'Neuron 036': 189, 'Neuron 045': 177, 'Neuron 020': 171, 'Neuron 095': 170, 'Neuron 025': 162, 'Neuron 035': 148, 'Neuron 101': 147, 'Non-neuron 002': 146, 'Neuron 060': 138, 'Neuron 046': 137, 'Neuron 013': 130, 'Neuron 085': 130, 'Neuron 100': 129, 'Neuron 026': 120, 'Neuron 030': 113, 'Neuron 108': 102, 'Neuron 078': 101, 'Neuron 089': 99, 'Neuron 033': 99, 'Neuron 056': 98, 'Neuron 028': 96, 'Neuron 070': 94, 'Neuron 053': 93, 'Neuron 021': 89, 'Neuron 027': 89, 'Neuron 073': 88, 'Neuron 024': 87, 'Neuron 029': 85, 'Neuron 042': 83, 'Neuron 043': 82, 'Neuron 097': 82, 'Neuron 044': 82, 'Neuron 039': 81, 'Neuron 061': 78, 'Non-neuron 003': 75, 'Neuron 072': 71, 'Non-neuron 007': 70, 'Neuron 004': 69, 'Neuron 091': 69, 'Neuron 109': 68, 'Neuron 054': 66, 'Neuron 105': 66, 'Neuron 064': 58, 'Neuron 019': 57, 'Neuron 002': 54, 'Neuron 011': 52, 'Neuron 012': 49, 'Neuron 104': 44, 'Neuron 010': 43, 'Neuron 099': 40, 'Neuron 014': 37, 'Neuron 037': 35, 'Neuron 106': 35, 'Neuron 022': 34, 'Neuron 007': 33, 'Non-neuron 009': 32, 'Neuron 103': 28, 'Neuron 110': 27, 'Neuron 050': 26, 'Neuron 018': 25, 'Neuron 031': 25, 'Neuron 102': 20, 'Neuron 017': 19, 'Non-neuron 008': 11, 'Neuron 008': 10}
    cell_type_designation_color = "cell_type_designation_color"
    cell_type_designation_order = "cell_type_designation_order"
    external_donor_name_label = "external_donor_name_label"  # {'H200.1030': 19575, 'H200.1023': 18518, 'H200.1025': 11324}
    external_donor_name_color = "external_donor_name_color"
    external_donor_name_order = "external_donor_name_order"
    specimen_type = "specimen_type"  # {'nucleus': 49417}
    full_genotype_label = "full_genotype_label"  # ?
    full_genotype_color = "full_genotype_color"
    full_genotype_order = "full_genotype_order"
    outlier_call = "outlier_call"  # {False: 47432, True: 1985}
    outlier_type = "outlier_type"  # {nan: 47432, 'Outlier L1-3 SST OR2AD1P': 516, 'Outlier L2-3 SST CALB1': 347, 'Outlier L5-6 SST KLHDC8A': 198, 'Outlier L6 FEZF2 RPS17P9': 157, 'Outlier L1-6 RORB PTP4A2P1': 150, 'Outlier L4-6 SLC17A7 CALM2P3': 111, 'Outlier L1-6 THEMIS CARM1P1': 74, 'Donor L6 SLC17A7 RARRES2': 68, 'Outlier L3 SLC17A7 RIN1': 68, 'Outlier L4-6 RORB TSHZ1': 59, 'Outlier L1-4 SST DNAJC19P8': 43, 'Outlier L1-2 THEMIS IGH': 42, 'Outlier L6 SLC17A7 GAL3ST4': 38, 'Outlier L1-6 THEMIS PARD3B': 31, 'Outlier L5-6 FEZF2 NFATC2IP': 19, 'Outlier L1-3 FGFR3 RLBP1': 16, 'Outlier L6 PDGFRA NKAIN4': 16, 'Outlier L5-6 THEMIS RUNDC1': 12, 'Outlier L2-6 GAD1 BOD1': 9, 'Donor L4 RORB LINC00702': 6, 'Donor L2-4 GAD1 GPR17': 5}

    assert (set(df_meta.columns) == set(v for v in locals() if not v.startswith("__")))
    assert df_meta.cluster_label.equals(df_meta.cell_type_alias_label)


celltypes = df_meta.cluster_label.str.split(" ", n=2, expand=True).rename(
    columns={0: 'type', 1: 'layer', 2: 'markers'}
)

# sample_name           type layer           markers
# F2S4_160113_027_A01    NaN   NaN               NaN
# F2S4_160113_027_B01    Inh  L2-5          VIP TOX2
# F2S4_160113_027_C01    Inh    L1       LAMP5 GGT8P
# F2S4_160113_027_D01    Inh    L1        LAMP5 NDNF
# F2S4_160113_027_E01    Inh  L1-3      VIP ZNF322P1
# ...                    ...   ...               ...
# F2S4_190227_100_C01  Astro  L1-6      FGFR3 ETNPPL
# F2S4_190227_100_E01    Exc    L6  THEMIS LINC00343
# F2S4_190227_100_F01    Inh    L1        LAMP5 NDNF
# F2S4_190227_100_G01  Oligo  L4-6            OPALIN
# F2S4_190227_100_H01  Oligo  L4-6            OPALIN


markers = Counter(chain.from_iterable(celltypes.markers.dropna().str.split(' ')))
# {'RORB': 13896, 'LINC00507': 7201, 'RPL9P17': 6502, 'THEMIS': 5522, 'FEZF2': 4610, 'VIP': 3008, 'LINC00343': 2852, 'PVALB': 2801, 'SST': 2604, 'LAMP5': 2434, 'LINC01474': 2254, 'RPL31P31': 2016, 'OPALIN': 1676, 'GRIK1': 1556, 'LCN15': 1262, 'FGFR3': 1187, 'BHLHE22': 1171, 'FAM95C': 1086, 'MAFB': 1015, 'ETNPPL': 966, 'PRSS12': 893, 'EGR3': 885, 'C8orf4': 880, 'DUSP4': 859, 'MYT1': 773, 'VWA2': 768, 'C1QC': 750, 'MTHFD2P6': 734, 'NDNF': 704, 'KRT17': 682, 'CMAHP': 626, 'MYBPHL': 584, 'CNTNAP3P2': 570, 'GPR21': 495, 'ADARB2': 447, 'ZNF322P1': 443, 'SEMA6D': 433, 'ADAM33': 428, 'HNRNPA1P46': 428, 'CACNG5': 428, 'SFTA3': 423, 'SNHG7': 410, 'IL7R': 390, 'HOMER3': 379, 'FOLH1B': 362, 'PSRC1': 360, 'CARTPT': 357, 'TMEM233': 349, 'CXCL14': 329, 'PAX6': 325, 'STON2': 325, 'CA13': 315, 'LINC01202': 310, 'CTXN3': 281, 'ANKRD20A1': 273, 'LINC00320': 268, 'ETV4': 259, 'TOX2': 258, 'CCDC168': 256, 'MOBP': 254, 'COL18A1': 254, 'HPCA': 240, 'HSPB3': 211, 'WFDC2': 202, 'CA4': 197, 'ISOC1': 191, 'GGH': 189, 'FAM150B': 177, 'KCNJ2': 171, 'CPZ': 170, 'RGS16': 162, 'CBLN1': 148, 'TBCC': 147, 'MT1G': 146, 'ASCL1': 138, 'AHR': 137, 'PRSS8': 130, 'UBE2F': 130, 'SLITRK6': 129, 'SSTR1': 120, 'CHRNA2': 113, 'CYP26B1': 102, 'CD24': 101, 'C6orf48': 99, 'DSEL': 99, 'RPS3P6': 98, 'ACHE': 96, 'ELOF1': 94, 'SCUBE3': 93, 'PCDH20': 89, 'RPL41P3': 89, 'OR1J1': 88, 'KCTD13': 87, 'LGI2': 85, 'TRIM67': 83, 'TH': 82, 'TBC1D26': 82, 'LHX6': 82, 'GLP1R': 82, 'KLHL14': 81, 'AIM2': 78, 'FOS': 75, 'PLA2G7': 71, 'CLDN5': 70, 'C1QL2': 69, 'THTPA': 69, 'RSAD2': 68, 'MFI2': 66, 'MORN2': 66, 'SLN': 58, 'PENK': 57, 'GGT8P': 54, 'NABP1': 52, 'RCN1': 49, 'SCN7A': 44, 'GRIP2': 43, 'P4HA3': 40, 'TNFAIP8L3': 37, 'NPY': 35, 'DYRK2': 35, 'PPAPDC1A': 34, 'SCGN': 33, 'MUSTN1': 32, 'DCN': 28, 'CABP7': 27, 'TAC1': 26, 'SOX11': 25, 'CCDC184': 25, 'ONECUT1': 20, 'DISP2': 19, 'CYP1B1': 11, 'ANKRD20A11P': 10}

(df_meta.cluster_label[df_meta.class_label == 'Non-neuronal'].value_counts())
# Oligo L4-6 OPALIN          1676
# Astro L1-6 FGFR3 ETNPPL     966
# OPC L1-6 MYT1               773
# Micro L1-6 C1QC             750
# Oligo L4-6 MOBP COL18A1     254
# Astro L1 FGFR3 MT1G         146
# Astro L1 FGFR3 FOS           75
# Endo L2-5 CLDN5              70
# Peri L1-6 MUSTN1             32
# VLMC L1-3 CYP1B1             11

print(df_meta.cluster_label[df_meta.class_label == "Glutamatergic"].value_counts())
assert all(df_meta.cluster_label[df_meta.class_label == "Glutamatergic"].str.slice(stop=3) == "Exc")
# Exc L2-3 LINC00507 RPL9P17    6502
# Exc L6 THEMIS LINC00343       2852
# Exc L4-5 RORB LINC01474       2254
# Exc L4-5 RORB RPL31P31        2016
# ...
# Exc L5 FEZF2 DYRK2              35
# Exc L3-5 FEZF2 DCN              28
# Exc L5-6 FEZF2 CABP7            27
# Exc L3-5 FEZF2 ONECUT1          20

print(df_meta.cluster_label[df_meta.class_label == "GABAergic"].value_counts())
assert all(df_meta.cluster_label[df_meta.class_label == "GABAergic"].str.slice(stop=3) == "Inh")
# Inh L3-5 SST MAFB           1015
# Inh L2-4 PVALB C8orf4        880
# Inh L1-4 LAMP5 DUSP4         859
# Inh L4-6 SST MTHFD2P6        734
# Inh L1 LAMP5 NDNF            704
# Inh L5 PVALB CNTNAP3P2       570
# Inh L1-3 VIP ZNF322P1        443
# ...
# Inh L1-3 VIP CCDC184          25
# Inh L1 VIP SOX11              25
# Inh L1 ADARB2 DISP2           19
# Inh L6 LAMP5 ANKRD20A11P      10
