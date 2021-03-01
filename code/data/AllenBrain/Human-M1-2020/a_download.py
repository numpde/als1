# RA, 2020-12-10

"""
This script only downloads the human M1 dataset [2].

From [1]:

Individual layers of cortex were dissected
from tissues covering the middle temporal gyrus (MTG),
anterior cingulate gyrus (CgGr), primary visual cortex (V1C),
primary motor cortex (M1C), primary somatosensory cortex (S1C) and
primary auditory cortex (A1C) derived from human brain,
and nuclei were dissociated and sorted using the neuronal marker NeuN.
Nuclei were sampled from postmortem and neurosurgical (MTG only) donor brains,
and expression was profiled with SMART-Seq v4 or 10x v3 RNA-sequencing.

M1: 10x v3 RNA-sequencing

[1]
http://portal.brain-map.org/atlases-and-data/rnaseq/protocols-human-cortex

[2]
https://portal.brain-map.org/atlases-and-data/rnaseq/human-m1-10x
"""

from pathlib import Path
from itertools import chain
from collections import Counter

import pandas as pd

from tcga.utils import download

URLS = {
    'expr': "https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/matrix.csv",
    'meta': "https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/metadata.csv",
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
    del df_expr

with download(URLS['meta']).now.open() as fd:
    df_meta = pd.read_csv(fd, sep=',', index_col=0)
    assert (df_meta.shape == (len(df_meta), 38))


class meta_cols:
    exp_component_name = "exp_component_name"
    cluster_label = "cluster_label"  # (!) Cell type and markers, e.g. "Exc L2 LINC00507 GLRA3"
    cluster_color = "cluster_color"
    cluster_order = "cluster_order"
    class_label = "class_label"  # {'Glutamatergic': 48536, 'GABAergic': 23992, 'Non-Neuronal': 4005}
    class_color = "class_color"
    class_order = "class_order"
    subclass_label = "subclass_label"  # {'VLMC': 40, 'L2/3 IT': 24231, 'L5 IT': 13834, 'Pvalb': 7782, 'Sst': 5936, 'Vip': 4858, 'Lamp5': 4454, 'L6 CT': 3734, 'Oligo': 2942, 'L6b': 2236, 'L6 IT': 1829, 'L5/6 NP': 1487, 'Sncg': 895, 'L5 ET': 858, 'Astro': 568, 'L6 IT Car3': 327, 'OPC': 283, 'Micro-PVM': 108, 'Sst Chodl': 67, 'Endo': 64}
    subclass_color = "subclass_color"
    subclass_order = "subclass_order"
    donor_sex_label = "donor_sex_label"  # {'F': 42728, 'M': 33805}
    donor_sex_color = "donor_sex_color"
    donor_sex_order = "donor_sex_order"
    region_label = "region_label"  # {'M1': 76533}
    region_color = "region_color"
    region_order = "region_order"
    cortical_layer_label = "cortical_layer_label"  # {'all': 76533}
    cortical_layer_color = "cortical_layer_color"
    cortical_layer_order = "cortical_layer_order"
    cell_type_accession_label = "cell_type_accession_label"  # {'CS1912131075': 12485, 'CS1912131090': 4631, 'CS1912131078': 4044, 'CS1912131080': 3331, 'CS1912131064': 3253, 'CS1912131101': 2602, 'CS1912131120': 2300, 'CS1912131086': 2134, 'CS1912131088': 1824, 'CS1912131108': 1801, 'CS1912131081': 1628, 'CS1912131003': 1539, 'CS1912131065': 1326, 'CS1912131091': 1231, 'CS1912131076': 1211, 'CS1912131046': 1154, 'CS1912131089': 1132, 'CS1912131077': 1088, 'CS1912131051': 1069, 'CS1912131006': 923, 'CS1912131094': 905, 'CS1912131082': 785, 'CS1912131042': 764, 'CS1912131116': 753, 'CS1912131087': 640, 'CS1912131029': 638, 'CS1912131062': 631, 'CS1912131111': 620, 'CS1912131060': 563, 'CS1912131083': 535, 'CS1912131119': 535, 'CS1912131020': 509, 'CS1912131096': 488, 'CS1912131005': 465, 'CS1912131019': 464, 'CS1912131063': 456, 'CS1912131095': 436, 'CS1912131047': 419, 'CS1912131004': 419, 'CS1912131021': 396, 'CS1912131043': 394, 'CS1912131028': 390, 'CS1912131124': 377, 'CS1912131052': 376, 'CS1912131039': 367, 'CS1912131014': 352, 'CS1912131097': 327, 'CS1912131050': 320, 'CS1912131102': 318, 'CS1912131025': 310, 'CS1912131092': 309, 'CS1912131026': 288, 'CS1912131072': 288, 'CS1912131118': 283, 'CS1912131113': 280, 'CS1912131098': 263, 'CS1912131073': 261, 'CS1912131061': 254, 'CS1912131084': 250, 'CS1912131002': 246, 'CS1912131045': 242, 'CS1912131100': 241, 'CS1912131071': 230, 'CS1912131114': 229, 'CS1912131030': 225, 'CS1912131018': 223, 'CS1912131044': 218, 'CS1912131085': 216, 'CS1912131048': 194, 'CS1912131068': 194, 'CS1912131015': 193, 'CS1912131057': 188, 'CS1912131016': 186, 'CS1912131009': 185, 'CS1912131074': 183, 'CS1912131103': 159, 'CS1912131099': 151, 'CS1912131110': 150, 'CS1912131093': 147, 'CS1912131035': 146, 'CS1912131008': 139, 'CS1912131069': 138, 'CS1912131038': 131, 'CS1912131049': 129, 'CS1912131023': 128, 'CS1912131037': 127, 'CS1912131109': 127, 'CS1912131066': 123, 'CS1912131122': 119, 'CS1912131117': 117, 'CS1912131034': 116, 'CS1912131024': 111, 'CS1912131070': 111, 'CS1912131107': 108, 'CS1912131011': 108, 'CS1912131127': 108, 'CS1912131059': 108, 'CS1912131115': 108, 'CS1912131067': 107, 'CS1912131105': 106, 'CS1912131079': 101, 'CS1912131056': 100, 'CS1912131027': 95, 'CS1912131001': 94, 'CS1912131055': 93, 'CS1912131036': 92, 'CS1912131031': 89, 'CS1912131112': 88, 'CS1912131033': 85, 'CS1912131022': 84, 'CS1912131058': 84, 'CS1912131007': 84, 'CS1912131041': 81, 'CS1912131053': 76, 'CS1912131106': 72, 'CS1912131123': 72, 'CS1912131032': 67, 'CS1912131040': 67, 'CS1912131125': 64, 'CS1912131010': 61, 'CS1912131013': 52, 'CS1912131012': 42, 'CS1912131126': 40, 'CS1912131017': 38, 'CS1912131054': 35, 'CS1912131104': 22, 'CS1912131121': 6}
    cell_type_accession_color = "cell_type_accession_color"
    cell_type_accession_order = "cell_type_accession_order"
    cell_type_alias_label = "cell_type_alias_label"  # {'Exc L2 LINC00507 GLRA3': 12485, 'Exc L3-5 RORB LNX2': 4631, 'Exc L3 LAMP5 CARM1P1': 4044, 'Exc L2-3 RORB CCDC68': 3331, 'Inh L2-5 PVALB RPH3AL': 3253, 'Exc L5-6 FEZF2 C9orf135-AS1': 2602, 'Oligo L2-6 OPALIN FTH1P3': 2300, 'Exc L3 RORB OTOGL': 2134, 'Exc L5 THEMIS SLC22A18': 1824, 'Exc L6 FEZF2 KLK7': 1801, 'Exc L2-3 RORB PTPN3': 1628, 'Inh L1-6 LAMP5 AARD': 1539, 'Inh L3 PVALB SAMD13': 1326, 'Exc L3-5 RORB RPRM': 1231, 'Exc L2-3 RORB RTKN2': 1211, 'Inh L3-5 SST GGTLC3': 1154, 'Exc L5-6 THEMIS TNFAIP6': 1132, 'Exc L2-3 LINC00507 DSG3': 1088, 'Inh L1-3 SST FAM20A': 1069, 'Inh L5-6 LAMP5 CRABP1': 923, 'Exc L6 THEMIS LINC00343': 905, 'Exc L3 THEMIS ENPEP': 785, 'Inh L5 SST RPL35AP11': 764, 'Exc L5-6 FEZF2 IFNG-AS1': 753, 'Exc L3-5 RORB LINC01202': 640, 'Inh L1-5 VIP SMOC1': 638, 'Inh L5 PVALB LRIG3': 631, 'Exc L3-5 FEZF2 ASGR2': 620, 'Inh L5-6 PVALB ZFPM2-AS1': 563, 'Oligo L3-6 OPALIN ENPP6': 535, 'Exc L3-5 RORB TNNT2': 535, 'Inh L1-2 VIP EXPH5': 509, 'Exc L6 THEMIS SNTG2': 488, 'Inh L1-6 LAMP5 CA1': 465, 'Inh L1-3 VIP CBLN1': 464, 'Inh L2-5 PVALB HHIPL1': 456, 'Exc L6 THEMIS SLN': 436, 'Inh L2-3 SST NMU': 419, 'Inh L1-6 LAMP5 NES': 419, 'Inh L1-3 VIP HSPB6': 396, 'Inh L5-6 SST ISX': 394, 'Inh L3-5 VIP IGDCC3': 390, 'Astro L1-6 FGFR3 PLCG1': 377, 'Inh L1-2 SST PRRT4': 376, 'Inh L3-5 VIP TAC3': 367, 'Inh L1 LAMP5 BMP2': 352, 'Exc L5-6 THEMIS SMYD1': 327, 'Inh L1-2 SST CCNJL': 320, 'Exc L5-6 FEZF2 FILIP1L': 318, 'Inh L1-3 VIP CHRNA2': 310, 'Exc L5 RORB MED8': 309, 'Inh L2 VIP SLC6A16': 288, 'Inh L1-6 PVALB COL15A1': 288, 'OPC L1-6 PDGFRA COL20A1': 283, 'Exc L5 FEZF2 PKD2L1': 280, 'Exc L5-6 FEZF2 OR1L8': 263, 'Exc L2 LAMP5 KCNG3': 261, 'Inh L3-5 PVALB ISG20': 254, 'Exc L3-5 RORB LAMA4': 250, 'Inh L1 LAMP5 RAB11FIP1': 246, 'Inh L3-5 SST OR5AH1P': 242, 'Exc L5 THEMIS RGPD6': 241, 'Inh L5-6 PVALB MEPE': 230, 'Exc L5 FEZF2 NREP-AS1': 229, 'Inh L3-5 VIP HS3ST3A1': 225, 'Inh L1 SST DEFB108B': 223, 'Inh L5-6 PVALB SST CRHR2': 218, 'Exc L5 THEMIS VILL': 216, 'Inh L1-2 SST CLIC6': 194, 'Inh L5-6 SST PAWR': 194, 'Inh L1 LAMP5 NMBR': 193, 'Inh L5-6 SST KLHL1': 188, 'Inh L1 PVALB SST ASIC4': 186, 'Inh L1 PAX6 CHRFAM7A': 185, 'Exc L2 LINC00507 ATP7B': 183, 'Exc L5-6 FEZF2 SH2D1B': 159, 'Exc L5 THEMIS LINC01116': 151, 'Exc L5 FEZF2 CSN1S1': 150, 'Exc L5 THEMIS FGF10': 147, 'Inh L3-6 VIP ZIM2-AS1': 146, 'Inh L1 PAX6 MIR101-1': 139, 'Inh L5-6 PVALB GAPDHP60': 138, 'Inh L3-6 VIP UG0898H09': 131, 'Inh L5-6 SST PIK3CD': 129, 'Inh L1-2 VIP SCML4': 128, 'Exc L6 FEZF2 POGK': 127, 'Inh L1-5 VIP LINC01013': 127, 'Inh L1-2 PVALB CDK20': 123, 'Astro L1-6 FGFR3 AQP1': 119, 'Exc L5-6 FEZF2 LPO': 117, 'Inh L2-5 VIP SOX11': 116, 'Inh L5-6 PVALB FAM150B': 111, 'Inh L1-3 VIP FNDC1': 111, 'Exc L5 FEZF2 RNF144A-AS1': 108, 'Exc L5-6 FEZF2 CFTR': 108, 'Micro L1-6 TYROBP CD74': 108, 'Inh L5-6 PVALB KCNIP2': 108, 'Inh L2 PAX6 FREM2': 108, 'Inh L2 PVALB FRZB': 107, 'Exc L6 FEZF2 FFAR4': 106, 'Oligo L2-6 OPALIN MAP6D1': 101, 'Inh L5-6 SST FBN2': 100, 'Inh L1 VIP KLHDC8B': 95, 'Inh L1 LAMP5 PVRL2': 94, 'Inh L5-6 SST C4orf26': 93, 'Inh L1-5 VIP PHLDB3': 92, 'Inh L5-6 VIP COL4A3': 89, 'Exc L3-5 FEZF2 LINC01107': 88, 'Inh L2-5 VIP BSPRY': 85, 'Inh L1-2 VIP PTGER3': 84, 'Inh L6 SST TH': 84, 'Inh L3-6 PAX6 LINC01497': 84, 'Inh L3-5 SST CDH3': 81, 'Inh L5-6 SST BEAN1': 76, 'Exc L6 FEZF2 PROKR2': 72, 'Astro L1 FGFR3 SERPINI2': 72, 'Inh L1-6 SST NPY': 67, 'Inh L1-5 VIP CD27-AS1': 67, 'Endo L2-5 NOSTRIN SRGN': 64, 'Inh L1-6 VIP SLC7A6OS': 61, 'Inh L1-2 VIP WNT4': 52, 'Inh L1-2 VIP HTR3A': 42, 'VLMC L1-5 PDGFRA COLEC12': 40, 'Inh L1 SST P4HA3': 38, 'Inh L5-6 SST DNAJC14': 35, 'Exc L6 FEZF2 PDYN': 22, 'Oligo L5-6 OPALIN LDLRAP1': 6}
    cell_type_alias_color = "cell_type_alias_color"
    cell_type_alias_order = "cell_type_alias_order"
    cell_type_alt_alias_label = "cell_type_alt_alias_label"  # {'L6 CT_2': 2602, 'Oligo_1': 535, 'Astro_2': 377, 'L6 IT Car3': 327, 'Chandelier': 288, 'OPC': 283, 'Micro-PVM': 108, 'Sncg_2': 108, 'Lamp5_1': 94, 'L5 ET_2': 88, 'Sst_7': 84, 'Sst Chodl': 67, 'Endo': 64, 'VLMC': 40}
    cell_type_alt_alias_color = "cell_type_alt_alias_color"
    cell_type_alt_alias_order = "cell_type_alt_alias_order"
    cell_type_designation_label = "cell_type_designation_label"  # {'Neuron 75': 12485, 'Neuron 90': 4631, 'Neuron 78': 4044, 'Neuron 80': 3331, 'Neuron 64': 3253, 'Neuron 101': 2602, 'Non-neuron 3': 2300, 'Neuron 86': 2134, 'Neuron 88': 1824, 'Neuron 108': 1801, 'Neuron 81': 1628, 'Neuron 3': 1539, 'Neuron 65': 1326, 'Neuron 91': 1231, 'Neuron 76': 1211, 'Neuron 46': 1154, 'Neuron 89': 1132, 'Neuron 77': 1088, 'Neuron 51': 1069, 'Neuron 6': 923, 'Neuron 94': 905, 'Neuron 82': 785, 'Neuron 42': 764, 'Neuron 116': 753, 'Neuron 87': 640, 'Neuron 29': 638, 'Neuron 62': 631, 'Neuron 111': 620, 'Neuron 60': 563, 'Non-neuron 2': 535, 'Neuron 83': 535, 'Neuron 20': 509, 'Neuron 96': 488, 'Neuron 5': 465, 'Neuron 19': 464, 'Neuron 63': 456, 'Neuron 95': 436, 'Neuron 47': 419, 'Neuron 4': 419, 'Neuron 21': 396, 'Neuron 43': 394, 'Neuron 28': 390, 'Non-neuron 8': 377, 'Neuron 52': 376, 'Neuron 39': 367, 'Neuron 14': 352, 'Neuron 97': 327, 'Neuron 50': 320, 'Neuron 102': 318, 'Neuron 25': 310, 'Neuron 92': 309, 'Neuron 26': 288, 'Neuron 72': 288, 'Non-neuron 1': 283, 'Neuron 113': 280, 'Neuron 98': 263, 'Neuron 73': 261, 'Neuron 61': 254, 'Neuron 84': 250, 'Neuron 2': 246, 'Neuron 45': 242, 'Neuron 100': 241, 'Neuron 71': 230, 'Neuron 114': 229, 'Neuron 30': 225, 'Neuron 18': 223, 'Neuron 44': 218, 'Neuron 85': 216, 'Neuron 68': 194, 'Neuron 48': 194, 'Neuron 15': 193, 'Neuron 57': 188, 'Neuron 16': 186, 'Neuron 9': 185, 'Neuron 74': 183, 'Neuron 103': 159, 'Neuron 99': 151, 'Neuron 110': 150, 'Neuron 93': 147, 'Neuron 35': 146, 'Neuron 8': 139, 'Neuron 69': 138, 'Neuron 38': 131, 'Neuron 49': 129, 'Neuron 23': 128, 'Neuron 109': 127, 'Neuron 37': 127, 'Neuron 66': 123, 'Non-neuron 6': 119, 'Neuron 117': 117, 'Neuron 34': 116, 'Neuron 24': 111, 'Neuron 70': 111, 'Neuron 59': 108, 'Neuron 11': 108, 'Neuron 107': 108, 'Non-neuron 11': 108, 'Neuron 115': 108, 'Neuron 67': 107, 'Neuron 105': 106, 'Neuron 0': 101, 'Neuron 56': 100, 'Neuron 27': 95, 'Neuron 1': 94, 'Neuron 55': 93, 'Neuron 36': 92, 'Neuron 31': 89, 'Neuron 112': 88, 'Neuron 33': 85, 'Neuron 58': 84, 'Neuron 7': 84, 'Neuron 22': 84, 'Neuron 41': 81, 'Neuron 53': 76, 'Non-neuron 7': 72, 'Neuron 106': 72, 'Neuron 32': 67, 'Neuron 40': 67, 'Non-neuron 9': 64, 'Neuron 10': 61, 'Neuron 13': 52, 'Neuron 12': 42, 'Non-neuron 10': 40, 'Neuron 17': 38, 'Neuron 54': 35, 'Neuron 104': 22, 'Non-neuron 4': 6}
    cell_type_designation_color = "cell_type_designation_color"
    cell_type_designation_order = "cell_type_designation_order"
    external_donor_name_label = "external_donor_name_label"  # {'H18.30.001': 42728, 'H18.30.002': 33805}
    external_donor_name_color = "external_donor_name_color"
    external_donor_name_order = "external_donor_name_order"
    specimen_type = "specimen_type"  # {'nucleus': 76533}
    full_genotype_label = "full_genotype_label"  # {}
    outlier_call = "outlier_call"  # {False: 76533}
    outlier_type = "outlier_type"  # {}

    assert (set(df_meta.columns) == set(v for v in locals() if not v.startswith("__")))
    assert df_meta.cluster_label.equals(df_meta.cell_type_alias_label)


celltypes = df_meta.cluster_label.str.split(" ", n=2, expand=True).rename(
    columns={0: 'type', 1: 'layer', 2: 'markers'}
)
# sample_name                           type layer          markers
# AAACCCAAGGATTTCC-LKTX_190129_01_A01    Inh  L1-2        SST CCNJL
# AAACCCAAGTATGGCG-LKTX_190129_01_A01    Exc  L5-6   FEZF2 IFNG-AS1
# AAACCCACAAAGTGTA-LKTX_190129_01_A01    Exc  L3-5   RORB LINC01202
# AAACCCACACTACTTT-LKTX_190129_01_A01    Exc    L2  LINC00507 GLRA3
# AAACCCACAGTGAGCA-LKTX_190129_01_A01  Oligo  L2-6    OPALIN FTH1P3
# ...                                    ...   ...              ...
# TTTGTTGAGATGGCGT-LKTX_190130_01_H01  Oligo  L2-6    OPALIN FTH1P3
# TTTGTTGCACAGCCAC-LKTX_190130_01_H01    Exc  L3-5        RORB LNX2
# TTTGTTGCAGAGACTG-LKTX_190130_01_H01    Exc  L2-3       RORB PTPN3
# TTTGTTGCATAATGAG-LKTX_190130_01_H01  Oligo  L2-6    OPALIN FTH1P3
# TTTGTTGTCTACTCAT-LKTX_190130_01_H01    Inh  L2-5     PVALB RPH3AL

markers = Counter(chain.from_iterable(celltypes.markers.str.split(' ')))
# {'RORB': 15900, 'LINC00507': 13756, 'GLRA3': 12485, 'LAMP5': 8536, 'PVALB': 7992, 'FEZF2': 7923, 'THEMIS': 6652, 'SST': 6644, 'VIP': 5013, 'LNX2': 4631, 'CARM1P1': 4044, 'CCDC68': 3331, 'RPH3AL': 3253, 'OPALIN': 2942, 'C9orf135-AS1': 2602, 'FTH1P3': 2300, 'OTOGL': 2134, 'SLC22A18': 1824, 'KLK7': 1801, 'PTPN3': 1628, 'AARD': 1539, 'SAMD13': 1326, 'RPRM': 1231, 'RTKN2': 1211, 'GGTLC3': 1154, 'TNFAIP6': 1132, 'DSG3': 1088, 'FAM20A': 1069, 'CRABP1': 923, 'LINC00343': 905, 'ENPEP': 785, 'RPL35AP11': 764, 'IFNG-AS1': 753, 'LINC01202': 640, 'SMOC1': 638, 'LRIG3': 631, 'ASGR2': 620, 'FGFR3': 568, 'ZFPM2-AS1': 563, 'ENPP6': 535, 'TNNT2': 535, 'PAX6': 516, 'EXPH5': 509, 'SNTG2': 488, 'CA1': 465, 'CBLN1': 464, 'HHIPL1': 456, 'SLN': 436, 'NES': 419, 'NMU': 419, 'HSPB6': 396, 'ISX': 394, 'IGDCC3': 390, 'PLCG1': 377, 'PRRT4': 376, 'TAC3': 367, 'BMP2': 352, 'SMYD1': 327, 'PDGFRA': 323, 'CCNJL': 320, 'FILIP1L': 318, 'CHRNA2': 310, 'MED8': 309, 'SLC6A16': 288, 'COL15A1': 288, 'COL20A1': 283, 'PKD2L1': 280, 'OR1L8': 263, 'KCNG3': 261, 'ISG20': 254, 'LAMA4': 250, 'RAB11FIP1': 246, 'OR5AH1P': 242, 'RGPD6': 241, 'MEPE': 230, 'NREP-AS1': 229, 'HS3ST3A1': 225, 'DEFB108B': 223, 'CRHR2': 218, 'VILL': 216, 'PAWR': 194, 'CLIC6': 194, 'NMBR': 193, 'KLHL1': 188, 'ASIC4': 186, 'CHRFAM7A': 185, 'ATP7B': 183, 'SH2D1B': 159, 'LINC01116': 151, 'CSN1S1': 150, 'FGF10': 147, 'ZIM2-AS1': 146, 'MIR101-1': 139, 'GAPDHP60': 138, 'UG0898H09': 131, 'PIK3CD': 129, 'SCML4': 128, 'LINC01013': 127, 'POGK': 127, 'CDK20': 123, 'AQP1': 119, 'LPO': 117, 'SOX11': 116, 'FAM150B': 111, 'FNDC1': 111, 'TYROBP': 108, 'CD74': 108, 'RNF144A-AS1': 108, 'KCNIP2': 108, 'CFTR': 108, 'FREM2': 108, 'FRZB': 107, 'FFAR4': 106, 'MAP6D1': 101, 'FBN2': 100, 'KLHDC8B': 95, 'PVRL2': 94, 'C4orf26': 93, 'PHLDB3': 92, 'COL4A3': 89, 'LINC01107': 88, 'BSPRY': 85, 'TH': 84, 'PTGER3': 84, 'LINC01497': 84, 'CDH3': 81, 'BEAN1': 76, 'PROKR2': 72, 'SERPINI2': 72, 'CD27-AS1': 67, 'NPY': 67, 'NOSTRIN': 64, 'SRGN': 64, 'SLC7A6OS': 61, 'WNT4': 52, 'HTR3A': 42, 'COLEC12': 40, 'P4HA3': 38, 'DNAJC14': 35, 'PDYN': 22, 'LDLRAP1': 6}

(df_meta.cluster_label[df_meta.class_label == 'Non-Neuronal'].value_counts())
# Oligo L2-6 OPALIN FTH1P3     2300
# Oligo L3-6 OPALIN ENPP6       535
# Astro L1-6 FGFR3 PLCG1        377
# OPC L1-6 PDGFRA COL20A1       283
# Astro L1-6 FGFR3 AQP1         119
# Micro L1-6 TYROBP CD74        108
# Oligo L2-6 OPALIN MAP6D1      101
# Astro L1 FGFR3 SERPINI2        72
# Endo L2-5 NOSTRIN SRGN         64
# VLMC L1-5 PDGFRA COLEC12       40
# Oligo L5-6 OPALIN LDLRAP1       6

(df_meta.cluster_label[df_meta.class_label == "Glutamatergic"].value_counts())
assert all(df_meta.cluster_label[df_meta.class_label == "Glutamatergic"].str.slice(stop=3) == "Exc")
# Exc L2 LINC00507 GLRA3         12485
# Exc L3-5 RORB LNX2              4631
# Exc L3 LAMP5 CARM1P1            4044
# Exc L2-3 RORB CCDC68            3331
#                                  ...
# Exc L3-5 FEZF2 LINC01107          88
# Exc L6 FEZF2 PROKR2               72
# Exc L6 FEZF2 PDYN                 22

(df_meta.cluster_label[df_meta.class_label == "GABAergic"].value_counts())
assert all(df_meta.cluster_label[df_meta.class_label == "GABAergic"].str.slice(stop=3) == "Inh")
# Inh L2-5 PVALB RPH3AL    3253
# Inh L1-6 LAMP5 AARD      1539
# Inh L3 PVALB SAMD13      1326
# Inh L3-5 SST GGTLC3      1154
#                           ...
# Inh L1-2 VIP HTR3A         42
# Inh L1 SST P4HA3           38
# Inh L5-6 SST DNAJC14       35
