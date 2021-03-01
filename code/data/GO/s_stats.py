# RA, 2021-03-16

from bugs import *
from tcga.utils import from_iterable

organism = "mouse"

file = unlist1(Path(__file__).parent.glob(f"a_*/goa_{organism}_sym2cat.txt.gz"))

df = pd.read_csv(file, sep='\t')
categories = df.symbol.groupby(df.cat).agg(list).sort_index()

print(f"Summary of GO category size ({organism}):", categories.apply(len).describe().to_markdown(), sep='\n')
