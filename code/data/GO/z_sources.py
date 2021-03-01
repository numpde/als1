# RA, 2021-03-15

from tcga.utils import download
from bugs import *

download = download.to(abs_path=(Path(__file__).with_suffix('') / "download_cache"))


class URLS:
    goa = pd.Series({
        # ftp://ftp.ebi.ac.uk/pub/databases/GO/goa
        'mouse': "ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/MOUSE/goa_mouse.gaf.gz",
    })

    # http://geneontology.org/docs/download-ontology/
    obo = "http://purl.obolibrary.org/obo/go.obo"


def url2file(url) -> Path:
    if url in set(URLS.goa):
        import re
        organism = unlist1(re.findall(pattern=r"go/goa/([a-z]+)/", string=str.lower(url)))
        return unlist1(Path(__file__).parent.glob(f"*/goa_{organism}_sym2cat.txt.gz"))

    if (url == URLS.obo):
        return unlist1(Path(__file__).parent.glob("*/obo_go.txt.gz"))

    raise KeyError("Unrecognized URL.")
