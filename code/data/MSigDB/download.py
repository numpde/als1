# RA, 2021-03-10

from bugs import *
from twig import log
from contextlib import redirect_stdout
from datetime import timezone, datetime

# About:
# https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp

urls = [
    "http://bioinf.wehi.edu.au/MSigDB/v7.1/Mm.c7.all.v7.1.entrez.rds",
    "http://bioinf.wehi.edu.au/MSigDB/v7.1/Mm.c2.cp.reactome.v7.1.entrez.rds",
]


def worker(url, local):
    # # Reference implementation
    # import urllib.request
    # from shutil import copyfileobj
    # from contextlib import closing
    # with closing(urllib.request.urlopen(url=url)) as rd:
    #     with local.open(mode='wb') as fd:
    #         copyfileobj(rd, fd)

    # The following is based on
    # https://stackoverflow.com/questions/15644964/python-progress-bar-and-downloads

    import requests
    from tqdm import tqdm

    response = requests.get(url, stream=True)
    length = int(response.headers.get('content-length', default=0))

    with tqdm(desc="Downloading", total=length, unit="B", unit_scale=True, unit_divisor=1024) as tqdm:
        with local.open(mode='wb') as fd:
            for chunk in response.iter_content(chunk_size=1024):
                tqdm.update(fd.write(chunk))


def download(url):
    local = mkdir(Path(__file__).with_suffix('')) / Path(url).name

    log.info(f"Downloading")
    log.info(f"   from: {url}")
    log.info(f"     to: {relpath(local)}")

    if local.exists():
        log.info("File exists; skipping.")
        return

    worker(url, local)

    with Path(str(Path(local)) + "--meta.txt").open(mode='w') as fd:
        with redirect_stdout(fd):
            print(f"Downloaded ({datetime.now(tz=timezone.utc).isoformat(sep=' ')}) from")
            print(url)

    log.info("Done.")


def main():
    for url in urls:
        download(url)


if __name__ == '__main__':
    main()
