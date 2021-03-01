# RA, 2021-03-16

"""
This file only prepares the downloader for the
“Mouse Whole Cortex and Hippocampus 10x” dataset
from Allen Brain, 2020.
https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-10x
https://celltypes.brain-map.org/rnaseq/mouse_ctx-hip_10x

The full download (done by a_download.py)
takes a few hours and requires over 4GB local storage.

10x protocol:
https://portal.brain-map.org/atlases-and-data/rnaseq/protocols-mouse-cortex-and-hippocampus#single_cell_sorting
https://www.protocols.io/view/10xv2-rnaseq-sample-processing-ynxfvfn/abstract
"""

from tcga.utils import download
from bugs import *

download = download.to(abs_path=(Path(__file__).with_suffix('') / "download_cache"))

URLS = {
    'expr': "https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hip_10x/matrix.csv",
    'meta': "https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hip_10x/metadata.csv",
}
