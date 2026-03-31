import sys
import os
from pathlib import Path

hmmer_hits="/data/joscha/output/hmmer_hits"
for hitfile in Path(hmmer_hits).iterdir():
    if hitfile.stat().st_size ==0 and Path(str(hitfile).replace("genomic","genomic_chunked")).is_file():
        print(hitfile)
        os.system(f"rm {hitfile}")
