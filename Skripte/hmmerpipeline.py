import sys
import os
from pathlib import Path


directory=sys.argv[1] #/bioinf/data/user/joscha/Downloads/Cyclostomata/ncbi_dataset/ncbi_dataset/data
hmmer_models=Path("/data/joscha/output/hmmer_models")
hmmer_hits="/data/joscha/output/hmmer_hits"
fna_files = list(Path(directory).rglob("*.fna"))
for file in fna_files:
    print(file)
    for model in hmmer_models.iterdir():
        if model.is_file():
            modelname=model.stem.split("strict")[0]
            hits_out=hmmer_hits+"/"+file.stem+modelname+".tbl"
            print(hits_out)
            os.system(f"esl-translate {file} | hmmsearch --tblout {hits_out} {model} -")
# for e in os.scandir(directory):
#     if e.is_dir():
        