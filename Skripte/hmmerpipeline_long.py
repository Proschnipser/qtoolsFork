import sys
import os
from pathlib import Path
from multiprocessing import Pool
import itertools

hmmer_models=Path("/data/joscha/output/hmmer_models_long")
# Types=["TANGO1","MIA","TALI","OTOR"]
# for protein in Types:
#     name=f"{hmmer_models}/{protein}strict_dedup_reduced_long_names.hmm /data/joscha/Data/{protein}strict_dedup_reduced_long_names.sto "
    
#     print(name)
#     os.system(f"hmmbuild {name}")
directory=sys.argv[1] #/bioinf/data/user/joscha/Downloads/Cyclostomata/ncbi_dataset/ncbi_dataset/data
hmmer_hits="/data/joscha/output/hmmer_hits_long"
fna_files = list(Path(directory).rglob("*.fna"))

def run_hmmsearch(hmmer_hits,file, model):
    modelname=model.stem.split("strict")[0]
    hits_out=hmmer_hits+"/"+file.stem+modelname+".tbl"
    std_out=hits_out.replace(".tbl", ".out")
    print(hits_out)
    os.system(f"esl-translate {file} | hmmsearch --tblout {hits_out} {model} - > {std_out}")

for file in fna_files:
    print(file)
    for model in hmmer_models.iterdir():
        if model.is_file():
            
# for e in os.scandir(directory):
#     if e.is_dir():
        