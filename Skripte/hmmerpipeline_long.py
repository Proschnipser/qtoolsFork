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
fna_files = list(Path(directory).rglob("*genomic.fna"))

def run_hmmsearch(args):
    hmmer_hits,file, model=args
    modelname=model.stem.split("strict")[0]
    hits_out=hmmer_hits+"/"+file.stem+modelname+".tbl"
    std_out=hits_out.replace(".tbl", ".out")
    print(hits_out)
    if not Path(hits_out).is_file() or Path(hits_out).stat().st_size ==0:
        os.system(f"esl-translate {file} | hmmsearch --tblout {hits_out} {model} - > {std_out}")
        print("Finished:", hits_out)

tasks=[]
commands=[]
for file in fna_files:
    print(file)
    chunked=Path(str(file).replace(".fna","_chunked.fna"))
    command=f"seqkit sliding -s 240000 -W 270000 {str(file)} > {chunked}"
    os.system(command)
    #commands.append(command)
    for model in hmmer_models.iterdir():
        if model.is_file():
            tasks.append((hmmer_hits, chunked,model))

# with Pool(processes=max(1, os.cpu_count() - 1)) as pool:
#     pool.map(os.system, commands)

with Pool(processes=max(1, os.cpu_count() - 1)) as pool:
    pool.map(run_hmmsearch, tasks)
            
# for e in os.scandir(directory):
#     if e.is_dir():
        