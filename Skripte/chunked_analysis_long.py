import sys
import os
from pathlib import Path
from multiprocessing import Pool

def run_hmmsearch(args):
    hmmer_hits, file,model, command, old_hitfile, chunked = args
    modelname=model.stem.split("strict")[0]
    hits_out=hmmer_hits+"/"+chunked.stem+modelname+".tbl"
    print(hits_out)
    std_out=hits_out.replace(".tbl", ".out")
    #os.system(command > chunked)
    os.system(f"{command}|esl-translate - | hmmsearch --tblout {hits_out} {model} - > {std_out}")
    return

#directory=sys.argv[1] #/bioinf/data/user/joscha/Downloads/Cyclostomata/ncbi_dataset/ncbi_dataset/data
cyclostomata="/data/joscha/Downloads/Cyclostomata/ncbi_dataset/data"
chondrichthyes="/data/joscha/Downloads/Chondrichthyes/ncbi_dataset/data"
hmmer_models=Path("/data/joscha/output/hmmer_models_long")
hmmer_hits="/data/joscha/output/hmmer_hits_long"
tasks=[]
for hitfile in Path(hmmer_hits).rglob("*.tbl"):
    if not hitfile.is_file() or hitfile.stat().st_size ==0:
        ID="_".join(hitfile.stem.split("_")[:2])
        print(ID, cyclostomata+"/"+ID)
        fna_files = list(Path(cyclostomata+"/"+ID).rglob("*.fna")) + list(Path(chondrichthyes+"/"+ID).rglob("*.fna"))
        print(fna_files)
        chunked=Path(str(fna_files[0]).replace(".fna","_chunked.fna"))
        command=f"seqkit sliding -s 240000 -W 270000 {str(fna_files[0])}"
        #print(command)
        #os.system(command)
        
        for file in fna_files:
            print(file)
            for model in hmmer_models.iterdir():
                if model.is_file():
                    tasks.append((hmmer_hits, file, model, command, hitfile,chunked))



with Pool(processes=max(1, 10)) as pool:
    pool.map(run_hmmsearch,tasks)