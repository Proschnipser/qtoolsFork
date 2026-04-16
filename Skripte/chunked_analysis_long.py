import sys
import os
from pathlib import Path
from multiprocessing import Pool

def run_hmmsearch(args):
    hmmer_hits, file,model, old_hitfile, chunked = args
    modelname=model.stem.split("strict")[0]
    hits_out=hmmer_hits+"/"+chunked.stem+modelname+".tbl"
    print(hits_out)
    std_out=hits_out.replace(".tbl", ".out")
    print(old_hitfile, chunked)
    os.system(f"esl-translate {chunked} | hmmsearch --tblout {hits_out} {model} - > {std_out}")
    os.system(f"rm {str(old_hitfile).replace('.tbl','')}*")
    return

def seqkit_chunk(fnafile):
    chunked=Path(fnafile.replace(".fna","_chunked.fna"))
    command= f"seqkit sliding -s 240000 -W 270000 --greedy {fnafile} > {chunked}"
    print(command)
    os.system(command)

#directory=sys.argv[1] #/bioinf/data/user/joscha/Downloads/Cyclostomata/ncbi_dataset/ncbi_dataset/data
#cyclostomata="/data/joscha/Downloads/Cyclostomata/ncbi_dataset/data"
chondrichthyes="/data/joscha/Downloads/Chondrichthyes/ncbi_dataset/data"
hmmer_models=Path("/data/joscha/output/hmmer_models_long")
hmmer_hits="/data/joscha/output/hmmer_hits_long"
tasks=[]
genome_set=set()
for hitfile in Path(hmmer_hits).rglob("*.tbl"):
    #print(hitfile.is_file(), repr(str(hitfile)))
    if not hitfile.is_file() or hitfile.stat().st_size ==0:
        ID="_".join(hitfile.stem.split("_")[:2])
        print(ID, chondrichthyes+"/"+ID)
        fna_files = list(Path(chondrichthyes+"/"+ID).rglob("*genomic.fna"))
        print(fna_files)
        chunked=Path(str(fna_files[0]).replace(".fna","_chunked.fna"))
        #command=f"seqkit sliding -s 240000 -W 270000 --greedy {str(fna_files[0])}"
        #print(command)
        #os.system(command)
        genome_set.add(str(fna_files[0]))
        for file in fna_files:
            print(file)
            for model in hmmer_models.iterdir():
                if model.is_file():
                    tasks.append((hmmer_hits, file, model, hitfile,chunked))
# for fnafile in genome_set:
#     chunked=Path(fnafile.replace(".fna","_chunked.fna"))
#     command= f"seqkit sliding -s 240000 -W 270000 --greedy {fnafile} > {chunked}"
#     print(command)
#     os.system(command)


with Pool(processes=max(1, os.cpu_count() - 1)) as pool:
    pool.map(seqkit_chunk,list(genome_set))


with Pool(processes=max(1, os.cpu_count() - 1)) as pool:
    pool.map(run_hmmsearch,tasks)