import pandas as pd
import sys
from pathlib import Path
from io import StringIO
from collections import defaultdict
import re

COLS = [
    'target_name', 'target_accession', 'query_name', 'query_accession',
    'e_value', 'score', 'bias',
    'e_value_best', 'score_best', 'bias_best',
    'exp', 'reg', 'clu', 'ov', 'env', 'dom', 'rep', 'inc',
    'description'
]
cyclostomata="/data/joscha/Downloads/Cyclostomata/ncbi_dataset/data"
chondrichthyes="/data/joscha/Downloads/Chondrichthyes/ncbi_dataset/data"
directory="/data/joscha/output/hmmer_hits_long/" #sys.argv[1]
dfdict=defaultdict(dict)
for filepath in Path(directory).rglob("*.tbl"):
    rows = []
    with open(filepath) as f:
        print(filepath)
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            rows.append(parts[:18] + [' '.join(parts[18:])])

    df = pd.DataFrame(rows, columns=COLS)

   

    # Extract full multi-word desc= value
    df['desc'] = df['description'].str.extract(r'desc=(.+)$')

    # Parse remaining key=value pairs (excluding desc=)
    def parse_kv(s):
        s = re.sub(r'\s*desc=.+$', '', s)
        return dict(re.findall(r'(\w+)=(\S+)', s))

    desc_df = df['description'].apply(parse_kv).apply(pd.Series)
    df = pd.concat([df.drop(columns='description'), desc_df], axis=1)

     # Cast float columns
    for col in ['e_value', 'score', 'bias', 'e_value_best', 'score_best', 'bias_best', 'exp']:
        df[col] = df[col].astype(float)
    for col in ['reg', 'clu', 'ov', 'env', 'dom', 'rep', 'inc', 'length']:
        if col in df.columns:
            df[col] = df[col].astype(int)
    splitname=filepath.stem.split("_")
    print("_".join(splitname[:2]),splitname[-1].replace("chunked","").replace("genomic","") )
    dfdict["_".join(splitname[:2])][filepath.stem]=df

orfdf=pd.DataFrame()
for k1,v1 in dfdict.items(): #iterate over genomes
    orfdict={}
    for k2, v2 in v1.items(): #iterate over protein types
        print(k2)
        for i, r in v2.iterrows():
            if r["e_value"] > 0.05:
                break
            if r["length"] > 99 and (not r["target_name"] in orfdict or r["e_value"] < orfdict[r["target_name"]][4]):
                print(list(r))
                orfdict[r["target_name"]]=list(r)
                orfdict[r["target_name"]].extend([k1,k2])
    print(orfdict)
    for orf, orflist in orfdict.items():
        fna_files = list(Path(cyclostomata+"/"+orflist[-2]).rglob("*.fna"))
        if fna_files == []:
                fna_files = list(Path(chondrichthyes+"/"+orflist[-2]).rglob("_chunked.fna"))
    
        


            
            

            
