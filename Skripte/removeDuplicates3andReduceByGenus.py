import Bio
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import os
import sys
filepath=sys.argv[1] #/data/joscha/Data/uniprot_2025_01_08_MIAs_02_MOTHs.csv
#df=pd.read_csv("/data/joscha/Data/TANGO1onlySP.csv")
df=pd.read_csv(filepath)
unique=set()
genuslist=set()
indices=[]
for i,r in df.iterrows():
    seq=r["Sequence"]
    if seq in unique or r["OS"] in genuslist: # check wether sequence is unique
        indices.append(i)
    else:
        unique.add(seq)
        genuslist.add(r["OS"])
df.drop(indices, inplace=True)
print(len(unique))
outputpath=os.path.splitext(filepath)[0]+"_dedup_reduced_long.csv"
print(outputpath)
df.to_csv(outputpath, index=False) # return deduplicated csv file.