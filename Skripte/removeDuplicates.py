import Bio
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import os
import sys
filepath=sys.argv[1] #/data/joscha/Data/uniprot_2025_01_08_MIAs_02_MOTHs.csv
#df=pd.read_csv("/data/joscha/Data/TANGO1onlySP.csv")
df=pd.read_csv(filepath)
unique={}
indices=[]
for i,r in df.iterrows():
    seq=r["Sequence cutted"]
    if seq in unique and r["OS"] == unique[seq]: # check wether sequence is unique and if not if its the same species it appears in.
        indices.append(i)
    else:
        unique[seq]=r["OS"]
df.drop(indices, inplace=True)
print(len(unique))
outputpath=os.path.splitext(filepath)[0]+"_dedup.csv"
print(outputpath)
df.to_csv(outputpath, index=False) # return deduplicated csv file.