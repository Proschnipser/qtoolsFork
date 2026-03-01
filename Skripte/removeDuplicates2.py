import Bio
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import os
import re

df=pd.read_csv("/data/joscha/Data/TANGO1onlySP.csv")
unique={}
indices=[]
for i,r in df.iterrows():
    seq=r["Sequence cutted"]
    if seq in unique and r["OC"] == unique[seq]: # check wether sequence is unique and if not if its the same genus it appears in.
        indices.append(i)
    else:
        unique[seq]=r["OC"]
df.drop(indices, inplace=True)
print(len(unique))
df.to_csv("/data/joscha/Data/TANGO1onlySP_dedup_genus.csv", index=False)