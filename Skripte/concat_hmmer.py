import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import pandas as pd
import os
import sys
import re
import random

max_no= max(pd.read_csv("/data/joscha/Data/uniprot_2025_01_08_MIAs_02_MOTHs.csv")["no."])+1
print(max_no)
chunked=pd.read_csv("/data/joscha/Data/hmmer_results_chunked.csv")
not_chunked=pd.read_csv("/data/joscha/Data/hmmer_results_not_chunked.csv")

Cyclo_df= pd.read_csv("/data/joscha/Downloads/Cyclostomata.tsv", sep="\t",header=0)
Chond_df= pd.read_csv("/data/joscha/Downloads/Chondrichthyes.tsv", sep="\t",header=0)

hmmer_df=pd.concat([chunked,not_chunked])
print(hmmer_df)
hmmer_df['OS'] = pd.Series(dtype='string')
hmmer_df['no.'] = range(max_no, max_no + len(hmmer_df))
for i,r in hmmer_df.iterrows():
    rows= Cyclo_df.loc[Cyclo_df["Assembly Accession"]==r["name"]]
    if len(rows)==0:
        rows= Chond_df.loc[Chond_df["Assembly Accession"]==r["name"]]
    if len(rows)==0:
        continue
    print(rows["Organism Name"], print(r["name"]))
    hmmer_df.loc[i, "OS"]=rows["Organism Name"].item()
    print(rows["Organism Name"].item())
print(hmmer_df["OS"])
unique=set([])
indices=[]

hmmer_df=hmmer_df.loc[hmmer_df["Sequence cutted"].notnull()]
for i,r in hmmer_df.iterrows():
    seq=r["Sequence cutted"]
    if seq+r["OS"].split()[0] in unique: # check wether sequence is unique
        indices.append(i)
        print("duplicate:", seq)
    else:
        unique.add(seq+r["OS"].split()[0])
print("Before:",len(hmmer_df))
print(hmmer_df.loc[indices])
hmmer_df.drop(indices, inplace=True)
print("After:",len(hmmer_df))
hmmer_df.to_csv("/data/joscha/Data/hmmer_Cyclostomata_Chondrithyes.csv",  index=False)

