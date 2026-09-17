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
from pathlib import Path 


def extract_seqrecords(df, threshold_length):
    """
    Extracts cutted Sequences, ID(no.) and Accession Codes (AC),Organism Classification (OC)from dataframe df into SeqRecord objects return as a list.
    @param df: Data frame from csv with exoasy annotation
    @type: pandas dataframe
    @return: seqrecords
    @rtype: list
    """
    global num
    seqrecords=[]
    for i,r in df.iterrows(): 
        seq=Seq(r["Sequence cutted"])
        if threshold_length <= len(seq):
            seqid=str(r["no."])
            descr=r["OS"].replace(" ", "_").replace("(", "<").replace(")",">").replace(",","")
            seqrec = Bio.SeqRecord.SeqRecord(seq, id=descr+seqid, 
                        description=seqid)
            seqrecords.append(seqrec)
    return seqrecords

def removeDuplicatesAndReduceByGenus(df):
    """
    Removes rows with duplicates in column "Sequence cutted" or removes
    @param df: df: Data frame from csv with exoasy annotation
    @type: pandas dataframe
    """
    unique=set()
    genuslist=set()
    indices=[]
    for i,r in df.iterrows():
        seq=r["Sequence cutted"]
        if seq in unique or (r["OC"],r["name"]) in genuslist: # check wether sequence is unique
            indices.append(i)
        else:
            unique.add(seq)
            genuslist.add((r["OC"],r["name"]))
    df.drop(indices, inplace=True)
    return df

directory="/data/joscha/Data/MOTH_uniprot_filtered_labeled"
df_dict = dict()
reg=re.compile("confident_([^_]+)_.*?(sheet.*)")
for file in Path(directory).glob("*.csv"):
    name=re.search(reg, file.stem)
    print(name)
    if name:
        name="_".join(name.groups())
        print(name)
        df_dict[name]  = pd.read_csv(file)
namelist = ['MIA_sheet','OTOR_sheet_no_outlier','TALI_sheet_no_outlier','TANGO1_sheet_no_outlier', 'TANGO1_sheet_Helix_residues_no_outlier'] 
df=pd.concat([df_dict[name] for name in namelist])
df_dict['TANGO1_sheet_Helix_residues_no_outlier'].rename(columns={df.columns[0]: "no."}, inplace=True)
df_dict['TANGO1_sheet_no_outlier'].rename(columns={df.columns[0]: "no."}, inplace=True)
set1=set(df_dict['TANGO1_sheet_no_outlier']['no.'])
set2=set(df_dict['TANGO1_sheet_Helix_residues_no_outlier']['no.'])
print("Symmetric difference length:",len(set1 ^ set2))
print("Difference length:",len(set1-set2))
for name in namelist:
    length=df_dict[name]["Sequence"].str.len()
    print(name, "Min:",length.min(), "Median:", length.median(), "Mean:", length.mean(),"Max:", length.max())
    print("Sequence cutted:")
    length=df_dict[name]["Sequence cutted"].str.len()
    print(name, "Min:",length.min(), "Median:", length.median(), "Mean:", length.mean(),"Max:", length.max())
sys.exit()
print(len(df))
df= removeDuplicatesAndReduceByGenus(df)
threshold_length= 99
df= df[df["Sequence cutted"].str.len()>=threshold_length] 
print(len(df))
print(df["Sequence"].str.len().min() )
df.rename(columns={df.columns[0]: "no."}, inplace=True)
print(df["no."])

fasta_out=os.path.splitext(file)[0]+"_names.fasta"
seqrecords=extract_seqrecords(df, threshold_length)
print(len(set([x.seq for x in seqrecords])), len(seqrecords))
sample_size=len(seqrecords)
SeqIO.write(seqrecords, fasta_out, "fasta")
aln_path=fasta_out.replace(".fasta",f"_thresh{threshold_length}aa_{sample_size}.fa")
aln_nex= aln_path.replace(".fa",".nex")
#aln_path_phy= aln_path.replace(".fa",".phy")
tree_path=aln_path.replace(".fa",".dnd")
print(f"clustalo -i {fasta_out} -o {aln_path} --outfmt=a2m --guidetree-out={tree_path} --force")
os.system(f"clustalo -i {fasta_out} -o {aln_path} --outfmt=a2m --guidetree-out={tree_path} --force")
records = list(SeqIO.parse(aln_path, "fasta"))
for record in records:
    record.annotations["molecule_type"] = "protein"
    record.id = re.sub(r"\<.*\>", "", record.id)
    record.id = re.sub(r"[^a-zA-Z0-9_]", "_", record.id)
    record.name = record.id
    record.description = ""
#records=random.sample(records, k=sample_size)
alignment = MultipleSeqAlignment(records)
print(aln_path)
print(aln_nex)
print(f"Alignment has {len(alignment)} sequences")
SeqIO.write(alignment, aln_nex, "nexus")