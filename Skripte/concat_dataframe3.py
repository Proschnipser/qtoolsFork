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



def extract_seqrecords(df, threshold_length, seqrecords):
    """
    Extracts cutted Sequences, ID(no.) and Accession Codes (AC),Organism Classification (OC)from dataframe df into SeqRecord objects return as a list.
    @param df: Data frame from csv with exoasy annotation
    @type: pandas dataframe
    @return: seqrecords
    @rtype: list
    """
    global num
    for i,r in df.iterrows(): 
        seq=Seq(r["Sequence cutted"])
        if threshold_length <= len(seq):
            seqid=str(r["no."])
            descr=r["OS"].replace(" ", "_").replace("(", "<").replace(")",">").replace(",","")
            seqrec = Bio.SeqRecord.SeqRecord(seq, id=descr+seqid, 
                        description=seqid)
            seqrecords.append(seqrec)
    return seqrecords



filepath=sys.argv[1]
df=pd.read_csv(filepath)
#df=pd.read_csv("/data/joscha/Data/TANGO1onlySP_dedup.csv")
sample_size=20
threshold_length=89
euteleostomi_OTOR=df[df["OC"].str.contains("Euteleostomi") & (df["name"]=="OTOR") & (df["Sequence cutted"].str.len()>=threshold_length)]
euteleostomi_OTOR=(
    euteleostomi_OTOR.groupby("Sequence cutted", group_keys=False)
    .apply(lambda g: g.sample(1))  # pick a random row per unique value
    .sample(n=sample_size, random_state=42)
)
print(len(euteleostomi_OTOR))
print(len(set(euteleostomi_OTOR["Sequence cutted"])))
euteleostomi_MIA=df[df["OC"].str.contains("Euteleostomi") & (df["name"]=="MIA") & (df["Sequence cutted"].str.len()>=threshold_length)]
euteleostomi_MIA=(
    euteleostomi_MIA.groupby("Sequence cutted", group_keys=False)
    .apply(lambda g: g.sample(1))  # pick a random row per unique value
    .sample(n=sample_size, random_state=42)
)
print(len(euteleostomi_MIA))
euteleostomi_TALI=df[df["OC"].str.contains("Euteleostomi") & (df["name"]=="TALI") & (df["Sequence cutted"].str.len()>=threshold_length)]
euteleostomi_TALI=(
    euteleostomi_TALI.groupby("Sequence cutted", group_keys=False)
    .apply(lambda g: g.sample(1))  # pick a random row per unique value
    .sample(n=sample_size, random_state=42)
)
print(len(euteleostomi_TALI))
print(len(set(euteleostomi_TALI["Sequence cutted"])))
Ecdysozoa=df[df["OC"].str.contains("Ecdysozoa") & (df["Sequence cutted"].str.len()>=threshold_length)]
Ecdysozoa=(
    Ecdysozoa.groupby("Sequence cutted", group_keys=False)
    .apply(lambda g: g.sample(1))  # pick a random row per unique value
    .sample(n=sample_size, random_state=42)
)
print(len(Ecdysozoa))
print(len(set(Ecdysozoa["Sequence cutted"])))

euteleostomi_TANGO1=df[df["OC"].str.contains("Euteleostomi") & (df["name"].str.contains("TANGO1")) & (df["Sequence cutted"].str.len()>=threshold_length)]
print("Length:", len(euteleostomi_TANGO1))
euteleostomi_TANGO1=(
    euteleostomi_TANGO1.groupby("Sequence cutted", group_keys=False)
    .apply(lambda g: g.sample(1))  # pick a random row per unique value
    .sample(n=sample_size, random_state=42)
)
print(len(euteleostomi_TANGO1))
print(len(set(euteleostomi_TANGO1["Sequence cutted"])))
Ecdysozoa=df[df["OC"].str.contains("Ecdysozoa") & (df["name"].str.contains("TANGO1")) & (df["Sequence cutted"].str.len()>=threshold_length)]
Ecdysozoa=(
    Ecdysozoa.groupby("Sequence cutted", group_keys=False)
    .apply(lambda g: g.sample(1))  # pick a random row per unique value
    .sample(n=min(len(set(Ecdysozoa["Sequence cutted"].values)),sample_size), random_state=42)
)
print(len(Ecdysozoa))
print(len(set(Ecdysozoa["Sequence cutted"])))

Spiralia=df[df["OC"].str.contains("Spiralia") & (df["name"].str.contains("TANGO1")) & (df["Sequence cutted"].str.len()>=threshold_length)]
Spiralia=(
    Spiralia.groupby("Sequence cutted", group_keys=False)
    .apply(lambda g: g.sample(1))  # pick a random row per unique value
    .sample(n=min(len(set(Spiralia["Sequence cutted"].values)),sample_size), random_state=42)
)
print(len(Spiralia))
Others= df[~df["OC"].str.contains("Spiralia|Ecdysozoa|Euteleostomi") & (df["name"].str.contains("TANGO1")) & (df["Sequence cutted"].str.len()>=threshold_length)]
Others=(
    Others.groupby("Sequence cutted", group_keys=False)
    .apply(lambda g: g.sample(1))  # pick a random row per unique value
    .sample(n=min(len(set(Others["Sequence cutted"].values)),sample_size), random_state=42)
)
print(len(Others))
df=pd.concat([euteleostomi_OTOR,euteleostomi_MIA,euteleostomi_TALI,euteleostomi_TANGO1, Ecdysozoa, Spiralia, Others])
print(len(df))
fasta_out=os.path.splitext(filepath)[0]+"_andHMMer"+"_names.fasta"
seqrecords=[]

seqrecords=extract_seqrecords(df, threshold_length, seqrecords)
print(len(set([x.seq for x in seqrecords])), len(seqrecords))
sample_size=len(seqrecords)
SeqIO.write(seqrecords, fasta_out, "fasta")
aln_path=fasta_out.replace(".fasta",f"_thresh{threshold_length}aa_{sample_size}_TANGO1.fa")
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

