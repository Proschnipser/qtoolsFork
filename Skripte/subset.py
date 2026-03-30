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



filepath=sys.argv[1]
df=pd.read_csv(filepath)
#df=pd.read_csv("/data/joscha/Data/TANGO1onlySP_dedup.csv")
fasta_out=os.path.splitext(filepath)[0]+"_names.fasta"
threshold_length=100
old_sampl_size=474
sample_size=150
aln_path=fasta_out.replace(".fasta",f"_thresh{threshold_length}aa_{old_sampl_size}.fa")
aln_nex= aln_path.replace(f"aa_{old_sampl_size}.fa",f"aa_{sample_size}alt.nex")
#aln_path_phy= aln_path.replace(".fa",".phy")
tree_path=aln_path.replace(".fa",".dnd")
records = list(SeqIO.parse(aln_path, "fasta"))
for record in records:
    record.annotations["molecule_type"] = "protein"
    record.id = re.sub(r"\<.*\>", "", record.id)
    record.id = re.sub(r"[^a-zA-Z0-9_]", "_", record.id)
    record.name = record.id
    record.description = ""
records=random.sample(records, k=sample_size)
alignment = MultipleSeqAlignment(records)
print(aln_path)
print(aln_nex)
print(f"Alignment has {len(alignment)} sequences")
SeqIO.write(alignment, aln_nex, "nexus")


# gap_open=[2.5,10,15,20,30,40,60]
# gap_ext= [0.2,0.5,1,2,4,8]
# for i in range(len(gap_ext)):
#     for j in range(len(gap_open)):
#         aln_path=fasta_out.replace(".fasta",f"_aln_open{gap_open[j]}_ext{gap_ext[i]}.fa")
#         tree_path=aln_path.replace(".fa",".dnd")
#         command=f"clustalw2 -ALIGN -INFILE={fasta_out} -OUTFILE={aln_path} -OUTPUT=FASTA -NEWTREE={tree_path} -OUTORDER=INPUT -GAPOPEN={str(gap_open[j])} -GAPEXT={str(gap_ext[i])} -quiet"
#         print(command)
#         os.system(command)
#Bio.AlignIO.convert(aln_path, "fasta", aln_path_phy, "phylip")
