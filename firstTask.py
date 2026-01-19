from pdbloader import *
import os
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio import SeqIO, Seq
from Bio.Align.Applications import ClustalwCommandline
from Bio import Align
import re
from ColToPosToCol import *
import random


# fasta_out="/home/joscha/DataQtools/output/bd85b277-d586-4eae-9e21-7d6f05162788.fasta"
pdb_dir="/data/joscha/Data/6df16995-6ece-4dc0-92c3-f776066b0bfb"
fasta_out=os.path.dirname(pdb_dir)+"/output/"+os.path.basename(pdb_dir)+".fasta"
print(fasta_out )
os.makedirs(os.path.dirname(fasta_out),exist_ok=True)

structures=[]
sequences=[]
#seqrecords=[]
pdbcodes=[]

regex= re.compile("(?<=AF-)[^-]*")
z=0

for file in os.scandir(pdb_dir):
    if file.path.endswith(".pdb"):
        print("File_Found:", file.path)
        structures.append(read_pdb(re.findall(regex,repr(file))[0],file.path))
        pdbcodes.append(re.findall(regex,repr(file))[0])
        #seqrecords.append(extract_seqrecords(re.findall(regex,repr(file))[0],structures[-1]))
        z+=1
print(pdbcodes)
seqrecords=extract_multi_seqrecords(pdbcodes,structures)
print(seqrecords)
SeqIO.write(seqrecords, fasta_out, "fasta")



clustalw_exe = "clustalw2"  # or "clustalw"
input_fasta = fasta_out

clustalw_cline = ClustalwCommandline(
    clustalw_exe,
    infile=input_fasta
)

stdout, stderr = clustalw_cline()
print(stdout)

alignments = Align.read(fasta_out.replace(".fasta",".aln"), "clustal")

# print(alignments.metadata)
# alignment = next(alignments)

print(alignments)
# posColIdx=[]
# for MSAseq in alignments:
#     posColIdx.append(IndexConverter(MSAseq))

# print(posColIdx[0])
