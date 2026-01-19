from pdbloader import *
import os
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio import SeqIO, Seq
from Bio.Align.Applications import ClustalwCommandline
import re



fasta_out="/home/joscha/DataQtools/output/bd85b277-d586-4eae-9e21-7d6f05162788.fasta"
structures=[]
sequences=[]
#seqrecords=[]
pdbcodes=[]

regex= re.compile("(?<=AF-)[^-]*")
z=0

for file in os.scandir("/home/joscha/DataQtools/bd85b277-d586-4eae-9e21-7d6f05162788"):
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

#extract_seqrecords("7R3M",struct)[1]

# struct= read_pdb("O14746","/home/joscha/DataQtools/bd85b277-d586-4eae-9e21-7d6f05162788/AF-O14746-F1-model_v6.pdb")
# print(extract_seqrecords("O14746",struct))
# struct.get_se
