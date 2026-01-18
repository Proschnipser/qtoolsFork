from pdbloader import *
import os
from Bio.PDB.MMCIFParser import MMCIFParser
structures=[]
sequences=[]
for file in os.scandir("/home/joscha/DataQtools/bd85b277-d586-4eae-9e21-7d6f05162788"):
    if file.path.endswith(".cif"):
        print("File_Found:", file.path)
        structures.append(read_cif("O14746",file.path))
# struct=MMCIFParser("O14746",files.path)
print(structures)
#extract_seqrecords("7R3M",struct)[1]

# struct= read_pdb("O14746","/home/joscha/DataQtools/bd85b277-d586-4eae-9e21-7d6f05162788/AF-O14746-F1-model_v6.pdb")
# print(extract_seqrecords("O14746",struct))
# struct.get_se
