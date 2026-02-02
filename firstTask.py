from pdbloader import *
import os
import Bio.PDB
from Bio import SeqIO, Seq
from Bio.Align.Applications import ClustalwCommandline
from Bio import Align, AlignIO
import re
from ColToPosToCol import *
import random
import numpy
from nexus import matrix2nexus

random.seed(10)

def euclidean(v1, v2):
    return sum((p-q)**2 for p, q in zip(v1, v2)) ** .5

def calc_evo_dstmtx(distance_matrices):
    n=len(distance_matrices)
    distance_matrix = np.zeros((n, n))
    for i in range(n):
            j=i+1
            while j<n:
                if j > i:
                    distance_matrix[i, j] = euclidean(distance_matrices[i], distance_matrices[j])
                    distance_matrix[j, i] = distance_matrix[i,j] 
                j+=1
    return distance_matrix
   

def ColPicker(alignments):
    length=len(alignments[0])
    columns=list(range(length))
    random.shuffle(columns)
    num_sampl=10
    rndColumns=[]
    i=0 
    j=0
    while(i<num_sampl):
        
        if not alignments[:,columns[j]].isalpha():
            print(alignments[:,columns[j]],columns[j],"UPPS")
            j+=1
            continue
        rndColumns.append(columns[j])
        i+=1
        j+=1
    return rndColumns

def build_dstmtxs(rndColumns,alignments, structures,pdbcodes):
    distance_mtxs=[]
    for i in range(len(alignments)):
        atoms=[]
        print("Dict:",pdbcodes[alignments[i].id])
        calphas=get_calphas(structures[pdbcodes[alignments[i].id]])
        for col in rndColumns:
            print(col)
            print(calphas[posColIdx[i][1][col]].get_parent())
            atoms.append(calphas[posColIdx[i][1][col]])
        distance_mtxs.append(calc_distance_matrix(atoms,flatten=True))
    return distance_mtxs

def IndexBuilder(alignments):
    posColIdx=[]
    for MSAseq in alignments:
        posColIdx.append(IndexConverter(MSAseq))
    return posColIdx

# fasta_out="/home/joscha/DataQtools/output/bd85b277-d586-4eae-9e21-7d6f05162788.fasta"
pdb_dir ="/data/joscha/Data/6df16995-6ece-4dc0-92c3-f776066b0bfb"
fasta_out=os.path.dirname(pdb_dir)+"/output/"+os.path.basename(pdb_dir)+".fasta"
print(fasta_out )
os.makedirs(os.path.dirname(fasta_out),exist_ok=True)

structures=[]
sequences=[]
pdbcodes=dict()

regex= re.compile("(?<=AF-)[^-]*")
z=0

for file in os.scandir(pdb_dir):
    if file.path.endswith(".pdb"):
        print("File_Found:", file.path)
        structures.append(read_pdb(re.findall(regex,repr(file))[0],file.path))
        pdbcodes[re.findall(regex,repr(file))[0]]=z
        z+=1
print(pdbcodes)
print(list(pdbcodes.keys()))
seqrecords=extract_multi_seqrecords(list(pdbcodes.keys()),structures)
print(seqrecords)
SeqIO.write(seqrecords, fasta_out, "fasta")



clustalw_exe = "clustalw2" 
input_fasta = fasta_out

clustalw_cline = ClustalwCommandline(
    clustalw_exe,
    infile=input_fasta
)

stdout, stderr = clustalw_cline()
print(stdout)

alignments = AlignIO.read(fasta_out.replace(".fasta",".aln"), "clustal")

print(type(alignments))
# alignment = next(alignments)






# length=len(alignments[0])
# columns=range(length)
# num_sampl=10
# rndColumns=random.sample(columns, num_sampl)
# columns=remove_Columns(rndColumns, columns)
# i=0 
# while(i<num_sampl):
#     for alig in alignments:
#         if alig[rndColumns[i]].isalpha():
#             continue
#         else:
#             rndColumns[i]=random.choice(columns)
#             columns[rndColumns[i]]=-1
#             break
#     i+=1
    


#print(get_calphas(structures[0])[posColIdx[0][1][296]].get_full_id())


# print(len(distance_mtxs), len(distance_mtxs[1]))
# print(distance_mtxs[0],"\n",distance_mtxs[1],"\n",distance_mtxs[2],"\n",distance_mtxs[3])
posColIdx=IndexBuilder(alignments)
rndColumns=ColPicker(alignments)
print(rndColumns)
distance_mtxs=build_dstmtxs(rndColumns,alignments,structures, pdbcodes)
print(distance_mtxs[0])
print(distance_mtxs[1])
evo_DST_MTX=calc_evo_dstmtx(distance_mtxs)
print(evo_DST_MTX)
ids=[i.id for i in alignments]
codes2speci={"O14746":"Homo sapiens" }
matrix2nexus(evo_DST_MTX,taxa=ids,nexusfile="./TERT.nex",imgformat='png',plot_now=True,cmdfile=False,cmdmode='w',splitstree_location='/splitstree/app/SplitsTree')