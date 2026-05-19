import pandas as pd
import  qtools.data_prepper as dp
import itertools
from pathlib import Path
from random import sample
import re
from natsort import natsorted, os_sorted, ns
from qtools.pdbloader import *
from qtools.firstTask import *
from Bio import SeqIO, Seq , Align, AlignIO
from qtools.ColToPosToCol import *

referencestructure= "/data/joscha/3D-Structures/boltz_results_MOTH_1_55_TGO1_HUMAN_Reviewed__1907_AA_TANGO1_homolog/predictions/MOTH_1_55_TGO1_HUMAN_Reviewed__1907_AA_TANGO1_homolog/MOTH_1_55_TGO1_HUMAN_Reviewed__1907_AA_TANGO1_homolog_model_0.cif"
pdb_files = "/data/joscha/3D-Structures/"
treefile = "/data/joscha/Downloads/SRw3UZCwUl830gEIOhHRkw_newick.tree"
ref_tree = dp.read_tree(treefile)
outfile_prefix = '/data/joscha/output/qtools/'+ str(Path(treefile).stem).replace(".tree","")+"/"
fasta_out=pdb_files+"output/"+Path(pdb_files).name+".fasta"
positions_of_interest=[x-23 for x in [54,56,63,85,98,100,101]]
print("Positions of interest:",positions_of_interest)

print(fasta_out )
Path(fasta_out).parent.mkdir(exist_ok=True)
Path(outfile_prefix).mkdir(parents=True, exist_ok=True)
for leaf in ref_tree.get_terminals():
    print(leaf.name)

# prepare names of species used for testing
train_names = sample([leaf.name for leaf in ref_tree.get_terminals()],10)
#_(\d+)_
IDs= [re.search(r'_(\d+)_', name).group(1) for name in train_names]

sorted_pairs = natsorted(zip(train_names, IDs), key=lambda pair: pair[1])
train_names, IDs = zip(*sorted_pairs)
print(train_names, IDs)

pattern= '|'.join(IDs)
regex = re.compile(rf'_({pattern})_\d+_')

matches = os_sorted([p for p in Path(pdb_files).iterdir() if regex.search(p.name)])
print("Matches:", matches, IDs)
structures=[]
names_struct={}
regx=re.compile(r"_(\d+)_\d+_")
z=0
structure_paths=""
for subdir in matches:
    for cif in subdir.rglob("*model_0.cif"):
        structure_paths+=str(cif)+" "
        identifier=  re.findall(regx,str(cif.name))[0]
        print("File_Found:", str(cif),identifier)
        structures.append(read_cif(identifier,str(cif)))
        names_struct[identifier]=z
        break
    z+=1
structure_paths+=referencestructure
print(names_struct)
print(list(names_struct.keys()))    


seqrecords=extract_multi_seqrecords(list(names_struct.keys()),structures)
print(seqrecords, len(seqrecords))
SeqIO.write(seqrecords, fasta_out, "fasta")
aln_pathclustal=fasta_out.replace(".fasta",f"_clustal.fa")
aln_path=fasta_out.replace(".fasta",f"_fdmason.fa")

tree_path=aln_path.replace(".fa",".dnd")
aln_path=fasta_out.replace(".fasta",f"_fdmason")
os.system(f"foldmason easy-msa {structure_paths} {aln_path}  tmpFolder --report-mode 2")
#os.system(f"clustalo -i {fasta_out} -o {aln_pathclustal} --outfmt=a2m --guidetree-out={tree_path} --force")

alignments = AlignIO.read(aln_path+"_aa.fa", "fasta")
print(type(alignments))
for alignment in alignments:
    alignment.id=re.findall(regx,str(alignment))[0]


sorted_records = natsorted(alignments, key=lambda r: r.id, reverse=True)
alignments = Align.MultipleSeqAlignment(sorted_records)
print(alignments, type(alignments))
posColIdx=IndexBuilder(alignments)

#rndColumns=ColPicker(alignments)
columnsofinterest=[posColIdx[-1][0][p] for p in positions_of_interest]
print(str(alignments[-1]))
print([alignments[-1][column] for column in columnsofinterest])
print(columnsofinterest)
alignments=alignments[:-1]
print(alignments)
mtx_vectors=build_dstmtxs(columnsofinterest,alignments,structures, names_struct)
print(mtx_vectors,len(mtx_vectors), IDs, train_names)
print(mtx_vectors)
vectors_df=pd.DataFrame({'spec': train_names, 'mtxvector': mtx_vectors})
vectors_df.to_csv(outfile_prefix+"vectors.csv", index=False)

# prune the reference tree to names from testing and write to file 
tree = dp.prune_tree(ref_tree, train_names)
dp.write_tree(tree, outfile_prefix+'tree.ph')

minibatches= dp.get_quartets_from_tree(tree)
print(minibatches)
minibatches.to_csv(outfile_prefix  + 'minibatches.csv')

siamesebatches = itertools.combinations(train_names, 2)
pd.DataFrame(siamesebatches).to_csv(outfile_prefix  + 'siamesebatches.csv')

distances = dp.get_edgelength(tree, train_names)
distances.to_csv(outfile_prefix  + 'patristic.csv')