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
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from datetime import datetime

def plot_PCA(mtx_vectors,labels, outfile_prefix):
    # Convert to array
    X = np.array(mtx_vectors)

    # Standardize
    X_scaled = StandardScaler().fit_transform(X)

    # Compute first 4 PCs
    n_pcs = 4
    pca = PCA(n_components=n_pcs)
    X_pca = pca.fit_transform(X_scaled)

    # Create pairwise plots
    pc_pairs = list(itertools.combinations(range(n_pcs), 2))

    fig, axes = plt.subplots(
        nrows=int(np.ceil(len(pc_pairs) / 2)),
        ncols=2,
        figsize=(12, 4 * int(np.ceil(len(pc_pairs) / 2)))
    )

    axes = np.array(axes).flatten()

    for ax, (i, j) in zip(axes, pc_pairs):
        for category in set(labels):
            mask = np.array(labels) == category
            ax.scatter(
                X_pca[mask, i],
                X_pca[mask, j],
                label=category,
                alpha=0.7
            )

        ax.set_xlabel(
            f"PC{i+1} ({pca.explained_variance_ratio_[i]:.1%})"
        )
        ax.set_ylabel(
            f"PC{j+1} ({pca.explained_variance_ratio_[j]:.1%})"
        )
        ax.grid(True)

    # Remove unused axes
    for ax in axes[len(pc_pairs):]:
        fig.delaxes(ax)

    handles, labels_ = axes[0].get_legend_handles_labels()
    fig.legend(
    handles,
    labels_,
    loc="center left",
    bbox_to_anchor=(0.85, 0.5)  # outside right edge
    )

    plt.tight_layout(rect=[0, 0, 0.85, 1])  # reserve 15% width for legend
    plt.savefig(outfile_prefix+"pca_plot.png", dpi=300, bbox_inches="tight")


referencestructure= "/data/joscha/3D-Structures/boltz_results_MOTH_1_55_TGO1_HUMAN_Reviewed__1907_AA_TANGO1_homolog/predictions/MOTH_1_55_TGO1_HUMAN_Reviewed__1907_AA_TANGO1_homolog/MOTH_1_55_TGO1_HUMAN_Reviewed__1907_AA_TANGO1_homolog_model_0.cif"
pdb_files = "/data/joscha/3D-Structures/"
treefile = "/data/joscha/Downloads/SRw3UZCwUl830gEIOhHRkw_newick.tree"
ref_tree = dp.read_tree(treefile)
timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M")
outfile_prefix = '/data/joscha/output/qtools/'+ str(Path(treefile).stem).replace(".tree","")+"/"+timestamp+"/"
fasta_out=pdb_files+"output/"+Path(pdb_files).name+".fasta"
positions_of_interest=[x-23 for x in [54,56,63,85,98,100,101]]
print("Positions of interest:",positions_of_interest)

print(fasta_out )
Path(fasta_out).parent.mkdir(exist_ok=True)
Path(outfile_prefix).mkdir(parents=True, exist_ok=True)
for leaf in ref_tree.get_terminals():
    print(leaf.name)

# prepare names of species used for testing
train_names = [leaf.name for leaf in ref_tree.get_terminals()]
#_(\d+)_
IDs= [re.search(r'_(\d+)_', name).group(1) for name in train_names]

sorted_pairs = natsorted(zip(train_names, IDs), key=lambda pair: pair[1])
train_names, IDs = zip(*sorted_pairs)
print(train_names, IDs)

pattern= '|'.join(IDs)
regex = re.compile(rf'_({pattern})_\d+_')

matches = os_sorted([p for p in Path(pdb_files).iterdir() if regex.search(p.name)])
print("Matches:", matches, IDs)
train_names= [name for name in train_names for _ in range(5)]
structures=[]
names_struct={}
regx=re.compile(r"_(\d+)_\d+_")
z=0
structure_paths=""
for subdir in matches:
    for cif in subdir.rglob("*.cif"):
        structure_paths+=str(cif)+" "
        print(cif.stem[-1])
        identifier=  re.findall(regx,str(cif.name))[0]+"_"+(cif.stem[-1])
        print("File_Found:", str(cif),identifier)
        structures.append(read_cif(identifier,str(cif)))
        names_struct[identifier]=z
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
command=f"foldmason easy-msa {structure_paths} {aln_path}  tmpFolder --report-mode 2"
print(command)
os.system(command)
#os.system(f"clustalo -i {fasta_out} -o {aln_pathclustal} --outfmt=a2m --guidetree-out={tree_path} --force")


alignments = AlignIO.read(aln_path+"_aa.fa", "fasta")
print(type(alignments))
for alignment in alignments:
    alignment.id=re.findall(regx,str(alignment.id))[0]+"_"+alignment.id[-1]


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
plot_PCA(mtx_vectors, train_names, outfile_prefix)