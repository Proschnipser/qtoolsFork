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
import plotly.express as px
import matplotlib
#matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.manifold import TSNE
#matplotlib.use('Agg')
import plotly
import webbrowser
import subprocess



def find_best_perplexity(X, perplexity_range=(5, 50), n_steps=6, n_iter=300, random_state=42):
    """
    Perplexity sweep.
    Tries different perplexities with the aim to select the best one based on lowest divergence
    @param X: scaled input X for t-SNE analysis
    @type X: numpy array
    """
    n_samples = X.shape[0]
    p_min = max(2, perplexity_range[0])
    p_max = min(perplexity_range[1], (n_samples - 1) / 3)

    perplexities = sorted(set(np.linspace(p_min, p_max, n_steps).astype(int)))

    print(f"Testing {len(perplexities)} perplexities: {perplexities}\n")
    print(f"{'Perplexity':>12}  {'KL Divergence':>14}")
    print("-" * 30)

    results = []
    for p in perplexities:
        kl = TSNE(perplexity=p, max_iter=n_iter, random_state=random_state).fit(X).kl_divergence_
        results.append((p, kl))
        print(f"{p:>12}  {kl:>14.4f}")

    best_perplexity = min(results, key=lambda x: x[1])[0]
    print(f"\nBest perplexity: {best_perplexity}")
    return best_perplexity

def plot_PCA_TSNE(mtx_vectors, names, outfile_prefix):
    # Convert to array
    X = np.array(mtx_vectors)

    # Standardize
    X_scaled = StandardScaler().fit_transform(X)

    # # Compute first 4 PCs
    # n_pcs = 4
    # pca = PCA(n_components=n_pcs)
    # X_pca = pca.fit_transform(X_scaled)

    # # Create pairwise plots
    # pc_pairs = list(itertools.combinations(range(n_pcs), 2))

    # fig, axes = plt.subplots(
    #     nrows=int(np.ceil(len(pc_pairs) / 2)),
    #     ncols=2,
    #     figsize=(12, 4 * int(np.ceil(len(pc_pairs) / 2)))
    # )

    # axes = np.array(axes).flatten()

    # for ax, (i, j) in zip(axes, pc_pairs):
    #     for category in set(labels):
    #         mask = np.array(labels) == category
    #         ax.scatter(
    #             X_pca[mask, i],
    #             X_pca[mask, j],
    #             label=category,
    #             alpha=0.7
    #         )

    #     ax.set_xlabel(
    #         f"PC{i+1} ({pca.explained_variance_ratio_[i]:.1%})"
    #     )
    #     ax.set_ylabel(
    #         f"PC{j+1} ({pca.explained_variance_ratio_[j]:.1%})"
    #     )
    #     ax.grid(True)

    # # Remove unused axes
    # for ax in axes[len(pc_pairs):]:
    #     fig.delaxes(ax)

    # handles, labels_ = axes[0].get_legend_handles_labels()
    # fig.legend(handles, labels_, loc="upper right")

    # plt.tight_layout()
    # plt.savefig(outfile_prefix+"pca_plot.png", dpi=300, bbox_inches="tight")
    # plt.close()

    # ------------------
    # t-SNE
    # ------------------

    # Adjust perplexity depending on dataset size

    tsne = TSNE(perplexity=10)#find_best_perplexity(X_scaled))

    X_tsne = tsne.fit_transform(X_scaled)

    print(X_tsne.shape)
    print(X_tsne[:5])
    df = pd.read_csv("/data/joscha/Data/uniprot_2025_01_08_MIAs_02_MOTHs.csv", index_col='no.')
    regex= re.compile(r"_(\d+)_")
    labels=names.copy()
    for i,label in enumerate(names):
        index= int(re.search(regex, label).group(1))
        #print( index)
        #print(";".join(df.iloc[index]["OC"].split(";")[2:5]))
        labels[i]= ";".join(df.iloc[index]["OC"].split(";")[2:6])


    labels = np.array(labels)
    print(labels.shape)
    print(np.unique(labels))

#Then inside the loop:

    for category in np.unique(labels):
        mask = labels == category
        print(category, np.sum(mask))

    plt.figure(figsize=(8, 6))
    plt.grid(True)

    for category in np.unique(labels):
        mask = labels == category

        plt.scatter(
            X_tsne[mask, 0],
            X_tsne[mask, 1],
            alpha=0.7,
            label=category
        )

    plt.title("t-SNE")
    plt.xlabel("t-SNE 1")
    plt.ylabel("t-SNE 2")
    
    plt.legend(
    bbox_to_anchor=(1.05, 1),
    loc="upper left",
    borderaxespad=0.
    )

    
    plt.savefig(outfile_prefix+"tsne.png", dpi=300, bbox_inches="tight")
    plt.close()
    plt.tight_layout(rect=[0, 0, 0.85, 1])  # reserve 15% width for legend
    plt.show()
    
    
    colors  = ['#534AB7', '#1D9E75', '#D85A30', '#BA7517']


    df = pd.DataFrame(X_tsne, columns=['t-SNE 1', 't-SNE 2'])
    df['Class'] = labels
    df['Names'] = names
    
    fig = px.scatter(df, x='t-SNE 1', y='t-SNE 2',
                     color='Class',
                     title='Interactive t-SNE',
                     hover_data={'t-SNE 1': ':.2f', 't-SNE 2': ':.2f'},
                     hover_name='Names')
    fig.update_traces(marker=dict(size=9, opacity=0.85))
    filepath=outfile_prefix+"tsne.html"
    fig.write_html(filepath)  # opens in browser; zoom, pan, hover all work out of the box
    subprocess.Popen([
        "firefox",
        "--new-tab",
        filepath
    ])


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
train_names = sample([leaf.name for leaf in ref_tree.get_terminals()],35)
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
    alignment.id=re.findall(regx,str(alignment.id))[0]


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
columnsofinterest=[i for i in range(alignments.get_alignment_length()) if "-" not in str(alignments[:,i])]
print("Columnsofinterest:", columnsofinterest)
mtx_vectors=build_dstmtxs(columnsofinterest,alignments,structures, names_struct)
#print(mtx_vectors,len(mtx_vectors), IDs, train_names)
print("mtx_vector:",len(mtx_vectors[0]), len(columnsofinterest))

plot_PCA_TSNE(mtx_vectors, list(train_names), outfile_prefix)

vectors_df=pd.DataFrame({'spec': train_names, 'mtxvector': mtx_vectors})
vectors_df.to_csv(outfile_prefix+"vectors.csv", index=False)

# prune the reference tree to names from testing and write to file 
tree = dp.prune_tree(ref_tree, train_names)
dp.write_tree(tree, outfile_prefix+'tree.ph')

minibatches= dp.get_quartets_from_tree(tree)
#print(minibatches)
minibatches.to_csv(outfile_prefix  + 'minibatches.csv')

siamesebatches = itertools.combinations(train_names, 2)
pd.DataFrame(siamesebatches).to_csv(outfile_prefix  + 'siamesebatches.csv')

distances = dp.get_edgelength(tree, train_names)
distances.to_csv(outfile_prefix  + 'patristic.csv')