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
from itertools import compress
import os,sys


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
    regex= re.compile(r"_(\d+)")
    labels=names.copy()
    for i,label in enumerate(names):
        index= int(re.search(regex, label).group(1))
        #print( index)
        #print(";".join(df.iloc[index]["OC"].split(";")[2:5]))
        #labels[i]= ";".join(df.iloc[index]["OC"].split(";")[2:6])
        labels[i]= df.iloc[index]["name"]


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
    #plt.show()
    
    
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
    """subprocess.Popen([
        "firefox",
        "--new-tab",
        filepath
    ])"""


def maxpairwise_branchlength_selector(tree, numberOfTaxa, pdb_files=None, threshold=False, max_attempts=1000):
    """
    @param tree: tree taxa should be sampled from
    @type tree: Bio.Phylo.Newick.Tree
    @param numberOfTaxa: total number of taxa to sample
    @type numberOfTaxa: int
    @param threshold: minimum pairwise distance allowed, or False/0 to auto-find max threshold
    @type threshold: float or bool
    """
    leaves = tree.get_terminals()
    num_leaves = len(leaves)
    
    if numberOfTaxa > num_leaves:
        raise ValueError(f"Cannot pick {numberOfTaxa} taxa; tree only has {num_leaves} leaves.")

    all_taxa = [leaf.name for leaf in leaves]
    dist_mat = np.zeros((num_leaves, num_leaves))
    distlist = []

    # Compute upper triangle only (faster since distance is symmetric)
    for i in range(num_leaves):
        for j in range(i + 1, num_leaves):
            d = tree.distance(leaves[i], leaves[j])
            dist_mat[i, j] = d
            dist_mat[j, i] = d
            distlist.append(d)

    distlist.sort(reverse=True)

    # Specific distance threshold provided
    if threshold:
        indices = list(range(num_leaves))
        for attempt in range(max_attempts):
            random.shuffle(indices)
            selected_idx = []
            
            for ca in indices:
                if all(dist_mat[ca, cb] >= threshold for cb in selected_idx):
                    selected_idx.append(ca)
                    if len(selected_idx) == numberOfTaxa:
                        return [all_taxa[idx] for idx in selected_idx]

        raise RuntimeError(
            f"Could not find {numberOfTaxa} taxa with minimum distance {threshold} "
            f"after {max_attempts} attempts. Try lowering threshold or numberOfTaxa."
        )

    # automatically find the largest possible threshold
    else:
        for num in distlist:
            # Note: This checks degree in distance graph, not strict pairwise clique
            fil = np.sum(dist_mat >= num, axis=0) >= (numberOfTaxa - 1)
            if np.sum(fil) >= numberOfTaxa:
                taxalist = list(compress(all_taxa, fil))
                print("Taxalist result", np.sum(fil), len(taxalist))
                return taxalist[:numberOfTaxa] # Trim to exact numberOfTaxa if needed

        raise RuntimeError(f"Could not find {numberOfTaxa} valid taxa.")


def workflow(numberOfTaxa):
    referencestructure= "/data/joscha/3D-Structures/boltz_results_MOTH_1_55_TGO1_HUMAN_Reviewed__1907_AA_TANGO1_homolog/predictions/MOTH_1_55_TGO1_HUMAN_Reviewed__1907_AA_TANGO1_homolog/MOTH_1_55_TGO1_HUMAN_Reviewed__1907_AA_TANGO1_homolog_model_0.cif"
    pdb_files = "/data/joscha/3D-Structures/"
    tree_file_name = sys.argv[1]
    tree_file_OC = sys.argv[1].replace("name","OC")
    ref_tree_OC= dp.read_tree(tree_file)
    ref_tree_name = dp.read_tree(tree_file)
    outfile_prefix = '/data/joscha/output/qtools/'+ str(Path(tree_file).stem).replace(".tree","")+f"_{numberOfTaxa}Test/"
    fasta_out=outfile_prefix+Path(pdb_files).name+".fasta"
    #positions_of_interest=[x-23 for x in [54,56,63,85,98,100,101]]
    #print("Positions of interest:", positions_of_interest)

    print(fasta_out )
    Path(fasta_out).parent.mkdir(exist_ok=True)
    Path(outfile_prefix).mkdir(parents=True, exist_ok=True)


    # prepare names of species used for testing
    #train_names = sample([leaf.name for leaf in ref_tree_name.get_terminals()],35)
    train_names = maxpairwise_branchlength_selector(ref_tree_name,numberOfTaxa,pdb_files, 0.14)
    print(len(train_names))
    #_(\d+)_
    IDs= [re.search(r'_(\d+)', name).group(1) for name in train_names]

    sorted_pairs = natsorted(zip(train_names, IDs), key=lambda pair: pair[1])
    train_names, IDs = zip(*sorted_pairs)
    print(train_names, IDs)

    pattern= '|'.join(IDs)
    regex = re.compile(rf'_({pattern})_\d+_')

    matches = os_sorted([p for p in Path(pdb_files).iterdir() if regex.search(p.name)])
    #print("misMatches:",len(IDs), IDs, matches)
    #print("Number of matches:", len(matches))

    structures=[]
    names_struct={}
    regx=re.compile(r"_(\d+)_\d+_")
    z=0
    structure_paths=""
    for subdir in matches:
        for cif in subdir.rglob("*model_0.cif"):
            structure_paths+=str(cif)+" "
            identifier=  re.findall(regx,str(cif.name))[0]
            #print("File_Found:", str(cif),identifier)
            structures.append(read_cif(identifier,str(cif)))
            names_struct[identifier]=z
            break
        z+=1
    structure_paths+=referencestructure
    #print(names_struct)
    #print(list(names_struct.keys()))    


    seqrecords=extract_multi_seqrecords(list(names_struct.keys()),structures)
    #print(seqrecords, len(seqrecords))
    SeqIO.write(seqrecords, fasta_out, "fasta")
    aln_pathclustal=fasta_out.replace(".fasta","_clustal.fa")
    aln_path=fasta_out.replace(".fasta","_fdmason.fa")

    tree_path=aln_path.replace(".fa",".dnd")
    aln_path=fasta_out.replace(".fasta","_fdmason")
    os.system(f"foldmason easy-msa {structure_paths} {aln_path}  tmpFolder --report-mode 2")
    #os.system(f"clustalo -i {fasta_out} -o {aln_pathclustal} --outfmt=a2m --guidetree-out={tree_path} --force")


    alignments = AlignIO.read(aln_path+"_aa.fa", "fasta")
    #print(type(alignments))
    for alignment in alignments:
        alignment.id=re.findall(regx,str(alignment.id))[0]


    sorted_records = natsorted(alignments, key=lambda r: r.id, reverse=True)
    alignments = Align.MultipleSeqAlignment(sorted_records)
    #print(alignments, type(alignments))
    posColIdx=IndexBuilder(alignments)


    #print(train_names)
    #print(str(alignments[-1]))
    alignments=alignments[:-1]
    #print(alignments)
    columnsofinterest=[i for i in range(alignments.get_alignment_length()) if "-" not in str(alignments[:,i])]
    #print("Columnsofinterest:", columnsofinterest)
    if len(columnsofinterest)< 80:
        sorted_indices = [(i,gaps) for i, gaps in sorted(((i, c.count('-')) for i, c in enumerate(zip(*alignments))), key=lambda x: x[1]) if gaps > 0]

        current=0
        while len(columnsofinterest)< 60 and current < len(sorted_indices):
            columnsofinterest.append(sorted_indices[current][0])
            current+=1
        print(sorted_indices)
        columnsofinterest.sort()
    mtx_vectors=build_dstmtxs(columnsofinterest,alignments,structures, names_struct)
    print(len(mtx_vectors), train_names)
    print("mtx_vector:",len(mtx_vectors[0]), len(columnsofinterest))

    plot_PCA_TSNE(mtx_vectors, list(train_names), outfile_prefix)

    vectors_df=pd.DataFrame({'spec': train_names, 'mtxvector': mtx_vectors})
    vectors_df.to_csv(outfile_prefix+"vectors.csv", index=False)

    # prune the reference tree to names from testing and write to file 
    tree_name = dp.prune_tree(ref_tree_name_name, train_names)
    treefile_name=outfile_prefix+'tree_name.ph'
    dp.write_tree(tree, treefile_name)
    tree_OC = dp.prune_tree(ref_tree_name_OC, train_names)
    treefile_OC=outfile_prefix+'tree_OC.ph'
    dp.write_tree(tree_OC, treefile_OC)

    minibatches= dp.get_quartets_from_tree(tree)
    #print(minibatches)
    minibatches.to_csv(outfile_prefix  + 'minibatches.csv')

    siamesebatches = itertools.combinations(train_names, 2)
    pd.DataFrame(siamesebatches).to_csv(outfile_prefix  + 'siamesebatches.csv')

    distances = dp.get_edgelength(tree, train_names)
    distances.to_csv(outfile_prefix  + 'patristic.csv')
    positionstxt=','.join([str(n) for n in columnsofinterest])
    with open(outfile_prefix+"position.txt", "w") as text_file:
        text_file.write(positionstxt)
    os.system(f"python3 /data/joscha/qtoolsFork/qtools/src/qtoolsstart.py {treefile_name} {positionstxt}")

def main():
    workflow(45)
    
if __name__ == '__main__':
    sys.exit(main())