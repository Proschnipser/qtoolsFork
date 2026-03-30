from ete3 import Tree
import pandas as pd
import re
import os
import sys

def get_terminal_species_for_edges(tree):
    """
    Parse a Newick file and find all terminal species for each edge.
    
    Returns a dictionary where keys are node names (or internal identifiers)
    and values are lists of terminal species names descended from that node.
    """
    # Load the tree
    
    
    edge_terminals = {}
    # Traverse all nodes in the tree
    for node in tree.traverse("postorder"):
        # Get all leaf (terminal) names descended from this node
        terminal_names = [leaf.name for leaf in node.get_leaves()]
        
        # Create a unique identifier for this edge/node,
        # You can use the node name if it exists, otherwise use id
        node_id = id(node)
        edge_terminals[node_id] = terminal_names
    
    return edge_terminals

def find_taxa(IDs, df):
    """
    Find common taxanomy in allsubnodes and returns the information as a String.
    @param IDs: IDs to find terminal nodes/species in dataframe
    @type IDs: list
    @param df: Data frame from csv with exoasy annotation
    @type df: pandas dataframe
    @return: String containing information about the taxonomy of the tree node
    @rtype: String
    """
    namelist= df["name"][IDs].tolist()
    print("Common names:",namelist)
    countdict={}
    for name in namelist:
        countdict[name] = countdict.get(name, 0) + 1
    percentages=""
    for k, v in sorted(countdict.items(), key=lambda item: item[1],  reverse=True):
        percentages += k+" "+ str(v/len(namelist)*100)+ "% "
    return percentages
    


def name_nodes_by_taxa(tree, edge_data, df):
    """
    Names nodes according to taxa present in all subnodes
    @param tree: tree object from ete3
    @type tree: ete3.Tree
    @param edge_data: Dictionary containing terminal names with node IDs as keys
    @type edge_data: dict
    @param df: Data frame from csv with exoasy annotation
    @type df: pandas dataframe
    @return: String containing information about the taxonomy of the tree node
    @rtype: String
    """
    regex = re.compile("(?<=_)\d+(?=_)")
    i=0
    for node in tree.traverse("postorder"):
        idlist=[]
        if not node.is_leaf():
            # Internal nodes get numbered names if they don't have one
            for name in edge_data[id(node)]:
                m=re.search(regex, name)
                if m:
                    idlist.append(int(m.group()))
            print(idlist)
            node.name = find_taxa(idlist, df)
            node.name="'"+str(i)+" "+node.name+"'"
            print(f"Named node: {node.name}")
            i+=1
    return tree

filepath=sys.argv[1]
df = pd.read_csv("/data/joscha/Data/uniprot_2025_01_08_MIAs_02_MOTHs.csv", index_col='no.')
newick_file = filepath #"/data/joscha/Data/output/TANGO1onlySP_dedup_names_aln.tree"
tree = Tree(newick_file, format=1)
edge_data = get_terminal_species_for_edges(tree)
tree= name_nodes_by_taxa(tree,edge_data, df)
# for node in tree.traverse("postorder"):
#     print(node)
output= filepath.replace("con.tree", "_annotatedname.tree").replace(".con.nwk", "_annotatedname.tree")
print(output)
tree.write(outfile=output, format=1)
# Print results
# for edge, terminals in edge_data.items():
#     print(f"Edge/Node: {edge}")
#     print(f"  Terminal species: {', '.join(terminals)}")
#     print()


