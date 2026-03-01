from lxml import etree
from collections import defaultdict

treefile="/data/joscha/Data/output/TANGO1onlySP_dedup_names_aln.nexml"

def get_terminal_species_for_edges(nexml_file):
    """
    Parse a neXML file and find all terminal species names for each edge.
    
    Args:
        nexml_file: Path to the neXML file
        
    Returns:
        dict: Dictionary mapping edge IDs to lists of terminal species names
    """
    # Parse the neXML file
    tree = etree.parse(nexml_file)
    root = tree.getroot()
    
    # Get the namespace
    ns = {'nex': 'http://www.nexml.org/2009'}
    
    # Find all trees in the file
    trees = root.findall('.//nex:tree', ns)
    
    results = {}
    
    for tree_elem in trees:
        tree_id = tree_elem.get('id')
        
        # Get all nodes and edges
        nodes = {}
        for node in tree_elem.findall('.//nex:node', ns):
            node_id = node.get('id')
            label = node.get('label', '')
            otu = node.get('otu', '')
            nodes[node_id] = {'label': label, 'otu': otu, 'is_leaf': otu != ''}
        
        # Build edge structure
        edges = {}
        children = defaultdict(list)
        
        for edge in tree_elem.findall('.//nex:edge', ns):
            edge_id = edge.get('id')
            source = edge.get('source')
            target = edge.get('target')
            edges[edge_id] = {'source': source, 'target': target}
            children[source].append(target)
        
        # Get OTU labels
        otu_labels = {}
        for otu in root.findall('.//nex:otu', ns):
            otu_id = otu.get('id')
            otu_label = otu.get('label', otu_id)
            otu_labels[otu_id] = otu_label
        
        # Function to find all terminal descendants of a node
        def get_terminals(node_id):
            if nodes[node_id]['is_leaf']:
                otu_id = nodes[node_id]['otu']
                return [otu_labels.get(otu_id, nodes[node_id]['label'])]
            
            terminals = []
            for child in children[node_id]:
                terminals.extend(get_terminals(child))
            return terminals
        
        # For each edge, find terminal species
        edge_terminals = {}
        for edge_id, edge_data in edges.items():
            target_node = edge_data['target']
            terminals = get_terminals(target_node)
            edge_terminals[edge_id] = terminals
        
        results[tree_id] = edge_terminals
    
    return results

# Example usage
if __name__ == '__main__':
    nexml_file = "/data/joscha/Data/output/TANGO1onlySP_dedup_names_aln.nexml"
    
    edge_terminals = get_terminal_species_for_edges(nexml_file)
    
    # Print results
    for tree_id, edges in edge_terminals.items():
        print(f"\nTree: {tree_id}")
        for edge_id, terminals in edges.items():
            print(f"  Edge {edge_id}: {', '.join(terminals)}")