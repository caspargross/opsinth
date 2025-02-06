import networkx as nx
import matplotlib.pyplot as plt
from dataclasses import dataclass
import json
import logging
from collections import defaultdict
import edlib

@dataclass
class SequenceNode:
    id: str
    sequence: str
    source: str  # 'anchor', 'paraphase', 'novel'
    metadata: dict = None
    # Reference mapping information
    ref_sequence: str = None
    ref_start: int = None
    ref_end: int = None
    edit_distance: int = None
    identity: float = None
    cigar: str = None

class SequenceGraph:
    def __init__(self):
        self.graph = nx.DiGraph()
        self.nodes = {}
        
    def add_node(self, node):
        if node.id not in self.nodes:
            self.nodes[node.id] = node
            self.graph.add_node(node.id, sequence=node.sequence, 
                              source=node.source, metadata=node.metadata)
            
    def add_edge(self, from_node, to_node, weight=1, reads=None):
        if reads is None:
            reads = []
        self.graph.add_edge(from_node, to_node, weight=weight, reads=reads)
    
    def find_connections(self, reads_aligned, max_edit_distance=0.1):
        """Find connections between nodes using aligned reads"""
        edge_support = defaultdict(list)
        
        for read_id, read_data in reads_aligned.items():
            matching_nodes = self._find_matching_nodes(read_data)
            
            if len(matching_nodes) >= 2:
                # Sort nodes by position in read
                sorted_nodes = sorted(matching_nodes, 
                                   key=lambda x: x[1])  # x[1] is start pos
                
                # Add edges between consecutive nodes
                for i in range(len(sorted_nodes)-1):
                    curr_node, _, _ = sorted_nodes[i]
                    next_node, _, _ = sorted_nodes[i+1]
                    edge = (curr_node, next_node)
                    edge_support[edge].append(read_id)
        
        # Add edges to graph
        for (from_node, to_node), reads in edge_support.items():
            self.add_edge(from_node, to_node, 
                         weight=len(reads), reads=reads)
    
    def _find_matching_nodes(self, read_data):
        """Find which nodes match parts of this read"""
        matches = []
        read_seq = read_data['seq']
        
        for node_id, node in self.nodes.items():
            # Try to align node sequence to read
            result = edlib.align(node.sequence, read_seq, 
                               task='locations',
                               mode='HW')
            
            if result['editDistance'] / len(node.sequence) <= 0.1:
                for start, end in result['locations']:
                    matches.append((node_id, start, end))
        
        return matches
    
    def discover_novel_sequences(self, reads_aligned, min_support=2):
        """Find novel sequences between known nodes"""
        novel_seqs = defaultdict(list)
        
        for read_id, read_data in reads_aligned.items():
            matching_nodes = self._find_matching_nodes(read_data)
            
            if len(matching_nodes) >= 2:
                sorted_nodes = sorted(matching_nodes, 
                                   key=lambda x: x[1])
                
                # Look for novel sequences between matched nodes
                for i in range(len(sorted_nodes)-1):
                    curr_node, _, curr_end = sorted_nodes[i]
                    next_node, next_start, _ = sorted_nodes[i+1]
                    
                    if next_start - curr_end > 20:  # Min gap size
                        novel_seq = read_data['seq'][curr_end:next_start]
                        key = (curr_node, next_node)
                        novel_seqs[key].append(
                            {'sequence': novel_seq, 'read': read_id})
        
        # Add well-supported novel sequences as nodes
        for (prev_node, next_node), seqs in novel_seqs.items():
            if len(seqs) >= min_support:
                # Create consensus from supporting sequences
                consensus_seq = self._make_consensus([s['sequence'] 
                                                   for s in seqs])
                
                node_id = f"novel_{prev_node}_{next_node}"
                novel_node = SequenceNode(
                    id=node_id,
                    sequence=consensus_seq,
                    source='novel',
                    metadata={'supporting_reads': [s['read'] for s in seqs]}
                )
                
                self.add_node(novel_node)
                # Add edges to and from the novel node
                self.add_edge(prev_node, node_id, 
                            weight=len(seqs), 
                            reads=[s['read'] for s in seqs])
                self.add_edge(node_id, next_node, 
                            weight=len(seqs),
                            reads=[s['read'] for s in seqs])
    
    def _make_consensus(self, sequences):
        """Simple majority consensus from multiple sequences"""
        # For now just return first sequence
        # TODO: implement proper consensus calling
        return sequences[0]
    
    def find_optimal_path(self, start_nodes=None, end_nodes=None):
        """Find highest-weight path through the graph"""
        if not start_nodes:
            start_nodes = [n for n in self.graph.nodes() 
                         if self.graph.in_degree(n) == 0]
        if not end_nodes:
            end_nodes = [n for n in self.graph.nodes() 
                        if self.graph.out_degree(n) == 0]
            
        best_path = None
        best_weight = -1
        
        for start in start_nodes:
            for end in end_nodes:
                try:
                    path = nx.shortest_path(
                        self.graph, start, end,
                        weight=lambda u, v, d: -d['weight']
                    )
                    weight = sum(self.graph[path[i]][path[i+1]]['weight']
                               for i in range(len(path)-1))
                    
                    if weight > best_weight:
                        best_path = path
                        best_weight = weight
                except nx.NetworkXNoPath:
                    continue
                    
        return best_path, best_weight
    
    def plot(self, output_path=None):
        pos = nx.spring_layout(self.graph)
        plt.figure(figsize=(12, 8))
        
        # Draw nodes with different colors based on source
        colors = {'anchor': 'lightblue', 'paraphase': 'lightgreen', 
                 'novel': 'lightgray'}
        for source in colors:
            nodes = [n for n, d in self.graph.nodes(data=True) 
                    if d['source'] == source]
            nx.draw_networkx_nodes(self.graph, pos, nodelist=nodes, 
                                 node_color=colors[source])
            
        # Draw edges with width proportional to weight
        edges = self.graph.edges(data=True)
        weights = [d['weight'] for _, _, d in edges]
        nx.draw_networkx_edges(self.graph, pos, width=weights)
        
        # Add labels
        nx.draw_networkx_labels(self.graph, pos)
        
        if output_path:
            plt.savefig(output_path)
        plt.close()

    def find_node_homology(self, ref_seq, min_identity=0.8):
        """
        Find homologous sequences in reference for each node.
        Updates nodes with reference mapping information.
        """
        logging.info("Finding homologous sequences for nodes")
        
        for node_id, node in self.nodes.items():
            if len(node.sequence) < 20:
                continue
                
            result = edlib.align(
                node.sequence,
                ref_seq,
                task='locations',
                mode='HW',
                additionalEqualities=[
                    ('N', 'A'), ('N', 'C'), ('N', 'G'), ('N', 'T')
                ]
            )
            
            identity = 1 - (result['editDistance'] / len(node.sequence))
            
            if identity >= min_identity and result['locations']:
                start, end = result['locations'][0]
                
                # Update node with reference information
                node.ref_start = start
                node.ref_end = end
                node.edit_distance = result['editDistance']
                node.identity = identity
                node.cigar = result['cigar']
                node.ref_sequence = ref_seq[start:end+1]
                
                logging.debug(
                    "Node %s mapped to ref pos %d-%d (identity: %.2f)",
                    node_id, start, end, identity
                )
            else:
                logging.debug(
                    "No significant homology found for node %s",
                    node_id
                )
    
    def get_reference_structure(self):
        """Return ordered list of nodes based on their reference positions."""
        mapped_nodes = []
        for node_id, node in self.nodes.items():
            if node.ref_start is not None:
                mapped_nodes.append({
                    'node_id': node_id,
                    'start': node.ref_start,
                    'end': node.ref_end,
                    'identity': node.identity
                })
        
        return sorted(mapped_nodes, key=lambda x: x['start'])

def parse_paraphase_json(json_path):
    """Parse Paraphase JSON output and extract haplotype sequences"""
    with open(json_path) as f:
        data = json.load(f)
    
    nodes = {}
    
    # Process OPN1LW haplotypes
    if 'opn1lw' in data:
        opn1lw = data['opn1lw']
        nodes['opn1lw'] = []
        
        # Extract haplotype sequences and their supporting reads
        for hap_seq, hap_id in opn1lw['final_haplotypes'].items():
            metadata = {
                'supporting_reads': opn1lw['unique_supporting_reads'].get(hap_seq, []),
                'annotated_type': opn1lw['annotated_haplotypes'].get(hap_id),
            }
            
            node = SequenceNode(
                id=hap_id,
                sequence=hap_seq,
                source='paraphase',
                metadata=metadata
            )
            nodes['opn1lw'].append(node)
            
def create_sequence_graph(results_ref, paraphase_nodes):
    """Create initial graph from anchors and Paraphase haplotypes"""
    graph = SequenceGraph()
    
    # Add anchor nodes
    anchors = results_ref['anchors']
    for anchor_id, seq in anchors.items():
        node = SequenceNode(
            id=f"anchor_{anchor_id}",
            sequence=seq,
            source='anchor'
            )
        graph.add_node(node)
    
    # Add Paraphase haplotype nodes
    for gene, nodes in paraphase_nodes.items():
        for node in nodes:
            graph.add_node(node)
            
    return graph