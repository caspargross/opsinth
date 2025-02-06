import networkx as nx
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import Dict, List, Set, Tuple
import json
import logging

@dataclass
class SequenceNode:
    id: str
    sequence: str
    ref_sequence: str
    source: str  # 'anchor', 'paraphase', 'novel'
    metadata: dict = None

class SequenceGraph:
    def __init__(self):
        self.graph = nx.DiGraph()
        self.nodes: Dict[str, SequenceNode] = {}
        
    def add_node(self, node: SequenceNode):
        if node.id not in self.nodes:
            self.nodes[node.id] = node
            self.graph.add_node(node.id, sequence=node.sequence, 
                              source=node.source, metadata=node.metadata)
            
    def add_edge(self, from_node: str, to_node: str, 
                 weight: int = 1, reads: List[str] = None):
        if reads is None:
            reads = []
        self.graph.add_edge(from_node, to_node, weight=weight, 
                           reads=reads)
        
    def plot(self, output_path: str = None):
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

def parse_paraphase_json(json_path: str) -> Dict[str, List[SequenceNode]]:
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
                'is_first_copy': hap_id in opn1lw.get('first_copies', []),
                'is_last_copy': hap_id in opn1lw.get('last_copies', [])
            }
            
            node = SequenceNode(
                id=hap_id,
                sequence=hap_seq,
                source='paraphase',
                metadata=metadata
            )
            nodes['opn1lw'].append(node)
            
def create_sequence_graph(results_ref: dict, paraphase_nodes: Dict[str, List[SequenceNode]]):
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