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
        self.edges = {}
        
    def add_node(self, node):
        if node.id not in self.nodes:
            self.nodes[node.id] = node
            self.graph.add_node(node.id, sequence=node.sequence, 
                              source=node.source, metadata=node.metadata)
            
    def add_edge(self, from_node, to_node, weight=1, reads=None):
        if reads is None:
            reads = []
        self.graph.add_edge(from_node, to_node, weight=weight, reads=reads)
    
    def find_connections(self, reads_aligned, anchors_on_reads=None, paraphase_data=None, max_edit_distance=0.1):
        """Find connections between nodes using aligned reads in three phases."""
        read_matches = defaultdict(list)
        edge_support = defaultdict(list)
        
        # Phase 1: Process pre-computed anchor alignments
        if anchors_on_reads:
            logging.info("Phase 1: Processing pre-computed anchor alignments")
            for read_id, anchor_data in anchors_on_reads.items():
                for anchor_id, aln in anchor_data.items():
                    if aln['editDistance'] / len(aln['query']) <= max_edit_distance:
                        node_id = f"anchor_{anchor_id}".replace('anchor_anchor', 'anchor')
                        if node_id in self.nodes:
                            for start, end in aln['locations']:
                                read_matches[read_id].append((node_id, start, end, aln['editDistance']))
                                logging.debug(f"Phase 1: Found anchor match: read {read_id} -> {node_id} at pos {start}-{end}")

        # Phase 2: Process paraphase connections
        if paraphase_data:
            logging.info("Phase 2: Processing paraphase connections")
            logging.debug(f"Paraphase data keys: {list(paraphase_data.keys())}")
            
            # Process non-unique supporting reads
            if 'nonunique_supporting_reads' in paraphase_data:
                logging.debug(f"Processing non-unique reads: {len(paraphase_data['nonunique_supporting_reads'])} reads")
                for read_id, haplotypes in paraphase_data['nonunique_supporting_reads'].items():
                    if read_id in reads_aligned:
                        read_seq = reads_aligned[read_id]['seq']
                        for hap_id in haplotypes:
                            if hap_id in self.nodes:
                                result = edlib.align(self.nodes[hap_id].ref_sequence, read_seq,
                                                  task='locations',
                                                  mode='HW')
                                
                                if result['editDistance'] / len(self.nodes[hap_id].ref_sequence) <= max_edit_distance:
                                    for start, end in result['locations']:
                                        read_matches[read_id].append((hap_id, start, end, result['editDistance']))
                                        logging.debug(f"Phase 2: Found non-unique match: read {read_id} -> {hap_id} at pos {start}-{end}")
            
            # Process unique supporting reads
            if 'unique_supporting_reads' in paraphase_data:
                logging.debug(f"Processing unique reads: {len(paraphase_data['unique_supporting_reads'])} reads")
                for read_id, hap_id in paraphase_data['unique_supporting_reads'].items():
                    if read_id in reads_aligned and hap_id in self.nodes:
                        read_seq = reads_aligned[read_id]['seq']
                        result = edlib.align(self.nodes[hap_id].ref_sequence, read_seq,
                                          task='locations',
                                          mode='HW')
                        
                        if result['editDistance'] / len(self.nodes[hap_id].ref_sequence) <= max_edit_distance:
                            for start, end in result['locations']:
                                read_matches[read_id].append((hap_id, start, end, result['editDistance']))
                                logging.debug(f"Phase 2: Found unique match: read {read_id} -> {hap_id} at pos {start}-{end}")

        # Log match statistics
        logging.info(f"Found matches for {len(read_matches)} reads")
        for read_id, matches in read_matches.items():
            logging.debug(f"Read {read_id} has {len(matches)} matches: {matches}")

        # Create edges from all collected matches
        logging.info("Creating edges from collected matches")
        for read_id, matches in read_matches.items():
            if len(matches) >= 2:
                # Sort matches by position in read
                sorted_matches = sorted(matches, key=lambda x: x[1])  # x[1] is start pos
                logging.debug(f"Sorted matches for read {read_id}: {sorted_matches}")
                
                # Add edges between consecutive non-overlapping matches
                for i in range(len(sorted_matches)-1):
                    curr_node, curr_start, curr_end, curr_dist = sorted_matches[i]
                    next_node, next_start, next_end, next_dist = sorted_matches[i+1]
                    
                    # Only create edge if positions don't overlap and nodes are different
                    if curr_end <= next_start and curr_node != next_node:
                        edge = (curr_node, next_node)
                        support_info = {
                            'read_id': read_id,
                            'from_start': curr_start,
                            'from_end': curr_end,
                            'from_dist': curr_dist,
                            'to_start': next_start,
                            'to_end': next_end,
                            'to_dist': next_dist
                        }
                        edge_support[edge].append(support_info)
                        logging.debug(f"Adding edge: {curr_node} -> {next_node} supported by read {read_id} "
                                    f"(positions: {curr_start}-{curr_end} -> {next_start}-{next_end})")

        # Add all discovered edges to graph
        edge_count = 0
        for (from_node, to_node), supports in edge_support.items():
            self.add_edge(from_node, to_node, weight=len(supports), reads=supports)
            edge_count += 1
            logging.info(f"Added edge {from_node} -> {to_node} with {len(supports)} supporting reads")
        
        logging.info(f"Created {edge_count} edges in total")

    def _find_node_matches_in_read(self, read_seq, node_id, node):
        """Find all occurrences of a node's sequence in a read.
        
        Args:
            read_seq (str): Read sequence
            node_id (str): Node identifier
            node (SequenceNode): Node object containing sequence
        
        Returns:
            list: List of tuples (node_id, start, end, edit_distance)
        """
        matches = []
        
        # Use reference sequence if available, otherwise use node sequence
        query_seq = node.ref_sequence if node.ref_sequence else node.sequence
        
        # Find all occurrences using edlib
        result = edlib.align(query_seq, read_seq,
                            task='locations',
                            mode='HW')
        
        if result['editDistance'] / len(query_seq) <= 0.1:  # 90% identity threshold
            for start, end in result['locations']:
                matches.append((node_id, start, end, result['editDistance']))
                logging.debug(f"Found match for {node_id} at {start}-{end} (edit distance: {result['editDistance']})")
        
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
        colors = {
            'anchor': 'lightblue', 
            'paraphase': 'lightgreen',
            'regulatory': 'lightpink',
            'novel': 'lightgray'
        }
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
        else:
            plt.show()
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

class SequenceManager:
    def __init__(self):
        self.sequences = {}
        
    def add_fasta_file(self, fasta_path):
        """Read sequences from FASTA file."""
        logging.info(f"Reading sequences from FASTA file: {fasta_path}")
        current_id = None
        current_seq = []
        
        with open(fasta_path) as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_id and current_seq:
                        self.sequences[current_id] = ''.join(current_seq)
                        logging.debug(f"Read FASTA sequence: {current_id} ({len(self.sequences[current_id])}bp)")
                    current_id = line[1:].split()[0]  # Remove '>' and take first word
                    logging.debug(f"Found FASTA header: {line} -> parsed ID: {current_id}")
                    current_seq = []
                elif line:
                    current_seq.append(line)
        
        if current_id and current_seq:
            self.sequences[current_id] = ''.join(current_seq)
            logging.debug(f"Read final FASTA sequence: {current_id} ({len(self.sequences[current_id])}bp)")
        
        logging.info(f"Added {len(self.sequences)} sequences from FASTA file")
        logging.debug(f"Available sequence IDs: {list(self.sequences.keys())}")
    
    def add_sequences(self, sequences):
        """Add sequences from dictionary."""
        logging.info(f"Adding {len(sequences)} sequences from dictionary")
        self.sequences.update(sequences)
        logging.debug(f"Available sequence IDs: {list(self.sequences.keys())}")
    
    def get_sequence(self, node_id):
        """Get sequence for a node ID."""
        if node_id.startswith('opn1mw'):
            fasta_id = 'OPN1MW'
        elif node_id.startswith('opn1lw'):
            fasta_id = 'OPN1LW'
        elif node_id.startswith('anchor_'):
            # Keep the full anchor ID for FASTA lookup
            fasta_id = node_id
        elif node_id == 'LCR':
            fasta_id = 'LCR'
        elif node_id == 'OPN1MW' or node_id == 'OPN1LW':
            # Handle standardized gene names
            fasta_id = node_id
        else:
            logging.warning(f"Unknown node ID format: {node_id}")
            return None
            
        if fasta_id in self.sequences:
            logging.debug(f"Found sequence for node {node_id} (FASTA ID: {fasta_id}, length: {len(self.sequences[fasta_id])}bp)")
            return self.sequences[fasta_id]
        else:
            logging.warning(f"No sequence found for node {node_id} (FASTA ID: {fasta_id})")
            return None

def parse_paraphase_nodes(paraphase_json):
    """Parse Paraphase output JSON into nodes.
    
    Args:
        paraphase_json (dict): Paraphase output JSON data
        
    Returns:
        dict: Dictionary mapping gene names to lists of nodes
        dict: Supporting reads data structure
    """
    nodes = defaultdict(list)
    supporting_reads = {
        'unique_supporting_reads': {},
        'nonunique_supporting_reads': {}
    }
    
    for gene, gene_data in paraphase_json.items():
        if 'haplotypes' not in gene_data:
            continue
            
        logging.debug(f"Processing gene {gene}")
        
        # Process haplotypes
        for hap_id, hap_data in gene_data['haplotypes'].items():
            sequence = hap_data.get('sequence', '')
            node = SequenceNode(
                id=hap_id,
                sequence=sequence,
                source='paraphase',
                metadata={'gene': gene}
            )
            nodes[gene].append(node)
            
        # Process unique supporting reads
        if 'unique_supporting_reads' in gene_data:
            for read_id, hap_id in gene_data['unique_supporting_reads'].items():
                supporting_reads['unique_supporting_reads'][read_id] = hap_id
                logging.debug(f"Added unique supporting read {read_id} -> {hap_id}")
                
        # Process non-unique supporting reads
        if 'nonunique_supporting_reads' in gene_data:
            for read_id, haplotypes in gene_data['nonunique_supporting_reads'].items():
                supporting_reads['nonunique_supporting_reads'][read_id] = haplotypes
                logging.debug(f"Added non-unique supporting read {read_id} -> {haplotypes}")
    
    logging.info(f"Parsed {sum(len(n) for n in nodes.values())} nodes from {len(nodes)} genes")
    logging.info(f"Found {len(supporting_reads['unique_supporting_reads'])} unique and "
                f"{len(supporting_reads['nonunique_supporting_reads'])} non-unique supporting reads")
    
    return nodes, supporting_reads

def create_sequence_graph(results_ref, paraphase_nodes, sequence_manager):
    """Create initial graph from anchors and Paraphase haplotypes"""
    graph = SequenceGraph()
    
    # Add anchor nodes with corrected IDs
    anchors = results_ref['anchors']
    for anchor_id, seq in anchors.items():
        # Remove potential double prefix
        node_id = f"anchor_{anchor_id}".replace('anchor_anchor', 'anchor')
        sequence = sequence_manager.get_sequence(node_id)
        node = SequenceNode(
            id=node_id,
            sequence=seq,
            source='anchor',
            ref_sequence=sequence
        )
        graph.add_node(node)
    
    # Add Paraphase haplotype nodes
    for gene, nodes in paraphase_nodes.items():
        for node in nodes:
            sequence = sequence_manager.get_sequence(node.id)
            node.ref_sequence = sequence
            graph.add_node(node)
    
    # Add LCR node
    sequence = sequence_manager.get_sequence('LCR')
    lcr_node = SequenceNode(
        id='LCR',
        sequence='',
        source='regulatory',
        ref_sequence=sequence
    )
    graph.add_node(lcr_node)
    
    return graph