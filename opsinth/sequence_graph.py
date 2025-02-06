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
                            # Use the actual alignment locations from anchor data
                            for start, end in aln['locations']:
                                read_matches[read_id].append((node_id, start, end, aln['editDistance']))
                                logging.debug(f"Phase 1: Found anchor match: read {read_id} -> {node_id} at pos {start}-{end}")

        # Phase 2: Process paraphase connections
        if paraphase_data:
            logging.info("Phase 2: Processing paraphase connections")
            for gene, gene_data in paraphase_data.items():
                if 'final_haplotypes' not in gene_data:
                    continue
                
                # Process unique supporting reads
                if 'unique_supporting_reads' in gene_data:
                    for hap_seq, hap_id in gene_data['final_haplotypes'].items():
                        supporting_reads = gene_data['unique_supporting_reads'].get(hap_seq, [])
                        for read_id in supporting_reads:
                            if read_id in reads_aligned:
                                read_seq = reads_aligned[read_id]['seq']
                                # Find actual alignment positions
                                result = edlib.align(hap_seq, read_seq, task='locations', mode='HW')
                                if result['editDistance'] / len(hap_seq) <= max_edit_distance:
                                    for start, end in result['locations']:
                                        read_matches[read_id].append((hap_id, start, end, result['editDistance']))
                                        logging.debug(f"Phase 2: Found unique match: read {read_id} -> {hap_id} at pos {start}-{end}")

                # Process non-unique supporting reads
                if 'nonunique_supporting_reads' in gene_data:
                    for hap_seq, hap_id in gene_data['final_haplotypes'].items():
                        nonunique_reads = gene_data['nonunique_supporting_reads'].get(hap_seq, [])
                        for read_id in nonunique_reads:
                            if read_id in reads_aligned:
                                read_seq = reads_aligned[read_id]['seq']
                                result = edlib.align(hap_seq, read_seq, task='locations', mode='HW')
                                if result['editDistance'] / len(hap_seq) <= max_edit_distance:
                                    for start, end in result['locations']:
                                        read_matches[read_id].append((hap_id, start, end, result['editDistance']))
                                        logging.debug(f"Phase 2: Found non-unique match: read {read_id} -> {hap_id} at pos {start}-{end}")

        # Log match statistics
        logging.info(f"Found matches for {len(read_matches)} reads")
        for read_id, matches in read_matches.items():
            logging.debug(f"Read {read_id} has {len(matches)} matches")

        # Create edges from all collected matches
        logging.info("Creating edges from collected matches")
        for read_id, matches in read_matches.items():
            if len(matches) >= 2:
                # Sort matches by position in read
                sorted_matches = sorted(matches, key=lambda x: x[1])  # x[1] is start pos
                
                # Add edges between consecutive matches that don't overlap
                for i in range(len(sorted_matches)-1):
                    for j in range(i+1, len(sorted_matches)):
                        curr_node, curr_start, curr_end, curr_dist = sorted_matches[i]
                        next_node, next_start, next_end, next_dist = sorted_matches[j]
                        
                        if curr_end <= next_start:
                            edge = (curr_node, next_node)
                            edge_support[edge].append({
                                'read_id': read_id,
                                'from_start': curr_start,
                                'from_end': curr_end,
                                'from_dist': curr_dist,
                                'to_start': next_start,
                                'to_end': next_end,
                                'to_dist': next_dist
                            })
                            logging.debug(f"Adding edge: {curr_node} -> {next_node} supported by read {read_id}")

        # Add all discovered edges to graph
        for (from_node, to_node), supports in edge_support.items():
            self.add_edge(from_node, to_node, weight=len(supports), reads=supports)
            logging.info(f"Added edge {from_node} -> {to_node} with {len(supports)} supporting reads")

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

    def read_sequences_from_fasta(self, fasta_path):
        """Read sequences from FASTA file and add them to corresponding nodes.
        
        Args:
            fasta_path (str): Path to FASTA file containing node sequences
        """
        logging.info(f"Reading sequences from {fasta_path}")
        
        current_id = None
        current_seq = []
        
        with open(fasta_path) as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # Save previous sequence if exists
                    if current_id and current_seq:
                        seq = ''.join(current_seq)
                        self._add_sequence_to_node(current_id, seq)
                    
                    # Start new sequence
                    current_id = line[1:].split()[0]  # Remove '>' and take first word
                    current_seq = []
                elif line:
                    current_seq.append(line)
        
        # Save last sequence
        if current_id and current_seq:
            seq = ''.join(current_seq)
            self._add_sequence_to_node(current_id, seq)

    def _add_sequence_to_node(self, seq_id, sequence):
        """Add reference sequence to corresponding node."""
        # Map FASTA ID to node ID
        if seq_id.startswith('anchor_'):
            # Remove double prefix if present
            node_id = seq_id.replace('anchor_anchor', 'anchor')
        elif seq_id.startswith('OPN1'):
            # Map paraphase haplotype IDs to standard gene names
            if seq_id.startswith('OPN1MW'):
                node_id = 'opn1mw_hap1'
            elif seq_id.startswith('OPN1LW'):
                node_id = 'opn1lw_hap1'
            else:
                logging.warning(f"Unknown OPN1 gene format: {seq_id}")
                return
        elif seq_id == 'LCR':
            node_id = 'LCR'
        else:
            logging.warning(f"Unknown sequence ID format: {seq_id}")
            return
        
        if node_id in self.nodes:
            self.nodes[node_id].ref_sequence = sequence
            logging.debug(f"Added reference sequence ({len(sequence)}bp) to node {node_id}")
        else:
            logging.warning(f"No matching node found for sequence {seq_id}")

def parse_paraphase_json(json_path):
    """Parse Paraphase JSON output and extract haplotype sequences"""
    with open(json_path) as f:
        data = json.load(f)
    
    nodes = {}
    logging.debug(f"Parsed json {json_path} with {len(data)} root nodes")
    
    # Process genes
    for gene in ['opn1lw']:
        if gene in data:
            gene_data = data[gene]
            nodes[gene] = []
            
            # Extract haplotype sequences and their supporting reads
            for hap_seq, hap_id in gene_data['final_haplotypes'].items():
                metadata = {
                    'supporting_reads': gene_data.get('unique_supporting_reads', {}).get(hap_seq, []),
                    'annotated_type': gene_data.get('annotated_haplotypes', {}).get(hap_id),
                }
                
                node = SequenceNode(
                    id=hap_id,
                    sequence=hap_seq,
                    source='paraphase',
                    metadata=metadata
                )
                nodes[gene].append(node)
    
    return nodes

def create_sequence_graph(results_ref, paraphase_nodes):
    """Create initial graph from anchors and Paraphase haplotypes"""
    graph = SequenceGraph()
    
    # Add anchor nodes with corrected IDs
    anchors = results_ref['anchors']
    for anchor_id, seq in anchors.items():
        # Remove potential double prefix
        node_id = f"anchor_{anchor_id}".replace('anchor_anchor', 'anchor')
        node = SequenceNode(
            id=node_id,
            sequence=seq,
            source='anchor'
        )
        graph.add_node(node)
    
    # Add Paraphase haplotype nodes with standardized names
    for gene, nodes in paraphase_nodes.items():
        for node in nodes:
            # Map haplotype IDs to standard gene names
            if node.id.startswith('opn1mw'):
                node.id = 'OPN1MW'
            elif node.id.startswith('opn1lw'):
                node.id = 'OPN1LW'
            graph.add_node(node)
    
    # Add LCR node
    lcr_node = SequenceNode(
        id='LCR',
        sequence='',  # Will be filled from FASTA file
        source='regulatory'
    )
    graph.add_node(lcr_node)
    
    return graph