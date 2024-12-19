# opsinth/analysis.py
import pysam
import edlib
from collections import Counter
from opsinth.utils import *

def read_files(bam_path, bed_path, ref_path, anchors_path):
    
    logging.info("Start reading input files")

    roi = read_bed_file(bed_path)
    bam = read_bam_file(bam_path, roi)
    seq_ref = read_reference_genome(ref_path, roi)
    anchors = read_anchors(anchors_path)
    reads = read_bam_file(bam_path, roi)
    
    logging.info("Finished reading input files")

    # Additional analysis logic can be added here
    return {
        'bam': bam,
        'roi': roi,
        'seq_ref': seq_ref,
        'anchors': anchors,
        'reads': reads
    }

def run_ref_analysis(**kwargs):
    
    logging.info("Start analysis")
    
    bam = kwargs.get('bam')
    roi = kwargs.get('roi')
    anchors = kwargs.get('anchors')
    reads = kwargs.get('reads')
    seq_ref = kwargs.get('seq_ref')
    
    # Align anchors to reference sequence
    anchors_on_ref = align_anchors_to_ref(anchors, seq_ref)

    # Find anchors on reads
    anchors_on_reads = find_anchors_on_reads(reads, anchors_on_ref)

    # Characterize reads into no anchor, unique anchor and double anchor reads
    (no_anchor_reads, unique_anchor_reads, double_anchor_reads) = characterize_read_anchors(reads, anchors_on_ref, anchors_on_reads)

    # Align reads to reference seq
    reads_aligned = align_reads_to_ref(reads, no_anchor_reads, unique_anchor_reads, double_anchor_reads, anchors_on_reads, seq_ref, anchors_on_ref)

    return{
        'seq_ref': seq_ref,
        'bam': bam,
        'roi': roi,
        'anchors': anchors,
        'anchors_on_ref': anchors_on_ref,
        'anchors_on_reads': anchors_on_reads,
        'no_anchor_reads' : no_anchor_reads,
        'unique_anchor_reads' : unique_anchor_reads,
        'double_anchor_reads' : double_anchor_reads,
        'reads_aligned': reads_aligned
    }

def run_denovo_analysis(results_ref, dataset):
    
    logging.info("Start denovo analysis")

    anchors = dataset['anchors']

    # Create draft reference sequence
    seq_draft = create_draft_ref(results_ref)

    # Align anchors to draft reference sequence
    anchors_on_ref = align_anchors_to_ref(anchors, seq_draft)

    # Find anchors on reads
    anchors_on_reads = find_anchors_on_reads(dataset['reads'], anchors_on_ref)

    # Characterize reads into no anchor, unique anchor and double anchor reads
    (no_anchor_reads, unique_anchor_reads, double_anchor_reads) = characterize_read_anchors(dataset['reads'], anchors_on_ref, anchors_on_reads)
    reads_aligned = align_reads_to_ref(dataset['reads'], no_anchor_reads, unique_anchor_reads, double_anchor_reads, anchors_on_reads, seq_draft, anchors_on_ref)

    
    # Update roi = max extend of anchors on draft seq
    min_pos = 1e8
    max_pos = 0
    for read in anchors_on_reads:
        for anchor in anchors_on_reads[read]:
            current_min = anchors_on_reads[read][anchor]['locations'][0][0]
            current_max = anchors_on_reads[read][anchor]['locations'][0][1]
            min_pos = min(min_pos, current_min)
            max_pos = max(max_pos, current_max)
    
    roi = [["denovo_ref", min_pos, max_pos]]
    
    return {
        'seq_denovo': seq_draft,
        'roi': roi,
        'anchors': anchors,
        'anchors_on_ref': anchors_on_ref,
        'anchors_on_reads': anchors_on_reads,
        'no_anchor_reads' : no_anchor_reads,
        'unique_anchor_reads' : unique_anchor_reads,
        'double_anchor_reads' : double_anchor_reads,
        'reads_aligned': reads_aligned
    }

def align_anchors_to_ref(fasta_anchors, seq_ref):
    
    """Sort anchor sequences by aligning them to the reference sequence.

    Args:
        fasta_anchors (dict): Dictionary mapping anchor names to sequences
        seq_ref (str): Reference sequence to align anchors against

    Returns:
        dict: Dictionary containing:
            - anchor_positions: Anchor alignments with start/end positions
            - forward: List of anchor names sorted by start position
            - reverse: List of anchor names sorted by end position (reversed)
    """
    anchors = {}
    for anchor in fasta_anchors:
        aln = edlib.align(fasta_anchors[anchor], seq_ref, task="path", mode="HW")
        if aln['editDistance'] > 0:
            logging.warning(f"No perfect match for {anchor}")
            logging.debug(f"Best anchor match: {aln}")

        anchors[anchor] = {
            'seq': fasta_anchors[anchor],
            'start': aln['locations'][0][0],
            'end': aln['locations'][0][1]
        }

    logging.debug(f"Processed anchors: {anchors}")
    return {
        'anchor_positions': anchors,
        'forward': [x[0] for x in sorted(anchors.items(), key=lambda x: x[1]['start'])],
        'reverse': [x[0] for x in sorted(anchors.items(), key=lambda x: x[1]['end'], reverse=True)]
    }


def find_anchors_on_reads(reads, anchors):
    """
    Align anchors to reads.

    Parameters:
    reads (dict): A dictionary containing read sequences and their metadata.
    anchors (dict): A dictionary containing anchor sequences and their positions on the reference genome.

    Returns:
    dict: A dictionary containing anchor alignments for each read.
    """
    anchor_alignments = {}

    for read in reads:
        seq = reads[read]['seq_query']
        anchor_alignments[read] = {}
        for anchor in anchors['anchor_positions']:
            aln = edlib.align(anchors['anchor_positions'][anchor]['seq'], seq, task="path", mode="HW")
            anchor_alignments[read][anchor] = aln

    logging.debug("Length of anchor table: %d", len(anchor_alignments))

    return anchor_alignments

def characterize_read_anchors(reads, anchors_on_ref, anchors_on_reads):
    """Identify and return unique anchor matches for each read while respecting read orientation.
    
    This function processes the anchor alignments for each read and ensures that the selected anchors
    are as close as possible to the start of the read. In cases where multiple anchors are found for a 
    read, the function selects the anchor with the lowest chromosome coordinate for forward reads and 
    the highest chromosome coordinate for reverse reads. Reads with anchors on both ends wirh a distance 
    greater than MAX_DISTANCE_BETWEEN_ANCHORS are considered double anchors.

    Args:
        reads (dict): A dictionary containing read sequences and their associated metadata.
        anchors (dict): A dictionary containing anchor sequences and their positions on the reference genome.
        anchor_alignments (dict): A dictionary containing the alignment results of anchors for each read.

    Returns:
        dict: 
    """
    EDIT_DISTANCE_THRESHOLD= 5
    MIN_DISTANCE_BETWEEN_ANCHORS = 40000 # At least 1 OPN1LW gene length

    valid_read_anchors = {}
    unique_anchor_reads = {}
    double_anchor_reads = {}
    no_anchor_reads = {}
    
    for read in anchors_on_reads:
        valid_read_anchors[read] = []
        for anchor in anchors_on_reads[read]:
            if anchors_on_reads[read][anchor]['editDistance'] < EDIT_DISTANCE_THRESHOLD:
                valid_read_anchors[read].append(anchor)

    for read in valid_read_anchors:
        if len(valid_read_anchors[read]) == 0:
            no_anchor_reads[read] = valid_read_anchors[read]
            logging.debug("Read %s has no valid anchors", read)
        elif len(valid_read_anchors[read]) == 1:
            unique_anchor_reads[read] = valid_read_anchors[read][0]
            logging.debug("Read %s has unique anchor: %s", read, unique_anchor_reads[read])
        elif len(valid_read_anchors[read]) > 1:
            
            distance = 0
            anchor_pair = (None, None)
            lowest_anchor_index = 1e7  # Large number
            
            # If read has multiple valid anchors, first check if largest between anchors distance > MIN_DISTANCE_BETWEEN_ANCHORS
            for anchor in valid_read_anchors[read]:
                
                # Caclulate largest pair distance
                for anchor2 in valid_read_anchors[read]:
                    new_distance =  abs(anchors_on_ref['anchor_positions'][anchor]['start'] - anchors_on_ref['anchor_positions'][anchor2]['start'])
                    if new_distance > distance:
                        distance = new_distance
                        anchor_pair = (anchor, anchor2)
            
            if distance >= MIN_DISTANCE_BETWEEN_ANCHORS:
                logging.debug("Read %s has anchors on both ends: %s", read, valid_read_anchors[read])
                double_anchor_reads[read] = anchor_pair    

            # If no pair distance > MIN_DISTANCE_BETWEEN_ANCHORS, find anchor with lowest chromosome coordinate
            else:
                for anchor in valid_read_anchors[read]:
                    # Find anchor with lowest chromosome coordinate
                    if reads[read]['strand'] == "+":
                        anchor_index = anchors_on_ref['forward'].index(anchor)
                    else :
                        anchor_index = anchors_on_ref['reverse'].index(anchor)            
                    
                    if anchor_index < lowest_anchor_index:
                        unique_anchor_reads[read] = anchor
                        logging.debug("Read %s has multiple valid anchors on one end, outmost possible anchor: %s", read, unique_anchor_reads[read])

            
    logging.info("Total reads: %d", len(reads))
    logging.info("Reads with no anchors: %d", len(no_anchor_reads))
    logging.info("Reads with single valid anchors: %d", len(unique_anchor_reads))
    logging.info("Reads with anchor on both ends: %d", len(double_anchor_reads))

    return(no_anchor_reads, unique_anchor_reads, double_anchor_reads)

def create_draft_ref(results_ref):
    # Different possible selection criteria for best spanning read
    # 2. Max aligned read length between anchors
    # 3. Highest scores of anchor hits  

    candidate_reads = {}
    double_anchor_reads = results_ref['double_anchor_reads']
    reads_aligned = results_ref['reads_aligned']
    anchors_on_reads = results_ref['anchors_on_reads']
    
    for read in double_anchor_reads:

        # Calculate anchor distances
        anchor1 = double_anchor_reads[read][0]
        anchor2 = double_anchor_reads[read][1]
        lower_anchor = anchor1 if anchors_on_reads[read][anchor1]['locations'][0][0] < anchors_on_reads[read][anchor2]['locations'][0][0] else anchor2
        upper_anchor = anchor2 if lower_anchor == anchor1 else anchor1
        pos_lower_anchor = anchors_on_reads[read][lower_anchor]['locations'][0][0]
        pos_upper_anchor = anchors_on_reads[read][upper_anchor]['locations'][0][1]
        
        candidate_reads[read] = {
            'anchor_distance': pos_upper_anchor - pos_lower_anchor,
            'lower_anchor': lower_anchor,
            'upper_anchor': upper_anchor,
            'pos_lower_anchor': pos_lower_anchor,
            'pos_upper_anchor': pos_upper_anchor,
            'edit_distance': anchors_on_reads[read][lower_anchor]['editDistance'] + anchors_on_reads[read][upper_anchor]['editDistance']
        }

    # Select best candidate read by minimum edit distance of both anchors, then by anchor distance
    best_read = min(candidate_reads, key=lambda r: (candidate_reads[r]['edit_distance'], candidate_reads[r]['anchor_distance']))

    logging.info("Best read: %s, edit distance: %d, anchor distance: %d", best_read, candidate_reads[best_read]['edit_distance'], candidate_reads[best_read]['anchor_distance'])

    seq = reads_aligned[best_read]['seq'][candidate_reads[best_read]['pos_lower_anchor']:candidate_reads[best_read]['pos_upper_anchor']]

    return seq


def align_reads_to_ref(reads, no_anchor_reads, unique_read_anchors, double_anchor_reads, anchor_alignments, seq_ref, anchors_on_ref):
    """Align reads to reference sequence starting from anchor positions.

    For each read with a unique anchor match, align the portion of the read after the anchor
    to the reference sequence starting from that anchor's position. For reverse reads, align
    the portion before the anchor to the reference sequence before the anchor position.

    Args:
        reads (dict): Dictionary containing read sequences and metadata
        no_anchor_reads (dict): Dictionary containing read sequences and metadata
        unique_read_anchors (dict): Dictionary mapping read IDs to their unique anchor match
        double_anchor_reads (dict): Dictionary mapping read IDs to their double anchor match
        anchor_alignments (dict): Dictionary containing anchor alignments for each read
        seq_ref (str): Reference sequence string
        anchors (dict): Dictionary containing anchor sequences and positions

    Returns:
        dict: Dictionary containing aligned read information including:
            - strand: Read orientation (+ or -)
            - aln: Alignment details from edlib
            - seq: Read sequence used for alignment
            - ref_length: Length of aligned reference sequence
            - query_qualities: Base quality scores for aligned sequence
    """
    reads_aligned = {}

    alignment_cases={'A':0, 'B':0, 'C':0, 'D':0, 'E':0, 'F':0, 'G':0, 'H':0}

    """
    A: Forward read, lower anchor
    B: Forward read, higher anchor
    C: Reverse read, lower anchor
    D: Reverse read, higher anchor
    E: Forward read, both anchors
    F: Reverse read, both anchors
    G: Forward read, no anchor
    H: Reverse read, no anchor
    """

    # Setup alignment function
    def align_read_to_ref(query, ref, mode):

        aln = edlib.align(query, ref, task = "path", mode = mode)
        aln_nice = edlib.getNiceAlignment(aln, query, ref)
        ref_length = len(aln_nice['target_aligned'].replace('-', ''))

        return {
            'aln' : aln,
            'seq' : query,
            'ref_length' : ref_length,
        }

    
    reads_aligned = {}

    # Align reads with no anchors (Infix Alignment)
    for read in no_anchor_reads:
        logging.debug("Aligning: %s : %s", read, no_anchor_reads[read])
        reads_aligned[read] = align_read_to_ref(reads[read]['seq_query'], seq_ref, "HW")
        
        reads_aligned[read]['query_qualities'] = reads[read]['query_qualities']
        reads_aligned[read]['strand'] = reads[read]['strand']

        if reads[read]['strand'] == "+": # + Forward read
            alignment_cases['G'] += 1
        else: # - Reverse read
            alignment_cases['H'] += 1

    # Align reads with unique anchors (PrefixAlignment)
    for read, anchor in unique_read_anchors.items():
        logging.debug("Aligning: %s : %s", read, anchor)

        is_lower_anchor = anchors_on_ref['forward'].index(anchor) < anchors_on_ref['reverse'].index(anchor)

        if reads[read]['strand'] == "+" : # + Forward read  
            if is_lower_anchor: # Case A
                seq_from_anchor = reads[read]['seq_query'][anchor_alignments[read][anchor]['locations'][0][0]:]
                qual_from_anchor = reads[read]['query_qualities'][anchor_alignments[read][anchor]['locations'][0][0]:]
                ref_from_anchor = seq_ref[anchors_on_ref['anchor_positions'][anchor]['start']:]    
                alignment_cases['A'] += 1
            else: # Case B
                seq_from_anchor = reads[read]['seq_query'][:anchor_alignments[read][anchor]['locations'][0][1]][::-1]
                qual_from_anchor = reads[read]['query_qualities'][:anchor_alignments[read][anchor]['locations'][0][1]][::-1]
                ref_from_anchor = seq_ref[:anchors_on_ref['anchor_positions'][anchor]['end']][::-1]
                alignment_cases['B'] += 1
        elif reads[read]['strand'] == "-": # - Reverse read
            if is_lower_anchor: # Case C
                seq_from_anchor = reads[read]['seq_query'][anchor_alignments[read][anchor]['locations'][0][0]:]
                qual_from_anchor = reads[read]['query_qualities'][anchor_alignments[read][anchor]['locations'][0][0]:]
                ref_from_anchor = seq_ref[anchors_on_ref['anchor_positions'][anchor]['start']]
                alignment_cases['C'] += 1
            else: # Case D
                seq_from_anchor = reads[read]['seq_query'][:anchor_alignments[read][anchor]['locations'][0][1]][::-1]
                qual_from_anchor = reads[read]['query_qualities'][:anchor_alignments[read][anchor]['locations'][0][1]][::-1]
                ref_from_anchor = seq_ref[:anchors_on_ref['anchor_positions'][anchor]['end']][::-1]
                alignment_cases['D'] += 1

        
        logging.debug(
            "Full read length: %d, Read length trimmed after anchor: %d",
            len(reads[read]['seq_query']),
            len(seq_from_anchor)
        )

        reads_aligned[read] = align_read_to_ref(seq_from_anchor, ref_from_anchor, "SHW")
        reads_aligned[read]['query_qualities'] = qual_from_anchor
        reads_aligned[read]['strand'] = reads[read]['strand']

        # Reverse high anchor reads
        if not is_lower_anchor:
            reads_aligned[read]['seq'] = seq_from_anchor[::-1]
            reads_aligned[read]['query_qualities'] = qual_from_anchor[::-1]

            high_anchor_end = anchors_on_ref['anchor_positions'][anchor]['end']
            new_aln_start = high_anchor_end - reads_aligned[read]['ref_length']
            rev_cigar = reverse_cigar(reads_aligned[read]['aln']['cigar'])

            aln_original = reads_aligned[read]['aln'].copy()
            reads_aligned[read]['aln'] = {
                'editDistance': aln_original['editDistance'],
                'alphabetLength': aln_original['alphabetLength'],
                'locations': [[new_aln_start, high_anchor_end]],
                'cigar': rev_cigar,
            }


        
    # Align reads with double anchors (Global Alignment)
    # Always start alignment from lower anchor
    
    for read, (anchor1, anchor2) in double_anchor_reads.items() :
        logging.debug("Aligning: %s : %s", read, double_anchor_reads[read])

        # Find low and high anchor on this read
        if anchor_alignments[read][anchor1]['locations'][0][0] < anchor_alignments[read][anchor2]['locations'][0][1]:
            anchor_low = anchor_alignments[read][anchor1]['locations'][0][0]
            anchor_high = anchor_alignments[read][anchor2]['locations'][0][1]
        else:
            anchor_low = anchor_alignments[read][anchor2]['locations'][0][0]
            anchor_high = anchor_alignments[read][anchor1]['locations'][0][1]

        if reads[read]['strand'] == "+" :   # + Forward read
            seq = reads[read]['seq_query'][anchor_low:anchor_high]
            qual = reads[read]['query_qualities'][anchor_low:anchor_high]
            ref = seq_ref[anchor_low:anchor_high]
            alignment_cases['E'] += 1
        else:                               # - Reverse read
            seq = reads[read]['seq_query'][anchor_low:anchor_high][::-1]
            qual = reads[read]['query_qualities'][anchor_low:anchor_high][::-1]
            ref = seq_ref[anchor_low:anchor_high][::-1]
            alignment_cases['F'] += 1

        reads_aligned[read] = align_read_to_ref(seq, ref, "NW")
        reads_aligned[read]['query_qualities'] = qual
        reads_aligned[read]['strand'] = reads[read]['strand']
        logging.debug(
            "Full read length: %d, Read length trimmed after anchor: %d",
            len(reads[read]['seq_query']),
            len(seq)
        )

    logging.debug("Alignment cases: %s", alignment_cases)   
    return reads_aligned