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
        'anchors': align_anchors(anchors, seq_ref),
        'reads': reads
    }

def run_ref_analysis(**kwargs):
    
    logging.info("Start analysis")
    
    bam = kwargs.get('bam')
    roi = kwargs.get('roi')
    anchors = kwargs.get('anchors')
    reads = kwargs.get('reads')
    seq_ref = kwargs.get('seq_ref')
    anchor_alignments = find_anchors_on_reads(reads, anchors)
    (unique_read_anchors, double_anchor_reads) = characterize_read_anchors(reads, anchors, anchor_alignments)
    reads_aligned = align_reads(reads, unique_read_anchors, anchor_alignments, seq_ref, anchors)

    return{
        'bam': bam,
        'roi': roi,
        'anchors': anchors,
        'unique_anchor_alignments' : unique_read_anchors,
        'double_anchor_reads' : double_anchor_reads,
        'reads_aligned': reads_aligned
    }

def run_denovo_analysis(double_anchor_reads, reads, anchors):
    
    logging.info("Start denovo analysis")
    
    draft_seq = create_draft_ref(double_anchor_reads, reads, anchors)
    anchor_alignments_draft = find_anchors_on_reads(reads, align_anchors(reads, draft_seq), draft_seq)

    return {}

def align_anchors(fasta_anchors, seq_ref):
    
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

def characterize_read_anchors(reads, anchors, anchor_alignments):
    """Return unique anchor matches for each read. Makes sure that read orientation is respected and
      anchors are as close as possible to the read start. When multiple anchors are found, the anchor
      with the lowest chromosome coordinate is selected.

    For reads that match multiple anchor sequences, select the correct anchor based on strand:
    - Forward read: Choose anchor with lowest coordinates 
    - Reverse read: Choose anchor with highest coordinates

    Args:
        reads (dict): Dictionary containing read sequences and metadata
        anchors (dict): Dictionary containing anchor sequences and positions
        anchor_alignments (dict): Dictionary containing anchor alignments for each read

    Returns:
        dict: Dictionary mapping read IDs to their unique anchor match
    """
    EDIT_DISTANCE_THRESHOLD= 5
    MIN_DISTANCE_BETWEEN_ANCHORS = 40000 # At least 1 OPN1LW gene length

    valid_read_anchors = {}
    unique_read_anchors = {}
    double_anchor_reads = {}
    
    for read in anchor_alignments:
        valid_read_anchors[read] = []
        for anchor in anchor_alignments[read]:
            if anchor_alignments[read][anchor]['editDistance'] < EDIT_DISTANCE_THRESHOLD:
                valid_read_anchors[read].append(anchor)

    for read in valid_read_anchors:
        # If read has multiple valid anchors, find the largest pair distance > MIN_DISTANCE_BETWEEN_ANCHORS
        if len(valid_read_anchors[read]) > 1:
            
            distance = 0
            anchor_pair = (None, None)
            lowest_anchor_index = 1e7  # Large number
            
            for anchor in valid_read_anchors[read]:
                # Find anchor with lowest chromosome coordinate
                if reads[read]['strand'] == "+":
                    anchor_index = anchors['forward'].index(anchor)
                else :
                    anchor_index = anchors['reverse'].index(anchor)            
                if anchor_index < lowest_anchor_index:
                    unique_read_anchors[read] = anchor

                # Caclulate largest pair distance
                for anchor2 in valid_read_anchors[read]:
                    new_distance =  abs(anchors['anchor_positions'][anchor]['start'] - anchors['anchor_positions'][anchor2]['start'])
                    if new_distance > distance:
                        distance = new_distance
                        anchor_pair = (anchor, anchor2)
            if distance >= MIN_DISTANCE_BETWEEN_ANCHORS:
                logging.debug("Read %s has multiple valid anchors: %s", read, valid_read_anchors[read])
                double_anchor_reads[read] = anchor_pair
        elif len(valid_read_anchors[read]) == 1:
            unique_read_anchors[read] = valid_read_anchors[read][0]
        else:
            logging.debug("Read %s has no valid anchors", read)

    logging.info("Total reads: %d", len(reads))
    logging.info("Reads with valid anchors: %d", len(unique_read_anchors))
    logging.info("Reads with double anchors: %d", len(double_anchor_reads))
    # For now unique_read_anchors also contains the double anchors. Might need changing later.
    return(unique_read_anchors, double_anchor_reads)

def create_draft_ref(double_anchor_reads, reads, anchors):
    max_qscore = 0
    for read in double_anchor_reads:
        # Calculate mean without depending on numpy
        mean_qscore = sum(reads[read]['query_qualities']) / len(reads[read]['query_qualities']) if reads[read]['query_qualities'] else 0
        if (mean_qscore > max_qscore):
            max_qscore = mean_qscore
            max_read = read
    logging.info("Max qscore: %d", max_qscore)
    logging.info("Max read: %s", max_read)  
    seq = reads[max_read]['seq_query']

    # TODO 
    # If no spannign reads are available, try to find overlap of reads

    return seq


def align_reads(reads, unique_read_anchors, anchor_alignments, seq_ref, anchors):
    """Align reads to reference sequence starting from anchor positions.

    For each read with a unique anchor match, align the portion of the read after the anchor
    to the reference sequence starting from that anchor's position. For reverse reads, align
    the portion before the anchor to the reference sequence before the anchor position.

    Args:
        reads (dict): Dictionary containing read sequences and metadata
        unique_read_anchors (dict): Dictionary mapping read IDs to their unique anchor match
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

    for read, anchor in unique_read_anchors.items():
        logging.debug("Aligning: %s : %s", read, anchor)

        is_lower_anchor = anchors['forward'].index(anchor) < anchors['reverse'].index(anchor)

        if reads[read]['strand'] == "+" :
            if is_lower_anchor: # Case A
                seq_from_anchor = reads[read]['seq_query'][anchor_alignments[read][anchor]['locations'][0][0]:]
                qual_from_anchor = reads[read]['query_qualities'][anchor_alignments[read][anchor]['locations'][0][0]:]
                ref_from_anchor = seq_ref[anchors['anchor_positions'][anchor]['start']:]    
                alignment_cases['A'] += 1
            else: # Case B
                seq_from_anchor = reads[read]['seq_query'][:anchor_alignments[read][anchor]['locations'][0][1]][::-1]
                qual_from_anchor = reads[read]['query_qualities'][:anchor_alignments[read][anchor]['locations'][0][1]][::-1]
                ref_from_anchor = seq_ref[:anchors['anchor_positions'][anchor]['end']][::-1]
                alignment_cases['B'] += 1
        elif reads[read]['strand'] == "-":
            if is_lower_anchor: # Case C
                seq_from_anchor = reads[read]['seq_query'][anchor_alignments[read][anchor]['locations'][0][0]:]
                qual_from_anchor = reads[read]['query_qualities'][anchor_alignments[read][anchor]['locations'][0][0]:]
                ref_from_anchor = seq_ref[anchors['anchor_positions'][anchor]['start']]
                alignment_cases['C'] += 1
            else: # Case D
                seq_from_anchor = reads[read]['seq_query'][:anchor_alignments[read][anchor]['locations'][0][1]][::-1]
                qual_from_anchor = reads[read]['query_qualities'][:anchor_alignments[read][anchor]['locations'][0][1]][::-1]
                ref_from_anchor = seq_ref[:anchors['anchor_positions'][anchor]['end']][::-1]
                alignment_cases['D'] += 1
        else:
            logging.error("Read %s has invalid strand: %s", read, reads[read]['strand'])
            continue
        
        logging.debug(
            "Full read length: %d, Read length trimmed after anchor: %d",
            len(reads[read]['seq_query']),
            len(seq_from_anchor)
        )

        aln = edlib.align(seq_from_anchor, ref_from_anchor, task = "path", mode = "SHW")
        aln_nice = edlib.getNiceAlignment(aln, seq_from_anchor, ref_from_anchor)
        ref_length = len(aln_nice['target_aligned'].replace('-', ''))

        reads_aligned[read] = {
            'strand' : reads[read]['strand'],
            'aln' : aln,
            'seq' : seq_from_anchor,
            'ref_length' : ref_length,
            'query_qualities' : qual_from_anchor,
        }
    
    logging.debug("Alignment cases: %s", alignment_cases)   
    return reads_aligned