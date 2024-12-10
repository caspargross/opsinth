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
        'anchors': sort_anchors(anchors, seq_ref),
        'reads': reads
    }

def run_analysis(**kwargs):
    
    logging.info("Start analysis")
    
    bam = kwargs.get('bam')
    roi = kwargs.get('roi')
    anchors = kwargs.get('anchors')
    reads = kwargs.get('reads')
    seq_ref = kwargs.get('seq_ref')
    anchor_alignments = find_anchors_on_reads(reads, anchors)
    unique_read_anchors = find_unique_read_anchors(reads, anchors, anchor_alignments)
    reads_aligned = align_reads(reads, unique_read_anchors, anchor_alignments, seq_ref, anchors)

    return{
        'bam': bam,
        'roi': roi,
        'anchors': anchors,
        'unique_anchor_alignments' : unique_read_anchors,
        'reads_aligned': reads_aligned
    }


def sort_anchors(fasta_anchors, seq_ref):
    
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

def find_unique_read_anchors(reads, anchors, anchor_alignments):
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

    unique_read_anchors = {}
    for read in anchor_alignments:
        lowest_anchor_index = 1e5  # Large number
        for anchor in anchor_alignments[read]:
            if anchor_alignments[read][anchor]['editDistance'] < EDIT_DISTANCE_THRESHOLD:
                
                if reads[read]['strand'] == "+":
                    anchor_index = anchors['forward'].index(anchor)
                else :
                    anchor_index = anchors['reverse'].index(anchor)            
                if anchor_index < lowest_anchor_index:
                    unique_read_anchors[read] = anchor
    # For reads with two good anchors, we could align them twice, this might skew results however.
    logging.debug("Anchors per read: ", Counter(unique_read_anchors.keys()))
    logging.debug("Reads per anchor: ", Counter(unique_read_anchors.values()))
    return(unique_read_anchors)


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

    for read, anchor in unique_read_anchors.items():
        logging.debug("Aligning: %s : %s", read, anchor)

        is_lower_anchor = anchors['forward'].index(anchor) < anchors['reverse'].index(anchor)

        if reads[read]['strand'] == "+" :
            if is_lower_anchor:
                seq_from_anchor = reads[read]['seq_query'][anchor_alignments[read][anchor]['locations'][0][0]:]
                qual_from_anchor = reads[read]['query_qualities'][anchor_alignments[read][anchor]['locations'][0][0]:]
                ref_from_anchor = seq_ref[anchors['anchor_positions'][anchor]['start']:]
            else: # is_higher_anchor
                seq_from_anchor = reads[read]['seq_query'][anchor_alignments[read][anchor]['locations'][0][1]:]
                qual_from_anchor = reads[read]['query_qualities'][anchor_alignments[read][anchor]['locations'][0][1]:]
                ref_from_anchor = seq_ref[:anchors['anchor_positions'][anchor]['end']]
        elif reads[read]['strand'] == "-":
            if is_lower_anchor:
                seq_from_anchor = reads[read]['seq_query'][:anchor_alignments[read][anchor]['locations'][0][1]][::-1]
                qual_from_anchor = reads[read]['query_qualities'][:anchor_alignments[read][anchor]['locations'][0][1]][::-1]
                ref_from_anchor = seq_ref[:anchors['anchor_positions'][anchor]['end']][::-1]
            else: # is_higher_anchor
                seq_from_anchor = reads[read]['seq_query'][:anchor_alignments[read][anchor]['locations'][0][0]][::-1]
                qual_from_anchor = reads[read]['query_qualities'][:anchor_alignments[read][anchor]['locations'][0][0]][::-1]
                ref_from_anchor = seq_ref[:anchors['anchor_positions'][anchor]['start']][::-1]
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
    return reads_aligned