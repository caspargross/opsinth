# opsin_analysis/analysis.py
import pysam
import edlib
from collections import Counter
from opsin_analysis.utils import *

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
    anchor_alignments = align_anchors(reads, anchors)
    unique_read_anchors = find_read_anchors(reads, anchors, anchor_alignments)
    reads_aligned = align_reads(reads, unique_read_anchors, anchor_alignments, seq_ref, anchors)

    return{
        'bam': bam,
        'roi': roi,
        'anchors': anchors,
        'unique_anchor_alignments' : unique_read_anchors,
        'reads_aligned': reads_aligned
    }


def sort_anchors(fasta_anchors, seq_ref):
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
        'as_ref': anchors,
        'forward': [x[0] for x in sorted(anchors.items(), key=lambda x: x[1]['start'])],
        'reverse': [x[0] for x in sorted(anchors.items(), key=lambda x: x[1]['end'], reverse=True)]
    }


def align_anchors(reads, anchors):
    anchor_alignments = {}

    for read in reads:
        seq = reads[read]['seq_query']
        anchor_alignments[read] = {}
        for anchor in anchors['as_ref']:
            aln = edlib.align(anchors['as_ref'][anchor]['seq'], seq, task="path", mode="HW")
            anchor_alignments[read][anchor] = aln

    logging.debug("Length of anchor table: %d", len(anchor_alignments))

    # Reads that match more then one anchor sequence span multiple anchors. To get the correct start position we need to chose the correct anchor sequence as follows: 
    #  - Take  the anchor with _lowest_ coordinates when on the forward strand
    #  - Select the anchor with _highest_ coordinates when on the reverse strand

    # In order to get correct coordinates: 
    #  1) First align anchors to ref genome
    #  2) Extract coordinates
    #  3) Sort anchors matches ascending start
    #  4) Sort anchors matches by descending end
    return anchor_alignments

def find_read_anchors(reads, anchors, anchor_alignments):
    DISTANCE_THRESHOLD= 5

    unique_read_anchors = {}
    for read in anchor_alignments:
        lowest_anchor_index = 1e5  # Large number
        for anchor in anchor_alignments[read]:
            if anchor_alignments[read][anchor]['editDistance'] < DISTANCE_THRESHOLD:
                
                if reads[read]['strand'] == "+":
                    anchor_index = anchors['forward'].index(anchor)
                else :
                    anchor_index = anchors['reverse'].index(anchor)            
                if anchor_index < lowest_anchor_index:
                    unique_read_anchors[read] = anchor

    logging.debug("Anchors per read: ", Counter(unique_read_anchors.keys()))
    logging.debug("Reads per anchor: ", Counter(unique_read_anchors.values()))
    return(unique_read_anchors)


def align_reads(reads, unique_read_anchors, anchor_alignments, seq_ref, anchors):
    
    reads_aligned = {}

    for read, anchor in unique_read_anchors.items():
        logging.debug("Aligning: %s : %s", read, anchor)

        if reads[read]['strand'] == "+":
            seq_from_anchor = reads[read]['seq_query'][anchor_alignments[read][anchor]['locations'][0][0]:]
            qual_from_anchor = reads[read]['query_qualities'][anchor_alignments[read][anchor]['locations'][0][0]:]
            ref_from_anchor = seq_ref[anchors['as_ref'][anchor]['start']:]
        else:
            seq_from_anchor = reads[read]['seq_query'][:anchor_alignments[read][anchor]['locations'][0][1]][::-1]
            qual_from_anchor = reads[read]['query_qualities'][:anchor_alignments[read][anchor]['locations'][0][1]][::-1]
            ref_from_anchor = seq_ref[:anchors['as_ref'][anchor]['end']][::-1]
        
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