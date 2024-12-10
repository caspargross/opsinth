# opsin_analysis/utils.py
import pysam
import logging
import sys

def read_bed_file(bed_path, n_window=100):
    roi = []
    with open(bed_path) as f:
        for line in f:
            l = line.strip().split()
            l[1] = max(0, int(l[1]) - n_window)
            l[2] = int(l[2]) + n_window
            roi.append(l)
    return roi

def read_reference_genome(ref_path, roi):
    f = pysam.FastaFile(ref_path)
    seq_ref = f.fetch(reference=roi[0][0], start=roi[0][1], end=roi[0][2])
    f.close()
    return seq_ref

def read_anchors(fasta_anchors):
    anchors = {}
    f_anchors = pysam.FastxFile(fasta_anchors)
    for anchor in f_anchors:
        anchors[anchor.name] = anchor.sequence
    f_anchors.close()
    return anchors

def read_bam_file(bam_path, roi):
    b = pysam.AlignmentFile(bam_path, "rb")
    reads = {}
    for read in b.fetch(contig=roi[0][0], start=roi[0][1], end=roi[0][2]):
        if not read.is_supplementary and not read.is_secondary:
            reads[str(read.query_name)] = {
                'name': read.query_name,
                'reference_id': read.reference_id,
                'pos': read.get_reference_positions(full_length=True),
                'seq_query': read.query_sequence,
                'pairs': read.get_aligned_pairs(),
                'strand': "-" if read.is_reverse else "+",
                'tags': read.get_tags(),
                'query_qualities': read.query_qualities
            }
    b.close()
    return reads

def reverse_cigar(cigarstring): 
    cigar_instructions = []
    n = ""

    for c in cigarstring:
        if c.isnumeric():
            n += c
        else:
            cigar_instructions.append(n +c)
            n = ""

    return "".join(cigar_instructions[::-1])

def write_bam_file(results, reads, out_prefix, inbam, version):
    #nbam, outbam, outf, reads_aligned, anchors, roi, unique_anchor_alignments, reads
    unsortedbam = f"{out_prefix}aligned_opsin_reads.unsorted.bam"
    outbam = f"{out_prefix}aligned_opsin_reads.bam"
    reads_aligned = results['reads_aligned']
    anchors = results['anchors']
    roi = results['roi']
    unique_anchor_alignments = results['unique_anchor_alignments']

    b = pysam.AlignmentFile(inbam, "rb")
    header = b.header.to_dict()
    
    # Create a header with a @PG line
    pg_line = {
            'ID': 'opsin_analysis',  # Program ID
            'PN': 'opsin_analysis',  # Program name
            'VN': version,             # Program version
            'CL': ' '.join(sys.argv)  # Command line (you can customize this)
    }
    
    header['PG'].append(pg_line)
    logging.debug("header['PG']: %s", header['PG'])

    with pysam.AlignmentFile(unsortedbam, "wb", header = header) as outf:
        for read in reads_aligned:

            if reads_aligned[read]['strand'] == "+" :
                ref_start =  roi[0][1] + anchors['as_ref'][unique_anchor_alignments[read]]['start']
                cigarstring = reads_aligned[read]['aln']['cigar']
                query_sequence = reads_aligned[read]['seq']
            else:
                ref_start = roi[0][1] + anchors['as_ref'][unique_anchor_alignments[read]]['end'] - reads_aligned[read]['ref_length']
                cigarstring = reverse_cigar(reads_aligned[read]['aln']['cigar'])
                query_sequence = reads_aligned[read]['seq'][::-1]
            a = pysam.AlignedSegment()
            a.query_name = read
            a.query_sequence = query_sequence
            a.flag = 0 if reads_aligned[read]['strand'] == "+" else 16
            a.reference_id = reads[read]['reference_id']
            a.reference_start = ref_start
            a.mapping_quality = 30 #Could be adjusted by edit distance ranges
            a.cigarstring = cigarstring
            a.query_qualities = reads_aligned[read]['query_qualities']
            a.tags = reads[read]['tags']

            outf.write(a)
    
    # Sort the temporary BAM file and write to the final output BAM file
    sort_successful = False
    index_successful = False

    try:
        pysam.sort("-o", outbam, unsortedbam, add_pg=True)
        sort_successful = True
    except Exception as e:
        logging.error(f"Error sorting BAM file: {e}")

    # Index the sorted BAM file
    try:
        pysam.index(outbam)
        index_successful = True
    except Exception as e:
        logging.error(f"Error indexing BAM file: {e}")

    # Remove the temporary BAM file only if both operations were successful
    if sort_successful and index_successful:
        import os
        os.remove(unsortedbam)
        logging.info("Sorting and Indexing successfull, removed temp output")


    b.close()

def configure_logging(verbosity=0):
    """Configure logging based on verbosity level.
    
    Args:
        verbosity (int): Verbosity level (0=WARNING, 1=INFO, 2=DEBUG)
    """
    log_levels = {
        0: logging.WARNING,  # -v not specified (default)
        1: logging.INFO,     # -v
        2: logging.DEBUG,    # -vv
    }
    
    # Get the appropriate level or default to WARNING
    level = log_levels.get(verbosity, logging.DEBUG)
    
    # Configure logging
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    logging.debug(f"Logging level set to: {logging.getLevelName(level)}")
