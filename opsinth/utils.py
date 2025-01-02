# opsin_analysis/utils.py
from opsinth import VERSION
import pysam
import logging
import sys

def read_bed_file(bed_path):
    roi = []
    with open(bed_path) as f:
        for line in f:
            l = line.strip().split()
            l[1] = max(0, int(l[1]))
            l[2] = int(l[2])
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


def aln_length_cigar(cigarstring):
    length = 0
    num = ""
    # These operations consume reference sequence
    ref_ops = set(['M', '=', 'X', 'D', 'N'])
    
    for c in cigarstring:
        if c.isnumeric():
            num += c
        else:
            if c in ref_ops:
                length += int(num)
            num = ""
            
    return length

def write_fastq(reads, outfile):
    out_prefix=outfile.replace(".fastq", "")
    # Works on aligned reads
    with open(f"{out_prefix}.fastq", "w") as fastq_file:
        for readname in reads:
            read = reads[readname]
            fastq_file.write(f"@{readname}\n{read['seq']}\n+\n{pysam.qualities_to_qualitystring(read['query_qualities'])}\n")

def write_bam(reads_aligned, reads, roi, out, version=VERSION, template_bam = False , output_format_bam = True):

    out_prefix = out.replace(".sam", "").replace(".bam", "")
    if output_format_bam:
        outfile = f"{out_prefix}.unsorted.bam"
        outfile_final = f"{out_prefix}.bam"
    else:
        outfile = f"{out_prefix}.sam"

    # Define PG Line
    pg_line = {
            'ID': 'opsinth',  # Program ID
            'PN': 'opsinth',  # Program name
            'VN': version,             # Program version
            'CL': ' '.join(sys.argv)  # Command line (you can customize this)
    }

    # Take existing header from template bam
    # Create new header for denovo reference
    if template_bam:
        b = pysam.AlignmentFile(template_bam, "rb")
        header = b.header.to_dict()
        
        header['PG'].append(pg_line)
        logging.debug("keep header from template bam and add ['PG'] line: %s", header['PG'])
    
    else:
        # Create a new header for denovo reference
        header = {
            'HD': {
                'VN': str(version),
            },
            'SQ': [
                {
                    'SN': str(roi[0][0]),
                    'LN': roi[0][2] - roi[0][1]
                }
            ],
            'PG': [pg_line]
        }
    

    if output_format_bam:
        write_format = "wb"
    else:
        write_format = "w" # Write in SAM format

    with pysam.AlignmentFile(outfile, write_format, header = header) as outf:
        for read in reads_aligned:
            
            a = pysam.AlignedSegment()

            # Recreate old coordinates for genome references
            if template_bam:
                a.reference_id  = reads[read]['reference_id']
                #a.reference_start = roi[0][1] + reads_aligned[read]['aln']['locations'][0][0]
                a.tags = reads[read]['tags'] #TODO: Add edit distance, remove obsolete tags

            # Use new coordinates for denovo reference
            else:
                a.reference_id = 0
                #a.reference_start =  reads_aligned[read]['aln']['locations'][0][1]
                a.tags = reads[read]['tags'] #TODO: Add edit distance, remove obsolete tags

            # This is always the same
            a.query_name = read
            # Reference start is 1-based (POS)
            a.reference_start = roi[0][1] + reads_aligned[read]['reference_start']
            a.query_sequence = reads_aligned[read]['seq']
            a.flag = 0 if reads_aligned[read]['strand'] == "+" else 16
            a.mapping_quality = 30 #Could be adjusted by edit distance ranges
            a.cigarstring = reads_aligned[read]['aln']['cigar']
            a.query_qualities = reads_aligned[read]['query_qualities']

            outf.write(a)
    
    if output_format_bam:
        # Sort the temporary BAM file and write to the final output BAM file
        sort_successful = False
        index_successful = False

        try:
            pysam.sort("-o", outfile_final, outfile, add_pg=True)
            sort_successful = True
        except Exception as e:
            logging.error(f"Error sorting BAM file: {e}")

        # Index the sorted BAM file
        try:
            pysam.index(outfile_final)
            index_successful = True
        except Exception as e:
            logging.error(f"Error indexing BAM file: {e}")

        # Remove the temporary BAM file only if both operations were successful
        if sort_successful and index_successful:
            import os
            os.remove(outfile)
            logging.info("Sorting and Indexing successfull, removed temp output")


def write_fasta(seq, roi, outfile):
    out_prefix=outfile.replace(".fasta", "")
    with open(f"{out_prefix}.fasta", "w") as fasta_file:
        # Fasta header coordinates are 0-based
        fasta_file.write(f">{roi[0][0]} {roi[0][1]}-{roi[0][2]}\n")
        fasta_file.write(seq + "\n")
    logging.info(f"Draft sequence written to {out_prefix}.fasta")   

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

    # Disable matplotlib font logging
    logging.getLogger('matplotlib.font_manager').disabled = True
    
    logging.debug(f"Logging level set to: {logging.getLevelName(level)}")
