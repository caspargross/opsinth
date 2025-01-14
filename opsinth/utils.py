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

def read_fasta_seq(fasta_path: str) -> str:
    """
    Read sequence from a FASTA file (supports gzipped files).
    
    Args:
        fasta_path: Path to FASTA file (can be .gz)
        
    Returns:
        String containing the sequence (without header)
        
    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If file is empty or malformed
    """
    import gzip
    
    try:
        # Check if file is gzipped
        is_gzipped = fasta_path.endswith('.gz')
        opener = gzip.open if is_gzipped else open
        
        with opener(fasta_path, 'rt') as f:
            sequence = ''
            header_found = False
            
            for line in f:
                line = line.strip()
                if not line:  # Skip empty lines
                    continue
                    
                if line.startswith('>'):
                    if not header_found:
                        header_found = True
                        continue
                    else:
                        break  # Stop at second header if multiple sequences
                else:
                    if not header_found:
                        raise ValueError("FASTA file must start with '>'")
                    sequence += line
                    
            if not sequence:
                raise ValueError("No sequence found in FASTA file")
                
            return sequence
            
    except FileNotFoundError:
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")
    except Exception as e:
        raise ValueError(f"Error reading FASTA file {fasta_path}: {str(e)}")

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
            seq = read['seq']
            quals = read['query_qualities']
            
            # Reverse complement sequence and reverse qualities for reverse strand reads
            if read['strand'] == '-':
                seq = seq[::-1].translate(str.maketrans('ACGT', 'TGCA'))
                quals = quals[::-1]
                
            fastq_file.write(f"@{readname}\n{seq}\n+\n{pysam.qualities_to_qualitystring(quals)}\n")

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
    logging.info(f"Export sequence to {out_prefix}.fasta")   

def configure_logging(verbosity: int = 0):
    """
    Configure logging based on verbosity level.
    
    Args:
        verbosity: Number of -v flags (0 = WARNING, 1 = INFO, 2 = DEBUG)
    """
    # Force reconfiguration of logging
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    
    # Set log level based on verbosity
    if verbosity == 0:
        level = logging.WARNING
    elif verbosity == 1:
        level = logging.INFO
    else:
        level = logging.DEBUG
        
    # Configure logging with more detailed format
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        force=True  # Force reconfiguration
    )
    
    # Set levels for specific loggers
    logging.getLogger('matplotlib').setLevel(logging.WARNING)
    logging.getLogger('PIL').setLevel(logging.WARNING)
    
def convert_coordinate(pos: int, cigar: str, direction: str = "query_to_ref", 
                      query_start: int = 0, ref_start: int = 0) -> int:
    """
    Convert coordinates between query and reference positions using a CIGAR string.
    Handles spliced alignments (N in CIGAR).
    """
    if direction not in ["query_to_ref", "ref_to_query"]:
        raise ValueError("direction must be either 'query_to_ref' or 'ref_to_query'")
        
    # Check if position is before start
    start_pos = query_start if direction == "query_to_ref" else ref_start
    if pos < start_pos:
        raise ValueError(f"Position {pos} is before start position {start_pos}")
    
    # Parse CIGAR string
    import re
    cigar_ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    
    curr_query = query_start
    curr_ref = ref_start
    
    for length, op in cigar_ops:
        length = int(length)
        
        if op in 'M=X':  # Match or mismatch
            if direction == "query_to_ref":
                if pos <= curr_query + length:
                    return curr_ref + (pos - curr_query)
            else:  # ref_to_query
                if pos <= curr_ref + length:
                    return curr_query + (pos - curr_ref)
            curr_query += length
            curr_ref += length
            
        elif op == 'I':  # Insertion to reference
            if direction == "query_to_ref":
                if pos <= curr_query + length:
                    return curr_ref
            curr_query += length
            
        elif op == 'D':  # Deletion from reference
            if direction == "ref_to_query":
                if pos <= curr_ref + length:
                    return curr_query
            curr_ref += length
            
        elif op == 'N':  # Skipped region/intron
            curr_ref += length  # Only advance reference position
            
        elif op == 'S':  # Soft clipping
            if direction == "query_to_ref":
                if pos <= curr_query + length:
                    return curr_ref
            curr_query += length
            
    raise ValueError(f"Position {pos} is beyond the aligned region")

# Wrapper functions for backward compatibility
def query_to_ref_pos(query_pos: int, cigar: str, query_start: int = 0, ref_start: int = 0) -> int:
    """Convert query position to reference position"""
    return convert_coordinate(query_pos, cigar, "query_to_ref", query_start, ref_start)

def ref_to_query_pos(ref_pos: int, cigar: str, query_start: int = 0, ref_start: int = 0) -> int:
    """Convert reference position to query position"""
    return convert_coordinate(ref_pos, cigar, "ref_to_query", query_start, ref_start)

def create_igv_session(reference_name, bam_name, roi, output_prefix, template_path):
    """
    Create an IGV session file from template.
    
    Args:
        reference_path (str): Path to the reference FASTA file
        bam_path (str): Path to the BAM file
        output_path (str): Where to save the IGV session file
        template_path (str): Path to the IGV session template
    """
    import os
    output_path = f"{output_prefix}.igv_session.xml"

    with open(template_path) as f:
        template = f.read()
    
    # Replace placeholders with actual paths
    session = template.replace('{{reference_name}}', os.path.basename(reference_name))
    session = session.replace('{{roi}}', f"{roi[0][0]}:{roi[0][1]}-{roi[0][2]}")
    session = session.replace('{{bam_name}}', os.path.basename(bam_name))
    
    # Write the session file
    with open(output_path, 'w') as f:
        f.write(session)

    logging.info(f"IGV session file created at {output_path}")
