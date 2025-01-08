import mappy as mp
import logging
from pathlib import Path
from typing import List, Tuple, Optional
from opsinth.resources import *
from opsinth.utils import *
from opsinth.analysis_variants import genotype_known_variants, format_variant_classification
from opsinth.resources import *

# Create module logger
logger = logging.getLogger('opsinth.analysis_genes')

def run_find_genes(results, out_prefix):
    """
    Find and classify OPN1 genes in the polished sequence.
    
    Args:
        results: Dictionary containing analysis results including polished sequence
        out_prefix: Prefix for output files
        
    Returns:
        List of dictionaries containing genotyped alignments
    """
    logger.info("Starting OPN1 gene analysis")
    
    # Check if polished sequence exists
    if 'seq_polished' not in results:
        logger.error("No polished sequence found in results")
        return []
    
    # Write polished sequence to FASTA
    polished_fasta = Path(f"{out_prefix}.polished.fasta")
    with open(polished_fasta, 'w') as f:
        f.write(f">polished_sequence\n{results['seq_polished']}\n")
    logger.info(f"Wrote polished sequence to {polished_fasta}")
    
    # Get sequences
    seq = results.get('seq_polished')
    logger.info("Reading reference sequences")
    opn1mw_seq = read_fasta_seq(OPN1MW_REF)
    opn1lw_seq = read_fasta_seq(OPN1LW_REF)

    # Align sequences
    logger.info("Aligning OPN1MW reference")
    aln_opn1mw = align_sequences(str(polished_fasta), OPN1MW_REF, preset="splice:hq")
    logger.info("Aligning OPN1LW reference")
    aln_opn1lw = align_sequences(str(polished_fasta), OPN1LW_REF, preset="splice:hq")

    # Filter alignments
    aln_opn1mw_filtered = filter_alignments(aln_opn1mw)
    aln_opn1lw_filtered = filter_alignments(aln_opn1lw)

    logger.info(f"Found {len(aln_opn1mw_filtered)} MW and {len(aln_opn1lw_filtered)} LW alignments")

    # Analyze variants
    genotyped_alignments = []
    
    # Process MW alignments
    for aln in aln_opn1mw_filtered:
        logger.info(f"\nAnalyzing MW alignment {aln['id']} ({aln['ref_start']}-{aln['ref_end']})")
        
        # Analyze exon 5 color-determining variants
        logger.info("Analyzing exon 5 color-determining variants")
        exon5_vars = genotype_known_variants(aln, seq, OPSIN_EXON5_COLOR_VARIANTS)
        logger.info(format_variant_classification(exon5_vars))
        
        # Analyze exon 3 splicing variants
        logger.info("Analyzing exon 3 splicing variants")
        exon3_vars = genotype_known_variants(aln, seq, OPSIN_EXON3_SPLICING_VARIANTS)
        logger.info(format_variant_classification(exon3_vars))

        genotyped_alignments.append({
            'alignment': aln,
            'type': 'MW',
            'exon5_variants': exon5_vars,
            'exon3_variants': exon3_vars
        })
    
    # Process LW alignments
    for aln in aln_opn1lw_filtered:
        logger.info(f"\nAnalyzing LW alignment {aln['id']} ({aln['ref_start']}-{aln['ref_end']})")
        
        # Analyze exon 5 color-determining variants
        logger.info("Analyzing exon 5 color-determining variants")
        exon5_vars = genotype_known_variants(aln, seq, OPSIN_EXON5_COLOR_VARIANTS)
        logger.info(format_variant_classification(exon5_vars))
        
        # Analyze exon 3 splicing variants
        logger.info("Analyzing exon 3 splicing variants")
        exon3_vars = genotype_known_variants(aln, seq, OPSIN_EXON3_SPLICING_VARIANTS)
        logger.info(format_variant_classification(exon3_vars))

        genotyped_alignments.append({
            'alignment': aln,
            'type': 'LW',
            'exon5_variants': exon5_vars,
            'exon3_variants': exon3_vars
        })

    logger.info(f"Completed analysis of {len(genotyped_alignments)} alignments")
    return genotyped_alignments

def align_sequences(reference: str | Path, query: str | Path, preset: str = "map-ont") -> List[dict]:
    """
    Align sequences using minimap2. Can handle both file paths and sequences directly.
    
    Args:
        reference: Either path to reference FASTA file or nucleotide sequence
        query: Either path to query FASTA file or nucleotide sequence
        preset: Minimap2 preset (default: map-ont)
        
    Returns:
        List of dictionaries containing alignment information:
        {
            'reference_name': Name of reference sequence,
            'query_name': Name of query sequence,
            'ref_start': Start position on reference,
            'ref_end': End position on reference,
            'query_start': Start position on query,
            'query_end': End position on query,
            'mapping_quality': Mapping quality score (0-60),
            'cigar': CIGAR string describing alignment operations,
            'strand': Strand orientation (+ or -),
            'sequence_identity': Alignment identity score (0-1)
        }
    """
    try:
        # Convert paths to strings if they are Path objects
        reference = str(reference)
        query = str(query)
        
        # Create aligner object with eqx flag
        aligner = mp.Aligner(reference, preset=preset, extra_flags=0x4000000)
        if not aligner:
            raise Exception("Failed to load reference file")

        results = []
        
        # Read query sequences and perform alignment
        for query_name, seq, qual in mp.fastx_read(query):
            if not seq:
                continue
                
            # Perform alignment
            for hit in aligner.map(seq):
                results.append({
                    'reference_name': hit.ctg,
                    'query_name': query_name,
                    'ref_start': hit.r_st,
                    'ref_end': hit.r_en,
                    'query_start': hit.q_st,
                    'query_end': hit.q_en,
                    'mapping_quality': hit.mapq,
                    'cigar': hit.cigar_str,
                    'strand': '-' if hit.strand == -1 else '+',
                    'sequence_identity': hit.mlen / hit.blen if hit.blen > 0 else 0
                })
        
        return results
    
    except Exception as e:
        logging.error(f"Error during sequence alignment: {str(e)}")
        raise
        
def filter_alignments(alignments: List[dict]) -> List[dict]:
    """
    Filter and sort alignments, adding sequential IDs.
    
    Args:
        alignments: List of alignment dictionaries
        
    Returns:
        Filtered and sorted list of alignments with IDs
    """
    if not alignments:
        logging.warning("No alignments to filter")
        return []
    
    # Log the first alignment info for debugging
    first_aln = alignments[0]
    logging.info(f"Filtering alignments for {first_aln['query_name']} - {first_aln['reference_name']}")
    
    # Filter by sequence identity
    alignments_filtered = [a for a in alignments if a['sequence_identity'] > 0.9]
    logging.info(f"Kept {len(alignments_filtered)} alignments of {len(alignments)} after filtering")

    # Sort alignments by reference start position
    alignments_filtered.sort(key=lambda x: x['ref_start'])

    # Assign sequential IDs
    for i, alignment in enumerate(alignments_filtered):
        alignment['id'] = f"opn1_locus_{i + 1:02d}"

    # Log sorted alignments
    logging.debug("\nSorted alignments with IDs:")
    for aln in alignments_filtered:
        logging.debug(f"ID: {aln['id']}, Position: {aln['ref_start']}-{aln['ref_end']}")

    return alignments_filtered

def format_alignment_results(results: List[dict]) -> str:
    """
    Format alignment results into a readable string.
    
    Args:
        results: List of alignment result dictionaries
        
    Returns:
        Formatted string containing alignment information
    """
    output = []
    output.append("Alignment Results:")
    output.append("-----------------")
    
    for hit in results:
        output.append(f"Reference: {hit['reference_name']}")
        output.append(f"Query: {hit['query_name']}")
        output.append(f"Reference Position: {hit['ref_start']}-{hit['ref_end']}")
        output.append(f"Query Position: {hit['query_start']}-{hit['query_end']}")
        output.append(f"Mapping Quality: {hit['mapping_quality']}")
        output.append(f"Sequence Identity: {hit['sequence_identity']:.2%}")
        output.append(f"CIGAR: {hit['cigar']}")
        output.append(f"Strand: {hit['strand']}")
        output.append("-----------------")
    
    return "\n".join(output)

def write_bed_file(results: List[dict], output_prefix: str):
    """
    Write alignment results to BED format file.
    
    Args:
        results: List of alignment result dictionaries
        output_prefix: Prefix for output files
    """
    bed_path = output_prefix + '.bed'
    
    with open(bed_path, 'w') as f:
        for hit in results:
            # BED format: chrom start end name score strand
            bed_line = (f"{hit['reference_name']}\t{hit['ref_start']}\t{hit['ref_end']}\t"
                       f"{hit['query_name']}_{hit['query_start']}-{hit['query_end']}\t"
                       f"{hit['mapping_quality']}\t{hit['strand']}\n")
            f.write(bed_line)
    
    logging.info(f"BED file written to {bed_path}")



def format_pairwise_alignment(query_seq: str, ref_seq: str, cigar: str, 
                            query_start: int = 0, ref_start: int = 0,
                            line_width: int = 60) -> str:
    """
    Format a pairwise alignment in BLAST-like format.
    
    Args:
        query_seq: Query sequence
        ref_seq: Reference sequence
        cigar: CIGAR string from alignment
        query_start: Start position in query (default: 0)
        ref_start: Start position in reference (default: 0)
        line_width: Number of characters per line (default: 60)
        
    Returns:
        Formatted string showing the alignment
        
    Example:
        ref   : ACTG-TGCATGC
        match : |||| |||||||
        query : ACTGNTGCATGC
    """
    # Parse CIGAR string
    import re
    cigar_ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    
    # Initialize alignment strings
    ref_aln = []
    query_aln = []
    match_aln = []
    
    curr_query = query_start
    curr_ref = ref_start
    
    # Build alignment strings based on CIGAR operations
    for length, op in cigar_ops:
        length = int(length)
        
        if op in '=X':  # Match or mismatch
            q_seq = query_seq[curr_query:curr_query + length]
            r_seq = ref_seq[curr_ref:curr_ref + length]
            ref_aln.append(r_seq)
            query_aln.append(q_seq)
            # Add match symbols
            match_symbols = ''.join('|' if q == r else ' ' 
                                  for q, r in zip(q_seq, r_seq))
            match_aln.append(match_symbols)
            curr_query += length
            curr_ref += length
            
        elif op == 'M':  # Match or mismatch (old style)
            q_seq = query_seq[curr_query:curr_query + length]
            r_seq = ref_seq[curr_ref:curr_ref + length]
            ref_aln.append(r_seq)
            query_aln.append(q_seq)
            # Add match symbols
            match_symbols = ''.join('|' if q == r else ' ' 
                                  for q, r in zip(q_seq, r_seq))
            match_aln.append(match_symbols)
            curr_query += length
            curr_ref += length
            
        elif op == 'I':  # Insertion to reference
            ref_aln.append('-' * length)
            query_aln.append(query_seq[curr_query:curr_query + length])
            match_aln.append(' ' * length)
            curr_query += length
            
        elif op == 'D':  # Deletion from reference
            ref_aln.append(ref_seq[curr_ref:curr_ref + length])
            query_aln.append('-' * length)
            match_aln.append(' ' * length)
            curr_ref += length
            
        elif op == 'S':  # Soft clipping
            # Skip soft clipped bases
            curr_query += length
    
    # Join alignment strings
    ref_str = ''.join(ref_aln)
    match_str = ''.join(match_aln)
    query_str = ''.join(query_aln)
    
    # Format output with line wrapping
    output = []
    for i in range(0, len(ref_str), line_width):
        chunk_end = i + line_width
        ref_chunk = ref_str[i:chunk_end]
        match_chunk = match_str[i:chunk_end]
        query_chunk = query_str[i:chunk_end]
        
        # Add position numbers
        ref_pos = ref_start + i
        query_pos = query_start + i
        
        output.extend([
            f"ref  {ref_pos:>8}: {ref_chunk} :{ref_pos + len(ref_chunk)}",
            f"               {match_chunk}",
            f"query{query_pos:>8}: {query_chunk} :{query_pos + len(query_chunk)}",
            ""
        ])
    
    return "\n".join(output)


def get_alignment_sequence(alignment: dict, ref_seq: str, query_seq: str) -> str:
    """
    Get formatted pairwise alignment for an alignment result.
    
    Args:
        alignment: Dictionary containing alignment information
        ref_seq: Full reference sequence
        query_seq: Full query sequence
        
    Returns:
        Formatted alignment string
    """
    return format_pairwise_alignment(
        query_seq=query_seq,
        ref_seq=ref_seq,
        cigar=alignment['cigar'],
        query_start=alignment['query_start'],
        ref_start=alignment['ref_start']
    )

def main(reference_path: str, query_path: str, output_prefix: Optional[str] = None):
    """
    Main function to run sequence alignment and output results.
    
    Args:
        reference_path: Path to reference FASTA file
        query_path: Path to query FASTA file
        output_prefix: Optional prefix for output files
    """
    try:
        # Verify input files exist
        if not Path(reference_path).exists():
            raise FileNotFoundError(f"Reference file not found: {reference_path}")
        if not Path(query_path).exists():
            raise FileNotFoundError(f"Query file not found: {query_path}")
            
        # Perform alignment
        results = align_sequences(reference_path, query_path)
        
        # Read sequences
        ref_seq = read_fasta_seq(reference_path)
        query_seq = read_fasta_seq(query_path)
        
        # Format results
        formatted_output = format_alignment_results(results)
        
        # Output results
        if output_prefix:
            with open(output_prefix + '.alignments.txt', 'w') as f:
                f.write(formatted_output)
                f.write("\n\nDetailed Alignments:\n")
                f.write("===================\n\n")
                for aln in results:
                    f.write(f"Alignment for {aln['query_name']} to {aln['reference_name']}:\n")
                    f.write(get_alignment_sequence(aln, ref_seq, query_seq))
                    f.write("\n" + "=" * 80 + "\n\n")
                    
            logging.info(f"Results written to {output_prefix}.alignments.txt")
            write_bed_file(results, output_prefix)
        else:
            print(formatted_output)
            
    except Exception as e:
        logging.error(f"Error in main function: {str(e)}")
        raise

def classify_opsin_type(alignment: dict, ref_seq: str) -> dict:
    """
    Classify the opsin type based on color-determining variants.
    
    Args:
        alignment: Dictionary containing alignment information
        ref_seq: Reference sequence
        
    Returns:
        Dictionary containing:
        - id: Alignment ID
        - region: Tuple of (start, end) positions
        - variants: List of variant dictionaries
        - ref_count: Number of reference (MW) matches
        - alt_count: Number of alternative (LW) matches
        - opsin_type: Classification result ("OPN1MW (Green)", "OPN1LW (Red)", or "Undetermined")
    """
    variants = genotype_known_variants(alignment, ref_seq, OPSIN_EXON5_COLOR_VARIANTS)
    
    ref_count = sum(1 for var in variants if var['is_ref'])
    alt_count = sum(1 for var in variants if var['is_alt'])
    
    # Determine type based on majority of matches
    if ref_count >= 2:
        opsin_type = "OPN1MW (Green)"
    elif alt_count >= 2:
        opsin_type = "OPN1LW (Red)"
    else:
        opsin_type = "Undetermined"
    
    return {
        'id': alignment['id'],
        'region': (alignment['ref_start'], alignment['ref_end']),
        'variants': variants,
        'ref_count': ref_count,
        'alt_count': alt_count,
        'opsin_type': opsin_type
    }



if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Align sequences using minimap2")
    parser.add_argument("--reference", required=True, help="Path to reference FASTA file")
    parser.add_argument("--query", required=True, help="Path to query FASTA file")
    parser.add_argument("--output", help="Optional path to output file")
    
    args = parser.parse_args()
    
    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    main(args.reference, args.query, args.output)
