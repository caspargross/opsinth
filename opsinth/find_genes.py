import mappy as mp
import logging
from pathlib import Path
from typing import List, Tuple, Optional

def align_sequences(reference_path: str, query_path: str, preset: str = "map-ont") -> List[Tuple[str, str, int, int, int, int, int, str, str]]:
    """
    Align sequences using minimap2 and return alignment results.
    
    Args:
        reference_path: Path to reference FASTA file
        query_path: Path to query FASTA file
        
    Returns:
        List of tuples containing (reference_name, query_name, ref_start, ref_end, 
                                 query_start, query_end, mapping_quality, cigar, strand)
    """
    try:
        # Create aligner object
        aligner = mp.Aligner(reference_path, preset=preset)
        if not aligner:
            raise Exception("Failed to load reference file")

        results = []
        
        # Read query sequences and perform alignment
        for query_name, seq, qual in mp.fastx_read(query_path): # qual will be None for FASTA
            if not seq:
                continue
                
            # Perform alignment
            for hit in aligner.map(seq):
                results.append((
                    hit.ctg,                # Reference sequence name
                    query_name,             # Query sequence name
                    hit.r_st,              # Start position on reference
                    hit.r_en,              # End position on reference
                    hit.q_st,              # Start position on query
                    hit.q_en,              # End position on query
                    hit.mapq,              # Mapping quality score
                    hit.cigar_str,         # CIGAR string
                    '-' if hit.strand == -1 else '+',  # Strand
                ))
                
        return results
    
    except Exception as e:
        logging.error(f"Error during sequence alignment: {str(e)}")
        raise

def format_alignment_results(results: List[Tuple[str, str, int, int, int, int, int, str, str]]) -> str:
    """
    Format alignment results into a readable string.
    
    Args:
        results: List of alignment result tuples
        
    Returns:
        Formatted string containing alignment information
    """
    output = []
    output.append("Alignment Results:")
    output.append("-----------------")
    
    for ref_name, query_name, ref_start, ref_end, q_start, q_end, mapq, cigar, strand in results:
        output.append(f"Reference: {ref_name}")
        output.append(f"Query: {query_name}")
        output.append(f"Reference Position: {ref_start}-{ref_end}")
        output.append(f"Query Position: {q_start}-{q_end}")
        output.append(f"Mapping Quality: {mapq}")
        output.append(f"CIGAR: {cigar}")
        output.append(f"Strand: {strand}")
        output.append("-----------------")
    
    return "\n".join(output)

def write_bed_file(results: List[Tuple[str, str, int, int, int, int, int, str, str]], output_prefix: str):
    """
    Write alignment results to BED format file.
    
    Args:
        results: List of alignment result tuples
        output_prefix: Prefix for output files
    """
    # Remove .txt extension if present and add .bed
    bed_path = output_prefix + '.bed'
    
    with open(bed_path, 'w') as f:
        for ref_name, query_name, ref_start, ref_end, q_start, q_end, mapq, cigar, strand in results:
            # BED format: chrom start end name score strand
            bed_line = (f"{ref_name}\t{ref_start}\t{ref_end}\t"
                       f"{query_name}_{q_start}-{q_end}\t"
                       f"{mapq}\t{strand}\n")
            f.write(bed_line)
    
    logging.info(f"BED file written to {bed_path}")

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
        
        # Format results
        formatted_output = format_alignment_results(results)
        
        # Output results
        if output_prefix:
            with open(output_prefix + '.alignments.txt', 'w') as f:
                f.write(formatted_output)
            logging.info(f"Results written to {output_prefix}.alignments.txt")
            write_bed_file(results, output_prefix)
        else:
            print(formatted_output)
            
    except Exception as e:
        logging.error(f"Error in main function: {str(e)}")
        raise

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
