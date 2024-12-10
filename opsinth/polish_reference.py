import mappy as mp
import subprocess
import logging
import os

def run_minimap2(reads_file, ref_file, output_file, threads=4):
    """Map reads to a reference sequence using mappy (minimap2)."""
    try:
        # Open the reference file with mappy
        aligner = mp.Aligner(ref_file, preset='map-ont', n_threads=threads)
        if not aligner:
            raise Exception("Failed to load/build index for reference file.")

        # Open the reads file and output SAM file
        with open(reads_file, 'r') as reads, open(output_file, 'w') as out:
            for name, seq, qual in mp.fastx_read(reads):
                for hit in aligner.map(seq):
                    # Write each alignment to the SAM file
                    out.write(f"{name}\t{hit.flag}\t{hit.ctg}\t{hit.r_st + 1}\t{hit.mapq}\t{hit.cigar_str}\t*\t0\t0\t{seq}\t{qual}\n")
        
        logging.info("minimap2 (mappy) completed successfully.")
    except Exception as e:
        logging.error(f"minimap2 (mappy) failed: {str(e)}")
        raise

def run_racon(reads_file, sam_file, ref_file, output_file, threads=4):
    """Polish the alignment using racon."""
    try:
        cmd = [
            "racon",
            "-t", str(threads),
            reads_file,
            sam_file,
            ref_file,
            output_file
        ]
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        logging.info("racon completed successfully.")
    except subprocess.CalledProcessError as e:
        logging.error(f"racon failed: {e.stderr.decode()}")
        raise

def polish_reference(reads_file, ref_file, output_prefix, threads=4):
    """Interface to map reads and polish the reference sequence."""
    sam_file = f"{output_prefix}.sam"
    polished_ref_file = f"{output_prefix}_polished.fasta"

    logging.info("Starting minimap2 to map reads to reference.")
    run_minimap2(reads_file, ref_file, sam_file, threads)

    logging.info("Starting racon to polish the reference.")
    run_racon(reads_file, sam_file, ref_file, polished_ref_file, threads)

    logging.info(f"Polished reference sequence saved to {polished_ref_file}")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Polish a draft reference sequence using minimap2 and racon.")
    parser.add_argument('--reads', required=True, help='Path to reads file (FASTQ)')
    parser.add_argument('--ref', required=True, help='Path to draft reference sequence (FASTA)')
    parser.add_argument('--out', required=True, help='Output prefix for generated files')
    parser.add_argument('--threads', type=int, default=4, help='Number of threads to use')
    parser.add_argument('-v', '--verbose', action='count', default=0, help='Increase verbosity level (e.g., -v, -vv, -vvv)')

    args = parser.parse_args()

    # Configure logging
    log_level = logging.WARNING if args.verbose == 0 else logging.INFO if args.verbose == 1 else logging.DEBUG
    logging.basicConfig(level=log_level, format='%(asctime)s - %(levelname)s - %(message)s')

    try:
        polish_reference(args.reads, args.ref, args.out, args.threads)
    except Exception as e:
        logging.error(f"An error occurred: {e}") 