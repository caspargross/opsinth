# opsin_analysis/main.py
VERSION = "0.1" 

import argparse
import os
import logging
from opsin_analysis.analysis import *
from opsin_analysis.plotting import *
from opsin_analysis.utils import configure_logging

def main():
    parser = argparse.ArgumentParser(description="Opsin Analysis CLI")
    parser.add_argument('--bam', required=True, help='Path to BAM file')
    parser.add_argument('--bed', required=True, help='Path to BED file')
    parser.add_argument('--ref', required=True, help='Path to reference genome (FASTA)')
    parser.add_argument('--anchors', required=True, help='Path to anchors FASTA file')
    parser.add_argument('--out', required=True, help='Output prefix (can be folder)')
    parser.add_argument('-v', '--verbose', action='count', default=0,
                      help='Increase verbosity level (e.g., -v, -vv, -vvv)') 
    
    args = parser.parse_args()
    
    # Configure logging based on verbosity level
    configure_logging(args.verbose)
    
    logging.info("Starting Opsin Analysis")

    # Read files
    dataset = read_files(args.bam, args.bed, args.ref, args.anchors)

    # Run analysis
    results = run_analysis(**dataset)
    
    # Generate outputs
    
    # Create output directory from path component of output prefix
    output_dir = os.path.dirname(args.out)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    plot_coverage(results, output_dir)
    plot_alignment_quality(results, output_dir)
    write_bam_file(results, dataset.get('reads'), args.out, args.bam, VERSION)

    logging.info("Completed successfully")

if __name__ == "__main__":
    main()