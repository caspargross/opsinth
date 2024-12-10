# opsin_analysis/main.py
VERSION = "0.1" 

import argparse
import os
import logging
from opsin_analysis.analysis import *
from opsin_analysis.plotting import *
from opsin_analysis.utils import configure_logging
from opsin_analysis.igv_viewer import create_igv_html, open_igv_viewer

def main():
    parser = argparse.ArgumentParser(description="Opsin Analysis CLI")
    parser.add_argument('--bam', required=True, help='Path to BAM file')
    parser.add_argument('--bed', required=True, help='Path to BED file')
    parser.add_argument('--ref', required=True, help='Path to reference genome (FASTA)')
    parser.add_argument('--anchors', required=True, help='Path to anchors FASTA file')
    parser.add_argument('--out', required=True, help='Output prefix (can be folder)')
    parser.add_argument('-v', '--verbose', action='count', default=0,
                      help='Increase verbosity level (e.g., -v, -vv, -vvv)')
    parser.add_argument('--no-igv', action='store_true', help='Do not start IGV.js viewer')  # New option
    
    args = parser.parse_args()
    
    # Configure logging based on verbosity level
    configure_logging(args.verbose)
    
    logging.info("Starting Opsin Analysis")

    # Read files
    dataset = read_files(args.bam, args.bed, args.ref, args.anchors)

    # Run analysis
    results = run_analysis(**dataset)
        
    # Create output directory from path component of output prefix
    output_dir = os.path.dirname(args.out)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    plot_coverage(results, output_dir)
    plot_alignment_quality(results, output_dir)
    write_bam_file(results, dataset.get('reads'), args.out, args.bam, VERSION)

    if not args.no_igv:
        # Convert ROI list to string format for IGV.js
        roi_list = dataset.get('roi', [['chrX', 154143438, 154295680]])  # Default values if 'roi' is not set
        target_region = f"{roi_list[0][0]}:{roi_list[0][1]}-{roi_list[0][2]}"
        logging.debug(f"Formatted Target Region: {target_region}")
        
        # Create IGV HTML file
        create_igv_html(output_dir, args.out, dataset.get('ref', 'GRCh38'), target_region)

        # Open IGV viewer
        open_igv_viewer(output_dir)

    logging.info("Completed successfully")

if __name__ == "__main__":
    main()