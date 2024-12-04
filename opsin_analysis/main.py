# opsin_analysis/main.py
import argparse
import os
from opsin_analysis.analysis import *
from opsin_analysis.plotting import *

def main():
    parser = argparse.ArgumentParser(description="Opsin Analysis CLI")
    parser.add_argument('--bam', required=True, help='Path to BAM file')
    parser.add_argument('--bed', required=True, help='Path to BED file')
    parser.add_argument('--ref', required=True, help='Path to reference genome (FASTA)')
    parser.add_argument('--anchors', required=True, help='Path to anchors FASTA file')
    parser.add_argument('--out', required=True, help='Output prefix (can be folder)')
    
    args = parser.parse_args()
    
    # Read files
    dataset = read_files(args.bam, args.bed, args.ref, args.anchors)

    # Run analysis
    results = run_analysis(**dataset)
    
    # Generate outputs
    
    # Create output directory from path component of output prefix
    output_dir = os.path.dirname(args.out)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    plot_coverage(results, args.out)
    plot_alignment_quality(results, args.out)
    write_bam_file(results, dataset.get('reads'), args.out, args.bam)

    print("Finished application")

if __name__ == "__main__":
    main()