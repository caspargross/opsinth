# opsin_analysis/main.py
import argparse
from opsin_analysis.utils import *
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
    results = run_analysis(dataset)
    
    # Generate plots
    plot_alignment_quality(results, args.out)
    plot_coverage(results, args.out)

    # Write bamfile
    write_bam_file(results, args.out)

if __name__ == "__main__":
    main()