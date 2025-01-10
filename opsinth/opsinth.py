# opsinth/opsinth.py

import argparse
import os
import logging
import json
from opsinth.analysis_sequence import run_ref_analysis, run_denovo_analysis
from opsinth.polish import run_polish_denovo
from opsinth.analysis_genes import run_find_genes
from opsinth.plots import *
from opsinth.utils import configure_logging
from opsinth.igv_web import create_igv_html, open_igv_viewer
from opsinth.resources import *

class Opsinth:
    def __init__(self):
        self.parser = argparse.ArgumentParser(description="Opsin Analysis CLI")
        self.parser.add_argument('--bam', required=True, help='Path to BAM file')
        self.parser.add_argument('--bed', required=False, default=DEFAULT_BED,
                                help='Path to BED file (default: opsin_region.bed)')
        self.parser.add_argument('--ref', required=False, default=DEFAULT_REF,
                                help='Path to reference genome (FASTA) (default: grch38 reduced)')
        self.parser.add_argument('--anchors', required=False, default=DEFAULT_ANCHORS,
                                help='Path to anchors FASTA file (default: wgs_anchors)')
        self.parser.add_argument('--out', required=False, help='Output prefix (can be folder)', default = DEFAULT_OUT)
        self.parser.add_argument('-v', '--verbose', action='count', default=0,
                                 help='Increase verbosity level (e.g., -v, -vv, -vvv)')
        self.parser.add_argument('--no-igv', action='store_true', help='Do not start IGV.js viewer')


    def run(self):
        args = self.parser.parse_args()
        
        # Configure logging based on verbosity level
        configure_logging(args.verbose)
        
        logging.info(f"OPSINTH v{VERSION}")

        # Determine output directory and prefix        # Determine output directory and prefix
        output_dir = os.path.dirname(args.out)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            out_prefix = os.path.join(output_dir, "opsinth")
        else:
            out_prefix = args.out

        # Run analysis
        results_ref = run_ref_analysis(args.bam, args.bed, args.ref, args.anchors)
        results_denovo = run_denovo_analysis(results_ref.copy())
        results_polished = run_polish_denovo(results_denovo.copy(), out_prefix)
        results_genes = run_find_genes(results_polished.copy(), out_prefix)

        # Write results_genes to file
        with open(out_prefix + ".genes.json", "w") as f:
            json.dump(results_genes, f, indent=4)

        # Plot and write reference results
        plot_coverage(results_ref, (out_prefix + ".ref"))
        plot_alignment_quality(results_ref, (out_prefix + ".ref"))
        write_bam(results_ref.get('reads_aligned'), results_ref.get('reads'), results_ref.get('roi'), (out_prefix + ".ref"), VERSION, args.bam)

        # Plot and write denovo results
        write_fasta(results_denovo.get('seq_denovo'), results_denovo.get('roi'), (out_prefix + ".denovo"))
        plot_coverage(results_denovo, (out_prefix + ".denovo"))
        plot_alignment_quality(results_denovo, (out_prefix + ".denovo"))
        write_bam(results_denovo.get('reads_aligned'), results_denovo.get('reads'), results_denovo.get('roi'), (out_prefix + ".denovo"), VERSION)

        # Plot and write polished results
        write_fasta(results_polished.get('seq_polished'), results_polished.get('roi'), (out_prefix + ".denovo.polished"))
        plot_coverage(results_polished, (out_prefix + ".denovo.polished"))
        plot_alignment_quality(results_polished, (out_prefix + ".denovo.polished"))
        write_bam(results_polished.get('reads_aligned'), results_polished.get('reads'), results_polished.get('roi'), (out_prefix + ".denovo.polished"), VERSION)
    
        # #if not args.no_igv:
        #     # Convert ROI list to string format for IGV.js
        #     roi_list = dataset.get('roi', [['chrX', 154143438, 154295680]])  # Default values if 'roi' is not set
        #     target_region = f"{roi_list[0][0]}:{roi_list[0][1]}-{roi_list[0][2]}"
        #     logging.debug(f"Formatted Target Region: {target_region}")
            
        #     # Create IGV HTML file
        #     create_igv_html(output_dir, os.path.basename(out_prefix), dataset.get('ref', 'hg38'), target_region)

        #     # Open IGV viewer
        #     open_igv_viewer(output_dir)

        logging.info("Completed successfully")