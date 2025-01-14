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
from opsinth.config import DEFAULTS, VARIANTS

class Opsinth:
    def __init__(self):
        self.parser = argparse.ArgumentParser(description="Opsin Analysis CLI")
        self.parser.add_argument('--bam', required=True, 
                               help='Path to BAM file')
        self.parser.add_argument('--out', required=False, 
                               default=DEFAULTS['out'],
                               help='Output prefix (can be folder)')
        self.parser.add_argument('--bed', required=False, 
                               default=DEFAULTS['bed'],
                               help='Path to BED file')
        self.parser.add_argument('--ref', required=False, 
                               default=DEFAULTS['ref'],
                               help='Path to reference genome (FASTA)')
        self.parser.add_argument('--anchors', required=False, 
                               default=DEFAULTS['anchors'],
                               help='Path to anchors FASTA file')
        self.parser.add_argument('--racon', required=False, 
                               default=DEFAULTS['racon'],
                               help='Path to racon executable')
        self.parser.add_argument('--n_polish_rounds', required=False, 
                               default=DEFAULTS['n_polish_rounds'],
                               help='Number of polishing rounds')
        self.parser.add_argument('--debug_plots', 
                                default=DEFAULTS['debug_plots'],
                                help='Create debug plots')
        self.parser.add_argument('--export_unpolished',
                                default=DEFAULTS['export_unpolished'],
                                help="Export unpolished sequence and alignments to file")
        self.parser.add_argument('--export_reference_based',
                                default=DEFAULTS['export_refbased'],
                                help="Export reference based alignments (without denovo)")
        self.parser.add_argument('-v', '--verbose', action='count', 
                               default=1,
                               help='Set verbosity level (e.g., -v, -vv, -vvv) Default: -vv')
        self.parser.add_argument('--open-igv', action='store_true', 
                               help='Open IGV.js viewer')


    def run(self):
        args = self.parser.parse_args()
        
        # Configure logging based on verbosity level
        configure_logging(args.verbose)
        logger = logging.getLogger('opsinth')
        logger.info(f"OPSINTH v{VERSION}")

        # Determine output directory and prefix
        output_dir = os.path.dirname(args.out)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            out_prefix = args.out
        else:
            out_prefix = args.out

        # Step1: Refbased analysis
        results_ref = run_ref_analysis(args.bam, args.bed, args.ref, args.anchors)

        if args.export_reference_based:
            write_bam(results_ref.get('reads_aligned'), results_ref.get('reads'), results_ref.get('roi'), (out_prefix + ".refbased"), VERSION, args.bam)
            if args.debug_plots:
                plot_coverage(results_ref, (out_prefix + ".refbased"))
                plot_alignment_quality(results_ref, (out_prefix + ".refbased"))

        # Step2: Denovo analysis
        results_denovo = run_denovo_analysis(results_ref.copy())

        if args.export_unpolished:
            write_fasta(results_denovo.get('seq_denovo'), results_denovo.get('roi'), (out_prefix + ".unpolished"))
            write_bam(results_denovo.get('reads_aligned'), results_denovo.get('reads'), results_denovo.get('roi'), (out_prefix + ".unpolished"), VERSION)
            if args.debug_plots:
                plot_coverage(results_denovo, (out_prefix + ".unpolished"))
                plot_alignment_quality(results_denovo, (out_prefix + ".unpolished"))

        # Step3: Polishing
        results_polished = run_polish_denovo(results_denovo.copy(), out_prefix, n_polish_rounds=args.n_polish_rounds, racon_path=args.racon)

        write_fasta(results_polished.get('seq_polished'), results_polished.get('roi'), (out_prefix))
        write_bam(results_polished.get('reads_aligned'), results_polished.get('reads'), results_polished.get('roi'), (out_prefix), VERSION)
        if args.debug_plots:
            plot_coverage(results_polished, out_prefix)
            plot_alignment_quality(results_polished, out_prefix)
            plot_polish_stats(results_polished.get('polish_stats'), out_prefix)

        # Step4: Find genes and haplotypes
        results_genes = run_find_genes(results_polished.copy(), out_prefix)

        with open(out_prefix + ".genes.json", "w") as f:
            json.dump(results_genes, f, indent=4)

        create_igv_session(f"{out_prefix}.fasta", f"{out_prefix}.bam", results_polished.get('roi'), out_prefix, DEFAULTS['igv_session_template'])

        # #if not args.no_igv:
        #     # Convert ROI list to string format for IGV.js
        #     roi_list = dataset.get('roi', [['chrX', 154143438, 154295680]])  # Default values if 'roi' is not set
        #     target_region = f"{roi_list[0][0]}:{roi_list[0][1]}-{roi_list[0][2]}"
        #     logging.debug(f"Formatted Target Region: {target_region}")
            
        #     # Create IGV HTML file
        #     create_igv_html(output_dir, os.path.basename(out_prefix), dataset.get('ref', 'hg38'), target_region)

        #     # Open IGV viewer
        #     open_igv_viewer(output_dir)
    
        logger.info("Opsinth completed successfully")
