import sys
import os 
import logging
import json
import argparse

# Use absolute imports instead of relative imports (remove the dots)
from opsinth.analysis_genes import run_find_genes
from opsinth.igv_web import create_igv_html, open_igv_viewer
from opsinth.config import DEFAULTS, VARIANTS, VERSION
from opsinth.analysis_sequence import run_ref_analysis

import logging
from IPython.display import display, HTML

# Configure logging for notebook display
class NotebookHandler(logging.Handler):
    def emit(self, record):
        # Convert log record to HTML with different colors based on level
        level_colors = {
            'DEBUG': 'grey',
            'INFO': 'black',
            'WARNING': 'orange',
            'ERROR': 'red',
            'CRITICAL': 'darkred'
        }
        color = level_colors.get(record.levelname, 'black')
        msg = f'<pre style="margin:0; color:{color};">{self.format(record)}</pre>'
        display(HTML(msg))

# Set up the logger
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)  # Set to DEBUG if you want to see debug messages

# Remove any existing handlers
for handler in logger.handlers[:]:
    logger.removeHandler(handler)

# Add the notebook handler
notebook_handler = NotebookHandler()
formatter = logging.Formatter('%(levelname)s: %(message)s')
notebook_handler.setFormatter(formatter)
logger.addHandler(notebook_handler)

# Now your logging output will appear in the notebook
# Example usage:
logging.info("Starting analysis...")

out_prefix = "nb_test_output"

bam = "input/wgs_short.locus.bam"
bed = "input/opsin_region.bed"
ref = "../data/grch38.chrX.roi.fa.gz"
anchors = "input/anchors.fasta"
json_file = "input/wgs_short.paraphase.json"
fasta_seq_file = "input/opsin_nodes_seq.fa"


from opsinth.analysis_sequence import run_ref_analysis

# Step1: Refbased analysis
results_ref = run_ref_analysis(bam, bed, ref, anchors)

from opsinth.analysis_sequence import create_draft_from_overlaps, create_draft_from_graph
seq = create_draft_from_graph(results_ref)

