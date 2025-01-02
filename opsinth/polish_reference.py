import os
import subprocess
import logging
from opsinth.utils import *
from opsinth.analysis import *

def run_polish_denovo(results_denovo, out_prefix, n_polish_rounds=1, delete_intermediate_files=True):
    
    logging.info(f"Start {n_polish_rounds} rounds of polishing")

    roi = results_denovo.get('roi')
    reads_aligned = results_denovo.get('reads_aligned')
    reads = results_denovo.get('reads')
    seq_unpolished = results_denovo.get('seq_denovo')
      
    polish_stats = {}

    for i in range(n_polish_rounds):

        # Prepare input files for racon
        polish_alignment_file = f"{out_prefix}.polish_round_{i}.aln.sam"
        polish_ref_file = f"{out_prefix}.polish_round_{i}.seq.fasta"
        polish_query_file = f"{out_prefix}.polish_round_{i}.reads.fastq"

        write_bam(reads_aligned, reads, roi, polish_alignment_file, output_format_bam=False)
        write_fasta(seq_unpolished, roi, polish_ref_file)
        write_fastq(reads_aligned, polish_query_file)
        
        logging.debug(f"Running racon with {polish_query_file}, {polish_alignment_file}, {polish_ref_file}")

        # Run racon
        seq_polished = run_racon(polish_query_file, polish_alignment_file, polish_ref_file)
        write_fasta(seq_polished, roi, (out_prefix + ".polished"))

        polish_stats[i] = evaluate_racon_improvement(seq_unpolished, seq_polished)
        
        logging.info(f"Polish round {i} / {n_polish_rounds} completed")

        # Update input files for next round
        if i < n_polish_rounds - 1:
            seq_unpolished = seq_polished
            reads_aligned = align_reads_to_ref(results_denovo.get('reads_aligned'), results_denovo.get('no_anchor_reads'), results_denovo.get('unique_anchor_reads'), results_denovo.get('double_anchor_reads'), results_denovo.get('anchors_on_reads'), seq_polished, results_denovo.get('anchors_on_ref'))

        if delete_intermediate_files and i < n_polish_rounds - 1:
            os.remove(polish_alignment_file)
            os.remove(polish_ref_file)
            os.remove(polish_query_file)
    
    # Update results_denovo with polished sequence  
    results_denovo['seq_polished'] = seq_polished
    results_denovo['polish_stats'] = polish_stats
    results_denovo['reads_aligned'] = reads_aligned

    return results_denovo

def run_racon(reads_file, sam_file, ref_file, threads=1, racon_path="lib/racon"):
    """Polish the alignment using racon."""
    try:
        cmd = [
            racon_path,
            "-t", str(threads),
            reads_file,
            sam_file,
            ref_file
        ]
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        logging.debug(f"racon log: {result.stderr}")
        logging.info("racon completed successfully.")
        seq = result.stdout.decode().splitlines()[1]
        return seq

    except subprocess.CalledProcessError as e:
        logging.error(f"racon failed: {e.stderr.decode()}")
        raise

def evaluate_racon_improvement(seq_unpolished, seq_polished):
    pass