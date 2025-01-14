import os
import subprocess
import logging
from opsinth.utils import *
from opsinth.analysis_sequence import *
from opsinth.plots import *
def run_polish_denovo(results, out_prefix, racon_path, n_polish_rounds, delete_intermediate_files=True):
    
    logging.info(f"Start {n_polish_rounds} rounds of polishing")
    logging.debug(f"Racon executable: {racon_path}")

    roi = results.get('roi')
    reads_aligned = results.get('reads_aligned') 
    reads = results.get('reads')
    anchors = results.get('anchors')
    seq_unpolished = results.get('seq_denovo')
      
    polish_stats = {}

    for i in range(n_polish_rounds):

        # Prepare input files for racon
        polish_alignment_file = f"{out_prefix}.polish_round_{i}.aln.sam"
        polish_ref_file = f"{out_prefix}.polish_round_{i}.seq.fasta"
        polish_query_file = f"{out_prefix}.polish_round_{i}.reads.fastq"

        write_bam(reads_aligned, reads, roi, polish_alignment_file, output_format_bam=False)
        #write_bam(reads_aligned, reads, roi, polish_alignment_file)
        write_fasta(seq_unpolished, roi, polish_ref_file)
        write_fastq(reads_aligned, polish_query_file)
        
        logging.debug(f"Running racon with {polish_query_file}, {polish_alignment_file}, {polish_ref_file}")

        # Run racon
        seq_polished = run_racon(polish_query_file, polish_alignment_file, ref_file = polish_ref_file, racon_path=racon_path, threads=1,  )

        # Evaluate sequence improvements in last round
        polish_stats[i] = evaluate_racon_improvement(seq_unpolished, seq_polished)
        
        logging.info(f"Polish round {i+1} / {n_polish_rounds} completed")

        # Update input files for next round
        if i < n_polish_rounds:
            seq_unpolished = seq_polished
                    
            # Align anchors to draft reference sequence
            anchors_on_ref = align_anchors_to_ref(anchors, seq_unpolished)

            min_pos = min([anchor['start'] for anchor in anchors_on_ref['anchor_positions'].values()])
            max_pos = max([anchor['end'] for anchor in anchors_on_ref['anchor_positions'].values()])
            roi = [[f"ref_denovo_polished", min_pos, max_pos]]

            # Find anchors on reads
            anchors_on_reads = find_anchors_on_reads(reads, anchors_on_ref) #TODO: Is this redundant? Alreadty calculated in ref analysis already ??

            # Characterize reads into no anchor, unique anchor and double anchor reads
            (no_anchor_reads, unique_anchor_reads, double_anchor_reads) = characterize_read_anchors(reads, anchors_on_ref, anchors_on_reads) # Redundant, was calculated in ref analysis already 
            
            reads_aligned = align_reads_to_ref(reads, no_anchor_reads, unique_anchor_reads, double_anchor_reads, anchors_on_reads, seq_unpolished, anchors_on_ref)
            reads_aligned = filter_reads_by_edit_distance(reads_aligned)

        if delete_intermediate_files:
            logging.debug(f"Deleting intermediate files for polish round {i}")
            os.remove(polish_alignment_file)
            os.remove(polish_ref_file)
            os.remove(polish_query_file)

        # Update results_denovo with polished sequence  
        results['seq_polished'] = seq_polished
        results['polish_stats'] = polish_stats
        results['roi'] = roi
        results['reads_aligned'] = reads_aligned

    return results

def run_racon(reads_file, sam_file, ref_file, racon_path, threads):
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
    # Perform global alignment
    alignment = edlib.align(seq_unpolished, seq_polished, mode="NW", task="path")
    
    # Parse CIGAR string to count operations
    cigar = alignment['cigar']
    
    # Initialize counters
    matches = 0
    mismatches = 0
    insertions = 0
    deletions = 0
    
    # Parse CIGAR string with numbers
    current_number = ''
    for char in cigar:
        if char.isdigit():
            current_number += char
        else:
            count = int(current_number)
            if char == '=':
                matches += count
            elif char == 'X':
                mismatches += count
            elif char == 'I':
                insertions += count
            elif char == 'D':
                deletions += count
            current_number = ''
    
    stats = {
        'len_unpolished': len(seq_unpolished),
        'len_polished': len(seq_polished),
        'edit_distance': alignment['editDistance'],
        'matches': matches,
        'mismatches': mismatches,
        'insertions': insertions, 
        'deletions': deletions
    }
    
    logging.info(f"Racon stats: {stats}")
    return stats