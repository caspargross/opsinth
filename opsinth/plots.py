# opsin_analysis/plotting.py
from opsinth.utils import *
import matplotlib.pyplot as plt
import numpy as np
from opsinth.resources import *


def plot_coverage(results, output_prefix):
    logging.info("Starting coverage plot generation")
    
    coverages = results['anchors_on_ref']['anchor_positions']
    coverages['no_anchor'] = {}
    coverages['both_anchors'] = {}
    roi = results['roi']
    ref_length = roi[0][2] - roi[0][1]
    logging.debug(f"Processing coverage data for ROI: {roi}")

    for anchor in coverages:
        coverages[anchor]['cov'] = [0]*ref_length

    for read in results['reads_aligned']:
        # Get start and end coords of the first alignment
        start = results['reads_aligned'][read]['aln']['locations'][0][0]  # Take only the first alignment
        end = results['reads_aligned'][read]['aln']['locations'][0][1]

        # Add coverage to the appropriate category
        if read in results['no_anchor_reads']:
            coverages['no_anchor']['cov'][start:end] = [x + 1 for x in coverages['no_anchor']['cov'][start:end]]
        elif read in results['double_anchor_reads']:
            coverages['both_anchors']['cov'][start:end] = [x + 1 for x in coverages['both_anchors']['cov'][start:end]]
        else:
            anchor = results['unique_anchor_reads'][read]
            coverages[anchor]['cov'][start:end] = [x + 1 for x in coverages[anchor]['cov'][start:end]]

    x = range(results['roi'][0][1], results['roi'][0][2])
    values = []
    labels = []

    for i, type in enumerate(coverages.keys()):
        values.append(coverages[type]['cov'])
        labels.append(type)

    plt.stackplot(x, values, step='pre', labels=labels)
    plt.legend()
    plt.savefig(f"{output_prefix}.coverage_plot.png")
    plt.close()

    logging.info("Coverage plot saved successfully")

def plot_alignment_quality(results, output_prefix):
    edit_distances = []
    matched_query_lengths = []
    colors = []

    anchors = {anchor: {} for anchor in results['anchors_on_ref']['anchor_positions'].keys()}
    anchors['no_anchor'] = {}
    anchors['both_anchors'] = {}
    anchor_color_map = {anchor: plt.cm.tab10(i) for i, anchor in enumerate(anchors)}

    for read in results['reads_aligned']:
        edit_distances.append(results['reads_aligned'][read]['aln']['editDistance'])
        matched_query_lengths.append(len(results['reads_aligned'][read]['seq']))
        
        if read in results['no_anchor_reads']:
            anchor = 'no_anchor'
        elif read in results['double_anchor_reads']:
            anchor = 'both_anchors' 
        else:
            anchor = results['unique_anchor_reads'][read]
        colors.append(anchor_color_map.get(anchor, 'gray'))

    plt.scatter(matched_query_lengths, edit_distances, edgecolors='black', c=colors)

    handles = [plt.Line2D([0], [0], marker='o', color='w', label=anchor,
                          markerfacecolor=anchor_color_map[anchor]) for anchor in anchors]
    plt.legend(handles=handles, title="Anchors")

    plt.xlabel("Aligned read length (bp)")
    plt.ylabel("Edit distance to reference genome")
    plt.title("Edit Distance vs Aligned Read Length")
    plt.xlim(min(matched_query_lengths), max(matched_query_lengths))
    plt.savefig(f"{output_prefix}.alignment_quality_plot.png")
    plt.close()

def plot_polish_stats(polish_stats, output_prefix):
    """
    Create two plots showing polishing statistics across racon rounds.
    
    Args:
        results: Dictionary containing polish_stats from run_polish_denovo
        output_prefix: Prefix for output files
    """
    logging.info("Generating polish statistics plots")
    
    if not polish_stats:
        logging.warning("No polish statistics found to plot")
        return
        
    rounds = [round + 1 for round in polish_stats.keys()]
    
    # First plot: Sequence lengths and edit distance
    plt.figure(figsize=(10, 6))
    
    # Create primary y-axis for lengths
    ax1 = plt.gca()
    ax1.plot(rounds, [stats['len_unpolished'] for stats in polish_stats.values()], 
             'b-', label='Unpolished Length', marker='o')
    ax1.plot(rounds, [stats['len_polished'] for stats in polish_stats.values()], 
             'g-', label='Polished Length', marker='o')
    ax1.set_xlabel('Polish Round')
    ax1.set_ylabel('Sequence Length (bp)', color='k')
    
    # Create secondary y-axis for edit distance
    ax2 = ax1.twinx()
    ax2.plot(rounds, [stats['edit_distance'] for stats in polish_stats.values()], 
             'r-', label='Edit Distance', marker='o')
    ax2.set_ylabel('Edit Distance', color='r')
    
    # Combine legends
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right')
    
    plt.title('Sequence Lengths and Edit Distance by Polish Round')
    plt.tight_layout()
    plt.savefig(f"{output_prefix}.polish_stats_lengths.png")
    plt.close()
    
    # Second plot: Operation counts
    plt.figure(figsize=(10, 6))
    x = np.arange(len(rounds))
    width = 0.2
    
    # Extract operation counts
    matches = [stats['matches'] for stats in polish_stats.values()]
    mismatches = [stats['mismatches'] for stats in polish_stats.values()]
    insertions = [stats['insertions'] for stats in polish_stats.values()]
    deletions = [stats['deletions'] for stats in polish_stats.values()]
    
    # Create bars
    plt.bar(x - width*1.5, matches, width, label='Matches')
    plt.bar(x - width/2, mismatches, width, label='Mismatches')
    plt.bar(x + width/2, insertions, width, label='Insertions')
    plt.bar(x + width*1.5, deletions, width, label='Deletions')
    
    plt.xlabel('Polish Round')
    plt.ylabel('Count')
    plt.title('Alignment Errors Fixed')
    plt.xticks(x, rounds)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{output_prefix}.polish_stats_operations.png")
    plt.close()
    
    logging.info("Polish statistics plots saved successfully")
