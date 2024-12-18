# opsin_analysis/plotting.py
from opsinth.utils import *
import matplotlib.pyplot as plt

def plot_coverage(results, output_prefix):
    logging.info("Starting coverage plot generation")
    
    coverages = results['anchors']['anchor_positions']
    coverages['no_anchor'] = {}
    coverages['both_anchors'] = {}
    roi = results['roi']
    ref_length = roi[0][2] - roi[0][1]
    logging.debug(f"Processing coverage data for ROI: {roi}")

    for anchor in coverages:
        coverages[anchor]['cov'] = [0]*ref_length

    for read in results['reads_aligned']:
        # Get start and end coords of the first alignment
        if results['reads_aligned'][read]['strand'] == "+":
            start = results['reads_aligned'][read]['aln']['locations'][0][0]  # Take only the first alignment
            end = results['reads_aligned'][read]['aln']['locations'][0][1]
        else:
            start = ref_length - results['reads_aligned'][read]['aln']['locations'][0][1]  # Take only the first alignment
            end = ref_length - results['reads_aligned'][read]['aln']['locations'][0][0]
        
        # Add coverage to the appropriate category
        if read in results['no_anchor_reads']:
            coverages['no_anchor']['cov'][start:end] = [x + 1 for x in coverages['no_anchor']['cov'][start:end]]
        elif read in results['double_anchor_reads']:
            coverages['both_anchors']['cov'][start:end] = [x + 1 for x in coverages['both_anchors']['cov'][start:end]]
        else:
            anchor = results['unique_anchor_alignments'][read]
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

    anchors = {anchor: {} for anchor in results['anchors']['anchor_positions'].keys()}
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
            anchor = results['unique_anchor_alignments'][read]
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
