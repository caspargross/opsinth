# opsin_analysis/plotting.py
import matplotlib.pyplot as plt

def plot_coverage(results, output_dir):
    print("Plotting coverages")
    coverages = results['anchors']['as_ref']
    roi = results['roi']

    for anchor in coverages:
        coverages[anchor]['cov'] = [0]*(roi[0][2] - roi[0][1])

    for read in results['reads_aligned']:
        unique_anchors = results['unique_anchor_alignments'][read]
        if results['reads_aligned'][read]['strand'] == "+" :
            start =  coverages[unique_anchors]['start']
            end = start + results['reads_aligned'][read]['ref_length']
        else:
            end = coverages[unique_anchors]['end']
            start = end - results['reads_aligned'][read]['ref_length']    
        
        #print(coverages[unique_anchor_alignments[read]]['cov'])
        for i,n in enumerate(coverages[unique_anchors]['cov']):
            coverages[unique_anchors]['cov'][i] = n+1 if start <= i <= end else n


    x = range(results['roi'][0][1], results['roi'][0][2])
    values = []
    labels = []

    for i, anchor in enumerate(coverages.keys()):
        # Assuming 'cov' is calculated elsewhere and added to results
        values.append(coverages[anchor]['cov'])
        labels.append(anchor)

    plt.stackplot(x, values, step='pre', labels=labels)
    plt.legend()
    plt.savefig(f"{output_dir}coverage_plot.png")
    #plt.show()

def plot_alignment_quality(results, output_dir):
    edit_distances = []
    matched_query_lengths = []

    for read in results['reads_aligned']:
        edit_distances.append(results['reads_aligned'][read]['aln']['editDistance'])
        matched_query_lengths.append(len(results['reads_aligned'][read]['seq']))

    plt.scatter(matched_query_lengths, edit_distances, edgecolors='black')
    plt.xlabel("Aligned read length (bp)")
    plt.ylabel("Edit distance to reference genome")
    plt.title("Edit Distance vs Aligned Read Length")
    plt.savefig(f"{output_dir}alignment_quality_plot.png")
    #plt.show()
