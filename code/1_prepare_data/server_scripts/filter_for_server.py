import pandas as pd
from os.path import join, exists
import numpy as np

def test_overlap(pos1, pos2, cutoff):
    '''
    Find out whether two defined regions overlap or not.
    Overlap is defined True when at least one of the domains
    overlap with the other above a certain fraction of the whole.
    '''
    pos1 = [int(s) for s in pos1]
    pos2 = [int(s) for s in pos2]
    both_pairs = sorted([pos1, pos2])

    # simplest case, they do not overlap at all
    if both_pairs[0][1] < both_pairs[1][0]:
        return False

    else:
        # get overlap part
        overlap_region = both_pairs[0][1] - both_pairs[1][0]

        # fraction of first domain
        frac_first = overlap_region / (both_pairs[0][1] - both_pairs[0][0])

        # fraction of second domain
        frac_second = overlap_region / (both_pairs[1][1] - both_pairs[1][0])

        if max(frac_first, frac_second) > cutoff:
            #print(both_pairs, True, max(frac_first, frac_second))
            return True
        else:
            #print(both_pairs, False, max(frac_first, frac_second))
            return False




def remove_overlaps(data, cutoff):
    '''
    When domains overlap (as defined by the cutoff level),
    keep only the domain with the best e-value score.
    '''
    overlaps = [True]
    while any(overlaps):

        start_end_pairs = list(zip(data.gene_match_from.values, data.gene_match_to.values))
        e_vals = data.match_evalue.values
        e_vals = e_vals.astype(np.float)


        for i in range(len(start_end_pairs)):
            overlaps = []
            for j in range(len(start_end_pairs[i+1:])):

                overlaps.append(test_overlap(start_end_pairs[i], start_end_pairs[i+1+j], cutoff))

            # prepend true or false for the first sequence, depending on whether there was an overlap
            if any(overlaps):
                overlaps = [False for s in range(i)] + [True] + overlaps

                # from the ones with an overlap, what is the position of the one with the lowest e-value
                overlap_e_vals = e_vals[overlaps]

                idx = overlap_e_vals.argmin()
                #print(overlaps)
                #print(overlap_e_vals)
                #print(i+idx)


                # flip the one that said true but had the smallest e-value
                overlaps[i+idx] = False

                # invert the list to flip True vs False
                overlaps = [not i for i in overlaps]

                # filter the data to keep all the ones that originally said False, but now say true
                data = data[overlaps]
                #print(data)
                break

    return data



outfile = 'pfam_hmm_results.tsv'
df = pd.read_csv(outfile, sep='\t')


# test a range of cutoff values to figure out which one is best
cutoff_vals = [s/100 for s in range(5, 105, 5)]
motifs = []
for cutoff in cutoff_vals:
    print(cutoff)

    temp_data_filepath = 'motifs_result_{}.tsv'.format(cutoff)
    if not exists(temp_data_filepath):
        non_redundant_data = df.groupby('uid').apply(lambda x: remove_overlaps(x, cutoff=cutoff))
    else:
        non_redundant_data = pd.read_csv(temp_data_filepath, sep='\t')

    motifs.append(non_redundant_data.uid.count())


print('Unfiltered domain number: %s' % df.uid.count())
