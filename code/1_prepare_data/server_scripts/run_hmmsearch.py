import os
import subprocess
from os.path import join, exists, getsize
import multiprocessing
import pandas as pd




def worker(package):
    '''
    '''
    fi, inpath, outpath, pfam_hmms = package

    fasta_file = join(inpath, fi)
    hmm_outfile = join(outpath, fi+'.out')

    if not exists(hmm_outfile):
        mycmd = 'hmmsearch -E 1e-15 --domtblout %s %s %s > tmp.out' % (hmm_outfile, pfam_hmms, fasta_file)
        os.system(mycmd)




def search_all(ecnum, reverse=False):
    '''
    '''
    num_cores = multiprocessing.cpu_count()

    # define folders
    outpath = join('.', 'brenda_domains')
    inpath = join('.', 'data')
    pfam_hmms = join(outpath, 'Pfam-A.hmm')

    # run assemble a list of the file to run
    packages = []
    for fi in sorted(os.listdir(inpath), reverse=reverse):

        if not fi.startswith(ecnum):
            continue

        if not fi.endswith('.fasta'):
            continue

        packages.append((fi, inpath, outpath, pfam_hmms))

    # distribute work on the cores
    with multiprocessing.Pool(processes=num_cores) as pool:
        results = pool.map(worker, packages)





def parse_hmm_output(filepath, ec):
    '''
    Parse an hmm output file and return as pandas data frame.
    '''
    data = {'ec':[], 'uid':[], 'hmm_model':[], 'pfam':[], 'hmm_model_len':[], 'hmm_match_from':[],
            'hmm_match_to':[], 'hmm_match_coverage':[], 'match_evalue':[], 'gene_match_from':[],
            'gene_match_to':[]}

    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split()
            if len(parts) < 10:
                print(line)

            uid = parts[0].split(';')[0]

            dom_name = parts[3]
            pfam = parts[4]

            evalue = parts[6]

            qlen = parts[5]
            hmm_from = parts[15]
            hmm_to = parts[16]
            coverage = (int(hmm_to)-int(hmm_from))/int(qlen)

            ali_from = parts[17]
            ali_to = parts[18]

            if coverage >= 0.35:
                data['ec'].append(ec)
                data['uid'].append(uid)
                data['hmm_model'].append(dom_name)
                data['pfam'].append(pfam)
                data['hmm_model_len'].append(qlen)
                data['hmm_match_from'].append(hmm_from)
                data['hmm_match_to'].append(hmm_to)
                data['hmm_match_coverage'].append('%.2f' % float(coverage))
                data['match_evalue'].append(evalue)
                data['gene_match_from'].append(ali_from)
                data['gene_match_to'].append(ali_to)

    return pd.DataFrame(data)


if __name__ == '__main__':

    import sys

    # align all the matches
    search_all(ecnum=sys.argv[1], reverse=False)

    # parse the outfiles
    filepath = join('.', 'brenda_domains')
    all_genome_data = None
    for fi in sorted(os.listdir(filepath)):

        if not fi.endswith('.fasta.out'):
            continue

        if getsize(join(filepath, fi)) == 0:
            continue

        print(fi)

        # parse outputs
        ec = fi.split('_')[0]
        hmm_data = parse_hmm_output(join(filepath, fi), ec)

        # add to the main data
        if all_genome_data is None:
            all_genome_data = hmm_data
        else:
            all_genome_data = all_genome_data.append(hmm_data)

    all_genome_data.to_csv(join('.', 'pfam_hmm_results.tsv'), sep='\t', index=False)

    print('done')
