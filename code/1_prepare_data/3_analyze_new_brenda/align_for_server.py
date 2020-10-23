

## this is a version of the same script that's in my jupyter notebook ##
## it is copied here to be able to run the job on the server using multiple cores##
## place ec_data_uid_orgs_domain.tsv and fasta_data_ec_uid_orgs_domain_2.tsv in the same folder as script ##
## place all augmented fasta files in a sub-folder called data ##

import os
from os.path import join, dirname, basename, exists, isdir
import multiprocessing

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import io




def pairwise_align(id1, id2, seq1, seq2):
    '''
    Perform pairwise alignment.
    '''
    records = '>%s\n%s\n>%s\n%s' % (id1, seq1, id2, seq2) #prepare 'virtual' FASTA file
    records_handle = io.StringIO(records) #turn string into a handle
    tempdata = records_handle.getvalue()
    muscle_cline = MuscleCommandline(seqtype='protein')
    stdout, stderr = muscle_cline(stdin=tempdata)

    with io.StringIO(stdout) as fasta:
        aln = SeqIO.parse(fasta, "fasta")

        output = []
        for entry in aln:
            header = entry.description
            seq = entry.seq

            output.append(header)
            output.append(seq)

    return output




def pair_ident(Seq1, Seq2, single_gaps=True):
    '''
    Takes two aligned sequences and returns their percent identity.
    Assumes that Seq1 and Seq2 are sequence strings
    '''

    l=0.0 # counts alignment length, excluding identical gaps, but including single gaps
    n=0.0 # count number of single gaps
    i=0.0 # counts identity hits

    for j in range(len(Seq2)):
        if Seq1[j] == '-' and Seq2[j] == '-': #DON'T count identical gaps towards alignment length
            pass
        else:
            if Seq2[j] == Seq1[j]:
                i += 1 #count matches
            elif Seq2[j] == '-' or Seq1[j] == '-':
                n += 1 #count number of single gaps
            l += 1 #count total length with single gaps

    if single_gaps is True: #include single gaps
        percent = round(100*(i/l),1) #calculate identity

    elif single_gaps is False: #exclude single gaps
        if n >= l:
            percent = 0.0
        else:
            percent = round(100*(i/(l-n)),1) #calculate identity

    return percent

def get_ec_uid_and_seq(ec):
    '''
    Return nested dictionary with EC number primary keys and uid values.
    In the second level each uid holds its protein sequence.
    Only return for uniprot identifiers from a specific ec.
    '''
    data = {}

    filepath = 'data'
    fasta_filepath = join(filepath, '%s_clustered_sequences_90_augmented.fasta' % ec)

    # go through each entry in the fasta
    fasta = SeqIO.parse(fasta_filepath, "fasta")
    for record in fasta:
        header = record.description
        uid = header.split(';')[0].rstrip().upper()
        seq = str(record.seq)

        data[uid] = seq

    return data


def load_final_ktup_flatfile(infile):
    '''
    For loading data from each ec number, in case computation was interrupted in previous runs.
    '''
    data = {}
    with open(infile, 'r') as f:
        f.readline()
        for line in f:
            ec, org, uid, domain, best_match, identity = line.strip().split('\t')

            if data.get(ec) is None:
                data[ec] = {}

            if data[ec].get(org) is None:
                data[ec][org] = {'closest_match':[], 'ktuple_dist':[], 'uids':[], 'domain':domain}

            if uid != 'None':
                uid = uid.upper()

            if best_match != 'None':
                best_match = best_match.upper()

            data[ec][org]['uids'].append(uid)
            data[ec][org]['closest_match'].append(best_match)
            data[ec][org]['ktuple_dist'].append(identity)

    return data



def worker(package):
    '''
    '''
    index, org, id1, seq1, id2, seq2 = package

    # sometimes there is no closest match (i.e. no characterized sequence)
    if id2 == 'None' or seq2 == '':
        return (index, org, 0)
    # sequences that are too long make muscle fail, so ignore those
    elif len(seq1) > 10000 or len(seq2) > 10000:
        return (index, org, 0)

    id1 = '>%s'  %id1
    id2 = '>%s' % id2

    if seq1 == seq2:
        identity = 100.0
    else:
        output = pairwise_align(id1=id1, id2=id2, seq1=seq1, seq2=seq2)

        aln_seq1 = output[1]
        aln_seq2 = output[3]

        identity = pair_ident(Seq1=aln_seq1, Seq2=aln_seq2, single_gaps=True)

    return (index, org, identity)




def align_all(ecnum, reverse=False):
    '''
    '''

    data = load_final_ktup_flatfile(join('parsed_info', 'fasta_data_ec_uid_orgs_domain_2_kdist.tsv'))

    num_cores = multiprocessing.cpu_count()

    for ec in sorted(data.keys(), reverse=reverse):

        tempfile = join('aln_temp_data', '%s.tsv' % ec)

        if exists(tempfile):
            continue

        if not ec.startswith(ecnum):
            continue

        print(ec)

        seq_data = get_ec_uid_and_seq(ec)

        packages = []
        for org in sorted(data[ec].keys()):

            for i in range(len(data[ec][org]['uids'])):
                id1 = data[ec][org]['uids'][i]
                id2 = data[ec][org]['closest_match'][i]


                # sometimes there is no closest match (i.e. no characterized sequence), the worker will return an identity of 0
                if id2 == 'None':
                    seq1 = ''
                    seq2 = ''

                else:
                    seq1 = seq_data[id1]
                    seq2 = seq_data[id2]

                packages.append((i, org, id1, seq1, id2, seq2))

        with multiprocessing.Pool(processes=num_cores) as pool:
            results = pool.map(worker, packages)

        # build empty list for each org
        for res in results:
            index, org, ident = res

            if data[ec][org].get('aln_ident') is None:
                data[ec][org]['aln_ident'] = [None for s in data[ec][org]['uids']]

            # add in identity data
            data[ec][org]['aln_ident'][index] = ident


        # save to file
        write_aln_temp_flatfile(data[ec], tempfile)



def write_aln_temp_flatfile(data, outfile):
    '''
    For saving data from each ec number, in case computation does not finish.
    '''
    out_data = ['organism\tuid\tdomain\tbest_match\tktuple_dist\tidentity']
    for org in sorted(data.keys()):
        domain = data[org]['domain']

        for i in range(len(data[org]['uids'])):
            uid = data[org]['uids'][i]
            best_match = data[org]['closest_match'][i]
            ktup = data[org]['ktuple_dist'][i]
            identity = data[org]['aln_ident'][i]
            out_data.append('%s\t%s\t%s\t%s\t%s\t%s' % (org, uid, domain, best_match, ktup, identity))

    with open(outfile, 'w') as f:
        f.write('\n'.join(out_data))



def load_aln_temp_flatfile(infile):
    '''
    For loading data from each ec number, in case computation was interrupted in previous runs.
    '''
    print('Loading %s' % infile)
    data = {}
    with open(infile, 'r') as f:
        f.readline()
        for line in f:
            org, uid, domain, best_match, ktup, identity = line.strip().split('\t')

            if data.get(org) is None:
                data[org] = {'closest_match':[], 'ktuple_dist':[], 'uids':[], 'domain':domain, 'aln_ident':[]}

            data[org]['uids'].append(uid)
            data[org]['closest_match'].append(best_match)
            data[org]['ktuple_dist'].append(ktup)
            data[org]['aln_ident'].append(identity)

    return data



def write_final_aln_flatfile(data, outfile):
    '''
    '''
    out_data = ['ec\torganism\tuid\tdomain\tbest_match\tktuple_dist\tidentity']
    for ec in sorted(data.keys()):
        for org in sorted(data[ec].keys()):
            domain = data[ec][org]['domain']

            for i in range(len(data[ec][org]['uids'])):
                uid = data[ec][org]['uids'][i]
                best_match = data[ec][org]['closest_match'][i]
                ktup = data[ec][org]['ktuple_dist'][i]
                identity = data[ec][org]['aln_ident'][i]
                out_data.append('%s\t%s\t%s\t%s\t%s\t%s\t%s' % (ec, org, uid, domain, best_match, ktup, identity))

    with open(outfile, 'w') as f:
        f.write('\n'.join(out_data))







if __name__ == '__main__':

    import sys

    # align all the matches
    align_all(ecnum=sys.argv[1], reverse=True)

    # assemble all the individual files to one outfile
    outfile = join('parsed_info', 'fasta_data_ec_uid_orgs_domain_2_kdist_aln.tsv')

    if not exists(outfile):
        aln_data = {}
        folder_path = 'aln_temp_data'
        for fi in sorted(os.listdir(folder_path)):
            ec = '.'.join(fi.split('.')[:-1])
            aln_data[ec] = load_aln_temp_flatfile(join(folder_path, fi))

        # make a single outfile
        write_final_aln_flatfile(aln_data, outfile)

    print('done')
