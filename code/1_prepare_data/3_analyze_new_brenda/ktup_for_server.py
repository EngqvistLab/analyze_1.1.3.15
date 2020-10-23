
import os
from os.path import join, dirname, basename, exists, isdir

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import io
#
# from subprocess import Popen, PIPE
# import random
import numpy as np
# import pandas as pd
# from matplotlib import pyplot as plt
# from scipy.stats import spearmanr
import multiprocessing



# https://doi.org/10.1186/s13059-017-1319-7
# https://github.com/aziele/alfpy
# http://www.combio.pl/alfree
# pip install alfpy

from alfpy.utils import seqrecords
from alfpy import word_distance
from alfpy import word_pattern
from alfpy import word_vector
from alfpy.utils import distmatrix
from alfpy import ncd



def make_flatfile(data, header, outfile, has_uid):
    '''
    Generic function for saving the flatfiles
    '''
    out_data = [header]
    for ec in sorted(data.keys()):
        for org in sorted(data[ec].keys()):
            domain = data[ec][org]['domain']

            if domain is None:
                if has_uid is False:
                    out_data.append('%s\t%s' % (ec, org))

                elif has_uid is True:
                    for uid in data[ec][org]['uids']:
                            out_data.append('%s\t%s\t%s' % (ec, org, uid))

            else:
                if has_uid is False:
                    out_data.append('%s\t%s\t%s' % (ec, org, domain))

                elif has_uid is True:
                    for uid in data[ec][org]['uids']:
                            out_data.append('%s\t%s\t%s\t%s' % (ec, org, uid, domain))

    with open(outfile, 'w') as f:
        f.write('\n'.join(out_data))



def load_flatfile(infile, has_uid):
    '''
    Generic function for loading the flatfiles
    '''
    data = {}
    with open(infile, 'r') as f:
        f.readline()
        for line in f:

            if has_uid is True:
                if len(line.strip().split('\t')) == 4:
                    ec, org, uid, domain = line.strip().split('\t')
                else:
                    ec, org, uid = line.strip().split('\t')
                    domain = None

                if data.get(ec) is None:
                    data[ec] = {}

                if data[ec].get(org) is None:
                    data[ec][org] = {'uids':set([]), 'domain':None}

                data[ec][org]['uids'].add(uid.upper())
                data[ec][org]['domain'] = domain


            elif has_uid is False:
                if len(line.strip().split('\t')) == 3:
                    ec, org, domain = line.strip().split('\t')
                else:
                    ec, org = line.strip().split('\t')
                    domain = None

                if data.get(ec) is None:
                    data[ec] = {}

                if data[ec].get(org) is None:
                    data[ec][org] = {'uids':None, 'domain':None}

                data[ec][org]['domain'] = domain

            else:
                raise ValueError


    return data


def load_swissprot(data={}):
    '''
    Add in the SwissProt data regarding which sequences
    count as characterized.
    '''
    with open(join('parsed_info', 'SwissProt-2020_02-protein-evidence.tsv'), 'r') as f:
        f.readline()

        for line in f:
            uid, ec, org, orgid = line.strip().split('\t')
            domain = None

            # standardize orgname
            org = ' '.join(org.replace('.', '').split()[:2]).capitalize()

            if data.get(ec) is None:
                data[ec] = {}

            if data[ec].get(org) is None:
                data[ec][org] = {'uids':set([]), 'domain':None}

            data[ec][org]['uids'].add(uid.upper())
            data[ec][org]['domain'] = domain

    return data


def get_ec_uid_and_seq(ec):
    '''
    Return nested dictionary with EC number primary keys and uid values.
    In the second level each uid holds its protein sequence.
    Only return for uniprot identifiers from a specific ec.
    '''
    data = {}

    filepath ='data'
    fasta_filepath = join(filepath, '%s_clustered_sequences_90_augmented.fasta' % ec)

    # go through each entry in the fasta
    fasta = SeqIO.parse(fasta_filepath, "fasta")
    for record in fasta:
        header = record.description
        uid = header.split(';')[0].rstrip().upper()
        seq = record.seq

        data[uid] = seq

    return data



def ktup(ids, seqs):
    '''
    Perform the actual ktuple calculation.

    Moving these bits here really helps with keeping memory consumption down.
    Reason is unclear but I assume somehow Python is better able to claim
    back memory after each computation.
    '''
    seq_records = seqrecords.SeqRecords(id_list=ids, seq_list=seqs)
    p = word_pattern.create(seq_records.seq_list, word_size=3)

    counts = word_vector.Counts(seq_records.length_list, p)
    dist = word_distance.Distance(counts, 'google')
    matrix = distmatrix.create(seq_records.id_list, dist)

    # get the last column (which contains data for the last added identifier)
    return matrix.data[0, 1]



def calculate_ktup_identities(package):
    '''
    Loop through all the EC classes and compare characterized sequences with uncharacterized ones.
    '''
    ec, data, tempfile = package

    # get all the sequences in clustered fasta
    seq_data = get_ec_uid_and_seq(ec)

    # assemble a list of the characterized sequences
    char_ids = []
    char_seqs = []

    if ec_org_uid_data.get(ec) is None:
        pass
    else:
        for char_org in sorted(ec_org_uid_data[ec].keys()):
            for char_uid in ec_org_uid_data[ec][char_org]['uids']:
                # skip records where sequence is missing, or if it's too long
                if 0 < len(seq_data[char_uid]) <= 10000:
                    char_ids.append(char_uid)
                    char_seqs.append(str(seq_data[char_uid]))

    print('%s: %s characterized, %s sequences total' % (ec, len(char_ids), len(list(seq_data.keys()))))

#     print(char_ids)
#     print(char_seqs)

    # now go through the uncharacterized sequences one by one and get the closest characterized one
    counter = 0
    for organism in sorted(data.keys()):

        scores = []
        matches = []

        for uid in data[organism]['uids']:

            if counter % 500 == 0:
                print(ec, counter)
            counter += 1

            closest_match = 'None'
            best_score = 1.0

            # if none are characterized put 1 without retreiving sequence
            if char_ids == []:
                matches.append(closest_match)
                scores.append('%.3f' % best_score)
                continue

            # get sequence for this uid
            if seq_data.get(uid) is None:
                print(uid, ec)
            seq = seq_data[uid]

            # compare in turn with each of the characterized ones
            for i in range(0, len(char_ids)):
                tmp_ids = [char_ids[i]] + [uid]
                tmp_seqs = [char_seqs[i]] + [str(seq)]

                # make k-tupe calculation
                ktup_dist = ktup(tmp_ids, tmp_seqs)

                if ktup_dist < best_score:
                    best_score = ktup_dist
                    closest_match = char_ids[i]

            matches.append(closest_match)
            scores.append('%.3f' % best_score)

        data[organism]['closest_match'] = matches[:]
        data[organism]['ktuple_dist'] = scores[:]


    # now save the data to disk in case I have to re-run the computation
    write_temp_flatfile(data, tempfile)
    print('%s results written to file'  % ec)
    return None


def count_many(data, ecnum, reverse=False):
    '''
    Takes a list of many filepaths (full os path to file) and uses paralell processing to count frequencies in these.
    Returns a list of tuples with (filepath, occurance_dictionary).
    WARNING: The order of the files in the result may not be the same as in the input.
    '''
    num_cores = 64 #multiprocessing.cpu_count()

    # split up det data on a per ec basis
    packages = []
    for ec in sorted(data.keys(), reverse=reverse):

        tempfile = join('ktup_temp_data', '%s.tsv' % ec)

        if exists(tempfile):
            continue

        # if not ec.startswith(ecnum):
        #    continue

        packages.append((ec, data[ec], tempfile))

    with multiprocessing.Pool(processes=num_cores) as pool:
        results = pool.map(calculate_ktup_identities, packages)

    for res in results:
        pass
    return None


def write_temp_flatfile(data, outfile):
    '''
    For saving data from each ec number, in case computation does not finish.
    '''
    out_data = ['organism\tuid\tdomain\tbest_match\tktuple_dist']
    for org in sorted(data.keys()):
        domain = data[org]['domain']

        for i in range(len(data[org]['uids'])):
            uid = list(data[org]['uids'])[i]
            best_match = data[org]['closest_match'][i]
            ktup = data[org]['ktuple_dist'][i]
            out_data.append('%s\t%s\t%s\t%s\t%s' % (org, uid, domain, best_match, ktup))

    with open(outfile, 'w') as f:
        f.write('\n'.join(out_data))



def load_temp_flatfile(infile):
    '''
    For loading data from each ec number, in case computation was interrupted in previous runs.
    '''
    print('Loading %s' % infile)
    data = {}
    with open(infile, 'r') as f:
        f.readline()
        for line in f:
            org, uid, domain, best_match, ktup = line.strip().split('\t')

            if data.get(org) is None:
                data[org] = {'closest_match':[], 'ktuple_dist':[], 'uids':[], 'domain':domain}

            data[org]['uids'].append(uid)
            data[org]['closest_match'].append(best_match)
            data[org]['ktuple_dist'].append(ktup)

    return data




def write_final_ktup_flatfile(data, outfile):
    '''
    '''
    out_data = ['ec\torganism\tuid\tdomain\tbest_match\tktuple_dist']
    for ec in sorted(data.keys()):
        for org in sorted(data[ec].keys()):
            domain = data[ec][org]['domain']

            for i in range(len(data[ec][org]['uids'])):
                uid = data[ec][org]['uids'][i]
                best_match = data[ec][org]['closest_match'][i]
                identity = data[ec][org]['ktuple_dist'][i]
                out_data.append('%s\t%s\t%s\t%s\t%s\t%s' % (ec, org, uid, domain, best_match, identity))

    with open(outfile, 'w') as f:
        f.write('\n'.join(out_data))




if __name__ == '__main__':

    import sys

    # get the fasta data
    infile_fasta = join('parsed_info', 'fasta_data_ec_uid_orgs_domain_2.tsv')
    fasta_data = load_flatfile(infile_fasta, has_uid=True)

    # get information on the characterized from BRENDA
    infile_html = join('parsed_info', 'ec_data_uid_orgs_domain_2.tsv')
    ec_org_uid_data = load_flatfile(infile_html, has_uid=True)

    # get information on the characterized from SwissProt
    ec_org_uid_data = load_swissprot(ec_org_uid_data)

    # calculate the k-tuple measures
    count_many(fasta_data, ecnum=sys.argv[1], reverse=True)


    # load each of the individual files
    outfile = join('parsed_info', 'fasta_data_ec_uid_orgs_domain_2_kdist.tsv')

    if not exists(outfile):
        ktup_data = {}
        folder_path = 'ktup_temp_data'
        for fi in sorted(os.listdir(folder_path)):
            ec = '.'.join(fi.split('.')[:-1])
            ktup_data[ec] = load_temp_flatfile(join(folder_path, fi))

        # make a single outfile
        write_final_ktup_flatfile(ktup_data, outfile)

    print('done')
