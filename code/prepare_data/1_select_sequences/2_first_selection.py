#!/usr/bin/env python3
"""
This is the script that was used to make a sequence selection from the 1.1.3.n sequences

Copyright (C) 2017  Martin Engqvist Lab
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

__version__ = '0.0.1' # MAJOR.MINOR.PATCH http://semver.org/
__author__ = 'Martin Engqvist'

import os
from dotenv import load_dotenv, find_dotenv
from os.path import join, dirname, isdir

### Load environmental variables from the project root directory ###
# find .env automagically by walking up directories until it's found
dotenv_path = find_dotenv()

# load up the entries as environment variables
load_dotenv(dotenv_path)

# now you can get the variables using their names, for example:
# database_url = os.environ.get("DATABASE_URL")
# other_variable = os.environ.get("OTHER_VARIABLE")

# get the current working directory
CURRENT_DIR = os.getcwd()

# project root directory
PROJ_ROOT_DIR = dirname(dotenv_path)

# internal data raw directory
INT_RAW_DIR = join(PROJ_ROOT_DIR, 'data/raw_internal/')

# external data raw directory
EXT_RAW_DIR = join(PROJ_ROOT_DIR, 'data/raw_external/')

# intermediate data directory
INTERIM_DIR = join(PROJ_ROOT_DIR, 'data/intermediate/')

# final data directory
FINAL_DIR = join(PROJ_ROOT_DIR, 'data/final/')

# output directory
OUT_DIR = join(PROJ_ROOT_DIR, 'results/')

# figure output directory
FIG_DIR = join(PROJ_ROOT_DIR, 'results/figures/')

# picture output directory
PIC_DIR = join(PROJ_ROOT_DIR, 'results/pictures/')


#### Your code here ####

import fasta
import subprocess
import distance
import bioseq
import logging
import sys
import multiprocessing
import types

#setup logfile
logger = logging.getLogger('myapp')
hdlr = logging.FileHandler(join(PROJ_ROOT_DIR, 'data', 'clustering_and_sampling.log'))
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr)
logger.setLevel(logging.INFO)


def initial_info():
    '''Generate information on the script'''
    #get python version
    if sys.version_info[0] < 3:
        logger.error('Python2 detected at runtime. This script was written for Python3')
        raise "Script must be run using Python 3"
    else:
        logger.info('Python3 detected at runtime.')

    #get version of this script
    logger.info('Script version %s used for computation.' % __version__)

    #get command used to run script
    #here

    #get modules and their versions
    def imports():
        for name, val in globals().items():
            if isinstance(val, types.ModuleType) and val.__name__ != 'builtins':
                try:
                  version = val.__version__
                except:
                  version = 'unknown'
                yield val.__name__, version
    logger.info('Detecting loaded modules.\nModule   version\n----------------\n%s' % '\n'.join(['%s = %s' % (s, y) for s, y in list(imports())]) )



def setup_folders():
    '''Make sure the destination folders exist'''
    logging.info('Ensuring that the data output folders exist.')

    folder_name = 'brenda_2017_1'
    interim_folder_path = join(INTERIM_DIR, folder_name)
    final_folder_path = join(FINAL_DIR, folder_name)

    #interim destination folder
    logging.info(interim_folder_path)
    if isdir(interim_folder_path) is False:
        os.makedirs(interim_folder_path)

    #final data folder
    logging.info(final_folder_path)
    if isdir(final_folder_path) is False:
        os.makedirs(final_folder_path)



def make_fasta():
    '''Convert BENDA csv file to FASTA'''
    infile = join(EXT_RAW_DIR, 'brenda_2017_1', '1_1_3__BRENDA_sequences.csv')
    outfile = join(INTERIM_DIR, 'brenda_2017_1', '1_1_3__BRENDA_sequences_brenda_2017_1.fasta')
    logging.info('Converting BRENDA csv file "%s" to fasta file "%s"' % (infile, outfile))

    with open(infile, 'r') as f:
        f.readline()
        f.readline()
        in_file = f.readlines()

    with open(outfile, 'w') as f:
        for line in in_file:
            header = ';'.join(line.split('\t')[:-1])
            seq = line.split('\t')[-1].strip()
            f.write('>%s\n%s\n' % (header, seq))



def filter_fasta(min_len = 200, max_len = 580):
    '''Remove sequences that are too long or too short, also sequenes that have X in them'''
    infile = join(INTERIM_DIR, 'brenda_2017_1', '1_1_3__BRENDA_sequences_2017_1.fasta')
    outfile = join(INTERIM_DIR, 'brenda_2017_1', '1_1_3__BRENDA_sequences_filtered_2017_1.fasta')
    logging.info('Filtering fasta file "%s" to retain sequences between %s and %s AA in length. Also remove sequences with "X" in them. Output to "%s".' % (infile, min_len, max_len, outfile))

    retained = 0
    removed_len = 0
    removed_x = 0
    removed_m = 0
    with open(outfile, 'w') as f:
        for entry in fasta.parse_file(infile):
            header, seq = entry
            seq = seq.lower()
            if min_len <= len(seq) <= max_len and 'x' not in seq and seq[0] == 'm':
                retained += 1
                f.write('>%s\n%s\n' % (header, seq))
            elif 'x' in seq:
                removed_x += 1
            elif seq[0] != 'm':
                removed_m += 1
            else:
                removed_len += 1

    logger.info('%s sequences kept. %s removed because of length, %s removed because of X, %s removed because of no m' % (retained, removed_len, removed_x, removed_m))



def build_blast_db():
    '''Build blast database, requires NCBI BLAST (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/)'''
    logger.info('Building BLAST database.')
    logger.info('Version info for makeblastdb: %s' % str(subprocess.check_output(['makeblastdb', '-version']), 'utf-8'))

    #asseble command
    mycmd = ['makeblastdb', '-in', '../../data/intermediate/brenda_2017_1/1_1_3__BRENDA_sequences_filtered_2017_1.fasta', '-dbtype', 'prot', '-out', '../../data/intermediate/brenda_2017_1/my_prot_blast_db']
    logger.info(' '.join(mycmd))

    #build database
    subprocess.call(mycmd)



def blast_fasta_file():
    '''Blast the fasta file, all vs all, requires NCBI BLAST (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/). -outfmt 6 in BLAST+ is equivalent to legacy -m 8'''
    logger.info('Performing all vs all BLAST')
    logger.info('Version info for blastp: %s' % str(subprocess.check_output(['blastp', '-version']), 'utf-8'))

    #assemble command
    mycmd = ['blastp', '-db', '../../data/intermediate/brenda_2017_1/my_prot_blast_db', '-query', '../../data/intermediate/brenda_2017_1/1_1_3__BRENDA_sequences_filtered_2017_1.fasta', '-outfmt', '6', '-out', '../../data/intermediate/brenda_2017_1/all-vs-all.tsv', '-num_threads', '4']
    logger.info(' '.join(mycmd))

    #run blast
    subprocess.call(mycmd)



def blast_to_abc():
    '''Convert blast result to abc format, requires blast result to be in the correct format'''
    logger.info('Filtering BLAST output to make abc file (node, node, edge) using cut')
    logger.info('Version info for cut: %s' % str(subprocess.check_output(['cut', '--version']), 'utf-8'))

    #assemble command
    mycmd = 'cut -f 1,2,11 ../../data/intermediate/brenda_2017_1/all-vs-all.tsv > ../../data/intermediate/brenda_2017_1/all-vs-all.abc'
    logger.info(mycmd)

    #execute
    os.system(mycmd)



def build_network():
    '''Make network from abc file using mcxload'''
    logger.info('Build network using abc file')
    logger.info('Version info for mcxload: %s' % str(subprocess.check_output(['mcxload', '--version']), 'utf-8'))

    #assemble command
    mycmd = "mcxload -abc ../../data/intermediate/brenda_2017_1/all-vs-all.abc --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o ../../data/intermediate/brenda_2017_1/seq.mci -write-tab ../../data/intermediate/brenda_2017_1/seq.tab"
    logger.info(mycmd)

    #execute
    os.system(mycmd)



def cluster(granularity_levels = ('1.4', '2', '4', '6')):
    '''cluster sequences using mcl'''
    logger.info('Cluster sequences in network using mcl with the granularity levels %s' % ', '.join(["%s" % s for s in granularity_levels]))
    logger.info('Version info for mcl: %s' % str(subprocess.check_output(['mcl', '--version']), 'utf-8'))

    for level in granularity_levels:
        level_formated = int(float(level)*10)
        mycmd = 'mcl ../../data/intermediate/brenda_2017_1/seq.mci -I %s -o ../../data/intermediate/brenda_2017_1/out.seq.mci.I%s' % (level, level_formated)
        logger.info(mycmd)
        os.system(mycmd)

        mycmd = 'mcxdump -icl ../../data/intermediate/brenda_2017_1/out.seq.mci.I%s -tabr ../../data/intermediate/brenda_2017_1/seq.tab -o ../../data/final/brenda_2017_1/dump.seq.mci.I%s' % (level_formated, level_formated)
        logger.info(mycmd)
        os.system(mycmd)



def cluster_measures():
    '''Get measures of clustering'''
    logger.info('Comparing clusters')
    pass



def make_cluster_fasta_files():
    '''Use clustering data to make one fasta file for each cluster'''
    cluster_infile = '../../data/final/brenda_2017_1/dump.seq.mci.I14'
    sequence_infile = join(INTERIM_DIR, 'brenda_2017_1', '1_1_3__BRENDA_sequences_filtered_2017_1.fasta')
    logger.info('Parsing clusters from "%s" and sequences from "%s" to generate separate fasta files.' % (cluster_infile, sequence_infile))

    #get information on which sequence is in what cluster
    with open(cluster_infile, 'r') as f:
        in_data = f.read().split('\n')

    cluster_data = {}
    cl = 0
    for cluster in in_data:
        for entry in cluster.split():
            cluster_data[entry.split(';')[0]] = cl
        cl += 1
    logger.info('%s total clusters' % cl)

    #load fasta file and connect sequence info with the cluster number
    with open(sequence_infile, 'r') as f:
        fasta_data = f.read()

    out_data = {}
    for entry in fasta.parse_string(fasta_data):
        header, seq = entry
        code = header.split(';')[0]
        cluster = cluster_data[code]

        if out_data.get(cluster) is None:
            out_data[cluster] = []

        out_data[cluster].append(entry)

    #write data to file
    for key in out_data.keys():
        cluster_outfile = join(FINAL_DIR, 'brenda_2017_1', 'cluster_%s_sequences.fasta' % key)
        logger.info('Writing cluster sequences to file "%s".' % cluster_outfile)
        with open(cluster_outfile, 'w') as f:
            for entry in out_data[key]:
                header, seq = entry
                f.write('>%s\n%s\n' % (header, seq))



def align_clusters():
    '''Make multiple sequence alignment of each cluster fasta file'''
    logger.info('Making multiple sequence alignment of each cluster using MUSCLE')
    logger.info('Version info for muscle: %s' % str(subprocess.check_output(['muscle', '-version']), 'utf-8'))

    files = os.listdir(join(FINAL_DIR, 'brenda_2017_1'))
    for fi in sorted(files):
        if fi.endswith('sequences.fasta'):
            mycmd = 'muscle -in ../../data/final/brenda_2017_1/%s -out ../../data/final/brenda_2017_1/%s' % (fi, fi.replace('sequences.fasta', 'alignment.fasta'))
            logger.info(mycmd)
            os.system(mycmd)



def worker(args):
    '''helper function to make paralell processing possible'''
    headers, sequences, outfile = args
    distance.greedy(headers, sequences, outfile)



def compute_mi():
    '''find informative sequences for each cluster that captures x% of the information therein'''
    logger.info('Using greedy maximization algorithm to compute which sequences most contribute to mutual information.')
    folder = join(FINAL_DIR, 'brenda_2017_1')
    files = os.listdir(folder)

    #collect data
    data = []
    for fi in sorted(files):
        if fi.endswith('alignment.fasta'):
            headers = []
            sequences = []
            fasta_data = fasta.parse_file(join(folder, fi))
            for entry in fasta_data:
                header, seq = entry
                headers.append(header)
                sequences.append(seq)
            data.append( (headers, sequences, join(folder, fi.replace('alignment.fasta', 'mutual_information.tsv'))) )

    #make a process pool that can be used to distribute work on all cores
    num_cores = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(num_cores)
    pool.map(worker, data)



def validate_insert(insert):
    '''Make sure that the insert is ok, no stop codons etc.'''
    if type(insert) == str:
        insert = bioseq.DNA(insert).lower()
    else:
        insert = insert.lower()

    #make sure it has atg
    assert insert.startswith('atg')

    #make sure it is full-length
    assert len(insert.translate()) == len(insert)/3

    #make sure there is no stop codon
    assert '*' not in insert.translate()

    #make sure it ends with histidines
    assert insert.translate().endswith('H'*4)



def validate_construct(construct, insert_start_pos, target_prot):
    '''check the sequences that they are ok, translation etc.'''
    if type(construct) == str:
        construct = bioseq.DNA(construct)
    else:
        construct = construct.lower()
    #logger.info('Validating construct integrity.')

    #translate from the insert starting point
    prot = construct.translate(start=insert_start_pos)

    #make sure there is a stop
    assert '*' in prot, prot

    #make sure the stop is at the right place
    prot_product = prot.split('*')[0]
    assert prot_product.endswith('H'*10), prot

    #make sure insert is inside
    target_prot = target_prot.upper()
    assert prot_product.startswith(target_prot), (prot_product, target_prot)



def make_order_file(mi_cutoff = 85):
    '''Generate file with sequences to order'''
    logger.info('Generate output file with sequences to order, using sequences that explain %s%% of the variation in a cluster.' % mi_cutoff)
    folder = join(FINAL_DIR, 'brenda_2017_1')
    selected_ids_outfile = join(folder, 'selected_ids.tsv')
    constructs_outfile = join(folder, 'sequence_order.tsv')
    utr_outfile = join(folder, 'utr_info.tsv')
    info_outfile = join(folder, 'temperature_info.tsv')

    files = os.listdir(folder)
    num_sequences = 0
    selected_sequences = []
    all_sequences = []

    #for all files in folder
    selected_ids = []
    for fi in sorted(files):

        #look at the ones with mutual information data
        if fi.endswith('mutual_information.tsv'):
            with open(join(folder, fi), 'r') as f:

                #skip the header and the first empty line, then read the data lines
                f.readline()
                f.readline()
                for line in f:
                    num, header, mi, entr, seq = line.split('\t')

                    #keep track of all sequences, both selected and not
                    all_sequences.append((header, seq))

                    #if the cumulative information in the sequences so far is below the cutoff
                    if float(mi) <= mi_cutoff:

                        #open fasta file and get the sequences
                        fasta_data = fasta.parse_file(join(folder, fi.replace('mutual_information.tsv', 'sequences.fasta')))
                        for entry in fasta_data:
                            head, seq = entry
                            if head == header:
                                selected_sequences.append((header, seq))
                                selected_ids.append(header.split(';')[0]) #save id for future use
                                num_sequences += 1
                                break

    #save selected ids to file
    with open(selected_ids_outfile, 'w') as f:
        f.write('\n'.join(selected_ids))

    logger.info('total sequences at cutoff level: %s', num_sequences)

    #reverse translate them and combine with backbone
    complete_backbone =  'tggcgaatgggacgcgccctgtagcggcgcattaagcgcggcgggtgtggtggttacgcgcagcgtgaccgctacacttgccagcgccctagcgcccgctcctttcgctttcttcccttcctttctcgccacgttcgccggctttccccgtcaagctctaaatcgggggctccctttagggttccgatttagtgctttacggcacctcgaccccaaaaaacttgattagggtgatggttcacgtagtgggccatcgccctgatagacggtttttcgccctttgacgttggagtccacgttctttaatagtggactcttgttccaaactggaacaacactcaaccctatctcggtctattcttttgatttataagggattttgccgatttcggcctattggttaaaaaatgagctgatttaacaaaaatttaacgcgaattttaacaaaatattaacgtttacaatttcaggtggcacttttcggggaaatgtgcgcggaacccctatttgtttatttttctaaatacattcaaatatgtatccgctcatgagacaataaccctgataaatgcttcaataatattgaaaaaggaagagtatgagtattcaacatttccgtgtcgcccttattcccttttttgcggcattttgccttcctgtttttgctcacccagaaacgctggtgaaagtaaaagatgctgaagatcagttgggtgcacgagtgggttacatcgaactggatctcaacagcggtaagatccttgagagttttcgccccgaagaacgttttccaatgatgagcacttttaaagttctgctatgtggcgcggtattatcccgtattgacgccgggcaagagcaactcggtcgccgcatacactattctcagaatgacttggttgagtactcaccagtcacagaaaagcatcttacggatggcatgacagtaagagaattatgcagtgctgccataaccatgagtgataacactgcggccaacttacttctgacaacgatcggaggaccgaaggagctaaccgcttttttgcacaacatgggggatcatgtaactcgccttgatcgttgggaaccggagctgaatgaagccataccaaacgacgagcgtgacaccacgatgcctgcagcaatggcaacaacgttgcgcaaactattaactggcgaactacttactctagcttcccggcaacaattaatagactggatggaggcggataaagttgcaggaccacttctgcgctcggcccttccggctggctggtttattgctgataaatctggagccggtgagcgtgggtctcgcggtatcattgcagcactggggccagatggtaagccctcccgtatcgtagttatctacacgacggggagtcaggcaactatggatgaacgaaatagacagatcgctgagataggtgcctcactgattaagcattggtaactgtcagaccaagtttactcatatatactttagattgatttaaaacttcatttttaatttaaaaggatctaggtgaagatcctttttgataatctcatgaccaaaatcccttaacgtgagttttcgttccactgagcgtcagaccccgtagaaaagatcaaaggatcttcttgagatcctttttttctgcgcgtaatctgctgcttgcaaacaaaaaaaccaccgctaccagcggtggtttgtttgccggatcaagagctaccaactctttttccgaaggtaactggcttcagcagagcgcagataccaaatactgtccttctagtgtagccgtagttaggccaccacttcaagaactctgtagcaccgcctacatacctcgctctgctaatcctgttaccagtggctgctgccagtggcgataagtcgtgtcttaccgggttggactcaagacgatagttaccggataaggcgcagcggtcgggctgaacggggggttcgtgcacacagcccagcttggagcgaacgacctacaccgaactgagatacctacagcgtgagctatgagaaagcgccacgcttcccgaagggagaaaggcggacaggtatccggtaagcggcagggtcggaacaggagagcgcacgagggagcttccagggggaaacgcctggtatctttatagtcctgtcgggtttcgccacctctgacttgagcgtcgatttttgtgatgctcgtcaggggggcggagcctatggaaaaacgccagcaacgcggcctttttacggttcctggccttttgctggccttttgctcacatgttctttcctgcgttatcccctgattctgtggataaccgtattaccgcctttgagtgagctgataccgctcgccgcagccgaacgaccgagcgcagcgagtcagtgagcgaggaagcggaagagcgcctgatgcggtattttctccttacgcatctgtgcggtatttcacaccgcatatatggtgcactctcagtacaatctgctctgatgccgcatagttaagccagtatacactccgctatcgctacgtgactgggtcatggctgcgccccgacacccgccaacacccgctgacgcgccctgacgggcttgtctgctcccggcatccgcttacagacaagctgtgaccgtctccgggagctgcatgtgtcagaggttttcaccgtcatcaccgaaacgcgcgaggcagctgcggtaaagctcatcagcgtggtcgtgaagcgattcacagatgtctgcctgttcatccgcgtccagctcgttgagtttctccagaagcgttaatgtctggcttctgataaagcgggccatgttaagggcggttttttcctgtttggtcactgatgcctccgtgtaagggggatttctgttcatgggggtaatgataccgatgaaacgagagaggatgctcacgatacgggttactgatgatgaacatgcccggttactggaacgttgtgagggtaaacaactggcggtatggatgcggcgggaccagagaaaaatcactcagggtcaatgccagcgcttcgttaatacagatgtaggtgttccacagggtagccagcagcatcctgcgatgcagatccggaacataatggtgcagggcgctgacttccgcgtttccagactttacgaaacacggaaaccgaagaccattcatgttgttgctcaggtcgcagacgttttgcagcagcagtcgcttcacgttcgctcgcgtatcggtgattcattctgctaaccagtaaggcaaccccgccagcctagccgggtcctcaacgacaggagcacgatcatgcgcacccgtggggccgccatgccggcgataatggcctgcttctcgccgaaacgtttggtggcgggaccagtgacgaaggcttgagcgagggcgtgcaagattccgaataccgcaagcgacaggccgatcatcgtcgcgctccagcgaaagcggtcctcgccgaaaatgacccagagcgctgccggcacctgtcctacgagttgcatgataaagaagacagtcataagtgcggcgacgatagtcatgccccgcgcccaccggaaggagctgactgggttgaaggctctcaagggcatcggtcgagatcccggtgcctaatgagtgagctaacttacattaattgcgttgcgctcactgcccgctttccagtcgggaaacctgtcgtgccagctgcattaatgaatcggccaacgcgcggggagaggcggtttgcgtattgggcgccagggtggtttttcttttcaccagtgagacgggcaacagctgattgcccttcaccgcctggccctgagagagttgcagcaagcggtccacgctggtttgccccagcaggcgaaaatcctgtttgatggtggttaacggcgggatataacatgagctgtcttcggtatcgtcgtatcccactaccgagatatccgcaccaacgcgcagcccggactcggtaatggcgcgcattgcgcccagcgccatctgatcgttggcaaccagcatcgcagtgggaacgatgccctcattcagcatttgcatggtttgttgaaaaccggacatggcactccagtcgccttcccgttccgctatcggctgaatttgattgcgagtgagatatttatgccagccagccagacgcagacgcgccgagacagaacttaatgggcccgctaacagcgcgatttgctggtgacccaatgcgaccagatgctccacgcccagtcgcgtaccgtcttcatgggagaaaataatactgttgatgggtgtctggtcagagacatcaagaaataacgccggaacattagtgcaggcagcttccacagcaatggcatcctggtcatccagcggatagttaatgatcagcccactgacgcgttgcgcgagaagattgtgcaccgccgctttacaggcttcgacgccgcttcgttctaccatcgacaccaccacgctggcacccagttgatcggcgcgagatttaatcgccgcgacaatttgcgacggcgcgtgcagggccagactggaggtggcaacgccaatcagcaacgactgtttgcccgccagttgttgtgccacgcggttgggaatgtaattcagctccgccatcgccgcttccactttttcccgcgttttcgcagaaacgtggctggcctggttcaccacgcgggaaacggtctgataagagacaccggcatactctgcgacatcgtataacgttactggtttcacattcaccaccctgaattgactctcttccgggcgctatcatgccataccgcgaaaggttttgcgccattcgatggtgtccgggatctcgacgctctcccttatgcgactcctgcattaggaagcagcccagtagtaggttgaggccgttgagcaccgccgccgcaaggaatggtgcatgcaaggagatggcgcccaacagtcccccggccacggggcctgccaccatacccacgccgaaacaagcgctcatgagcccgaagtggcgagcccgatcttccccatcggtgatgtcggcgatataggcgccagcaaccgcacctgtggcgccggtgatgccggccacgatgcgtccggcgtagaggatcgagatctcgatcccgcgaaattaatacgactcactataggggaattgtgagcggataacaattcccctctagaaataattttgtttaactttaagaaggagatatacatatggctagcatgactggtggacagcaaatgggtcgcggatccgaattcgagctccgtcgacaagcttgcggccgcactcgagcaccaccaccaccaccactgagatccggctgctaacaaagcccgaaaggaagctgagttggctgctgccaccgctgagcaataactagcataaccccttggggcctctaaacgggtcttgaggggttttttgctgaaaggaggaactatatccggat'

    backbone_start = 'tggcgaatgggacgcgccctgtagcggcgcattaagcgcggcgggtgtggtggttacgcgcagcgtgaccgctacacttgccagcgccctagcgcccgctcctttcgctttcttcccttcctttctcgccacgttcgccggctttccccgtcaagctctaaatcgggggctccctttagggttccgatttagtgctttacggcacctcgaccccaaaaaacttgattagggtgatggttcacgtagtgggccatcgccctgatagacggtttttcgccctttgacgttggagtccacgttctttaatagtggactcttgttccaaactggaacaacactcaaccctatctcggtctattcttttgatttataagggattttgccgatttcggcctattggttaaaaaatgagctgatttaacaaaaatttaacgcgaattttaacaaaatattaacgtttacaatttcaggtggcacttttcggggaaatgtgcgcggaacccctatttgtttatttttctaaatacattcaaatatgtatccgctcatgagacaataaccctgataaatgcttcaataatattgaaaaaggaagagtatgagtattcaacatttccgtgtcgcccttattcccttttttgcggcattttgccttcctgtttttgctcacccagaaacgctggtgaaagtaaaagatgctgaagatcagttgggtgcacgagtgggttacatcgaactggatctcaacagcggtaagatccttgagagttttcgccccgaagaacgttttccaatgatgagcacttttaaagttctgctatgtggcgcggtattatcccgtattgacgccgggcaagagcaactcggtcgccgcatacactattctcagaatgacttggttgagtactcaccagtcacagaaaagcatcttacggatggcatgacagtaagagaattatgcagtgctgccataaccatgagtgataacactgcggccaacttacttctgacaacgatcggaggaccgaaggagctaaccgcttttttgcacaacatgggggatcatgtaactcgccttgatcgttgggaaccggagctgaatgaagccataccaaacgacgagcgtgacaccacgatgcctgcagcaatggcaacaacgttgcgcaaactattaactggcgaactacttactctagcttcccggcaacaattaatagactggatggaggcggataaagttgcaggaccacttctgcgctcggcccttccggctggctggtttattgctgataaatctggagccggtgagcgtgggtctcgcggtatcattgcagcactggggccagatggtaagccctcccgtatcgtagttatctacacgacggggagtcaggcaactatggatgaacgaaatagacagatcgctgagataggtgcctcactgattaagcattggtaactgtcagaccaagtttactcatatatactttagattgatttaaaacttcatttttaatttaaaaggatctaggtgaagatcctttttgataatctcatgaccaaaatcccttaacgtgagttttcgttccactgagcgtcagaccccgtagaaaagatcaaaggatcttcttgagatcctttttttctgcgcgtaatctgctgcttgcaaacaaaaaaaccaccgctaccagcggtggtttgtttgccggatcaagagctaccaactctttttccgaaggtaactggcttcagcagagcgcagataccaaatactgtccttctagtgtagccgtagttaggccaccacttcaagaactctgtagcaccgcctacatacctcgctctgctaatcctgttaccagtggctgctgccagtggcgataagtcgtgtcttaccgggttggactcaagacgatagttaccggataaggcgcagcggtcgggctgaacggggggttcgtgcacacagcccagcttggagcgaacgacctacaccgaactgagatacctacagcgtgagctatgagaaagcgccacgcttcccgaagggagaaaggcggacaggtatccggtaagcggcagggtcggaacaggagagcgcacgagggagcttccagggggaaacgcctggtatctttatagtcctgtcgggtttcgccacctctgacttgagcgtcgatttttgtgatgctcgtcaggggggcggagcctatggaaaaacgccagcaacgcggcctttttacggttcctggccttttgctggccttttgctcacatgttctttcctgcgttatcccctgattctgtggataaccgtattaccgcctttgagtgagctgataccgctcgccgcagccgaacgaccgagcgcagcgagtcagtgagcgaggaagcggaagagcgcctgatgcggtattttctccttacgcatctgtgcggtatttcacaccgcatatatggtgcactctcagtacaatctgctctgatgccgcatagttaagccagtatacactccgctatcgctacgtgactgggtcatggctgcgccccgacacccgccaacacccgctgacgcgccctgacgggcttgtctgctcccggcatccgcttacagacaagctgtgaccgtctccgggagctgcatgtgtcagaggttttcaccgtcatcaccgaaacgcgcgaggcagctgcggtaaagctcatcagcgtggtcgtgaagcgattcacagatgtctgcctgttcatccgcgtccagctcgttgagtttctccagaagcgttaatgtctggcttctgataaagcgggccatgttaagggcggttttttcctgtttggtcactgatgcctccgtgtaagggggatttctgttcatgggggtaatgataccgatgaaacgagagaggatgctcacgatacgggttactgatgatgaacatgcccggttactggaacgttgtgagggtaaacaactggcggtatggatgcggcgggaccagagaaaaatcactcagggtcaatgccagcgcttcgttaatacagatgtaggtgttccacagggtagccagcagcatcctgcgatgcagatccggaacataatggtgcagggcgctgacttccgcgtttccagactttacgaaacacggaaaccgaagaccattcatgttgttgctcaggtcgcagacgttttgcagcagcagtcgcttcacgttcgctcgcgtatcggtgattcattctgctaaccagtaaggcaaccccgccagcctagccgggtcctcaacgacaggagcacgatcatgcgcacccgtggggccgccatgccggcgataatggcctgcttctcgccgaaacgtttggtggcgggaccagtgacgaaggcttgagcgagggcgtgcaagattccgaataccgcaagcgacaggccgatcatcgtcgcgctccagcgaaagcggtcctcgccgaaaatgacccagagcgctgccggcacctgtcctacgagttgcatgataaagaagacagtcataagtgcggcgacgatagtcatgccccgcgcccaccggaaggagctgactgggttgaaggctctcaagggcatcggtcgagatcccggtgcctaatgagtgagctaacttacattaattgcgttgcgctcactgcccgctttccagtcgggaaacctgtcgtgccagctgcattaatgaatcggccaacgcgcggggagaggcggtttgcgtattgggcgccagggtggtttttcttttcaccagtgagacgggcaacagctgattgcccttcaccgcctggccctgagagagttgcagcaagcggtccacgctggtttgccccagcaggcgaaaatcctgtttgatggtggttaacggcgggatataacatgagctgtcttcggtatcgtcgtatcccactaccgagatatccgcaccaacgcgcagcccggactcggtaatggcgcgcattgcgcccagcgccatctgatcgttggcaaccagcatcgcagtgggaacgatgccctcattcagcatttgcatggtttgttgaaaaccggacatggcactccagtcgccttcccgttccgctatcggctgaatttgattgcgagtgagatatttatgccagccagccagacgcagacgcgccgagacagaacttaatgggcccgctaacagcgcgatttgctggtgacccaatgcgaccagatgctccacgcccagtcgcgtaccgtcttcatgggagaaaataatactgttgatgggtgtctggtcagagacatcaagaaataacgccggaacattagtgcaggcagcttccacagcaatggcatcctggtcatccagcggatagttaatgatcagcccactgacgcgttgcgcgagaagattgtgcaccgccgctttacaggcttcgacgccgcttcgttctaccatcgacaccaccacgctggcacccagttgatcggcgcgagatttaatcgccgcgacaatttgcgacggcgcgtgcagggccagactggaggtggcaacgccaatcagcaacgactgtttgcccgccagttgttgtgccacgcggttgggaatgtaattcagctccgccatcgccgcttccactttttcccgcgttttcgcagaaacgtggctggcctggttcaccacgcgggaaacggtctgataagagacaccggcatactctgcgacatcgtataacgttactggtttcacattcaccaccctgaattgactctcttccgggcgctatcatgccataccgcgaaaggttttgcgccattcgatggtgtccgggatctcgacgctctcccttatgcgactcctgcattaggaagcagcccagtagtaggttgaggccgttgagcaccgccgccgcaaggaatggtgcatgcaaggagatggcgcccaacagtcccccggccacggggcctgccaccatacccacgccgaaacaagcgctcatgagcccgaagtggcgagcccgatcttccccatcggtgatgtcggcgatataggcgccagcaaccgcacctgtggcgccggtgatgccggccacgatgcgtccggcgtagaggatcgagatctcgatcccgcgaaattaatacgactcactataggggaattgtgagcggataacaattcccctctagaaataattttgtttaactttaagaaggagatatacat'

    linker = 'gcggccgcactcgagcatcatcaccat'.upper()

    backbone_end = 'caccaccaccaccaccactgagatccggctgctaacaaagcccgaaaggaagctgagttggctgctgccaccgctgagcaataactagcataaccccttggggcctctaaacgggtcttgaggggttttttgctgaaaggaggaactatatccggat'

    vector_name = 'pET21a'
    logger.info('Parent vector is %s' % vector_name)
    logger.info('Complete parent vector sequence with no insert: %s' % complete_backbone)
    logger.info('sequence upstream of insert: %s' % backbone_start)
    logger.info('linker downstream of insert: %s' % linker)
    logger.info('sequence downstream of insert: %s' % backbone_end)

    #assemble constructs
    constructs = []
    for entry in selected_sequences:
        header, prot_seq = entry
        seq = bioseq.Protein(prot_seq.upper()).reverse_translate().reverse_transcribe().upper()
        insert = seq + linker
        validate_insert(insert) #ensure that the insert is ok
        backbone_with_insert = backbone_start + insert + backbone_end
        validate_construct(construct=backbone_with_insert, insert_start_pos=len(backbone_start)+1, target_prot=prot_seq) #ensure that whole construct is ok
        constructs.append( ('%s-%s' % (vector_name, header.split(';')[0]), complete_backbone, backbone_with_insert, insert) )

    #write file
    logger.info('Writing output to "%s".' % constructs_outfile)
    with open(constructs_outfile, 'w') as f:
        f.write('name\tbackbone\tbackbone+insert\tinsert\n')
        for entry in constructs:
            f.write('%s\t%s\t%s\t%s\n' % entry)

    #make file with 5' and 3' UTR for codon fitting by TWIST
    logger.info('Writing output to "%s".' % utr_outfile)
    five_utr = 'ggggaattgtgagcggataacaattcccctctagaaataattttgtttaactttaagaaggagatatacat'
    three_utr = 'gatccggctgctaacaaagcccgaaaggaagctgagttggctgctgccaccgctgagcaataac'
    with open(utr_outfile, 'w') as f:
        f.write("5'UTR\tinsert\t3'UTR\n")
        for entry in selected_sequences:
            header, seq = entry
            seq = bioseq.Protein(seq.upper()).reverse_translate().reverse_transcribe().upper()
            insert = seq+linker
            validate_insert(insert) #ensure that the insert is ok
            f.write('%s\t%s\t%s\n' % (five_utr, insert, three_utr))


    #write file with temperature info for all sequences
    import orgtemp
    with open(info_outfile, 'w') as f:
        f.write("uniprot_id\tec\tdomain\torg\ttemp\tselected\n")
        for entry in all_sequences:
            header, seq = entry
            protid, name, ec, org, source, length = header.split(';')
            organism = ' '.join(org.lower().split()[:2])
            domain = 'not implemented'
            temp = orgtemp.retrieve(organism)
            selected = protid in selected_ids

            f.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (protid, ec, domain, organism, temp, selected))



def write_readme():
    '''Write readme file explaining what the outputs contain.'''
    pass



def make_gb_files():
    '''generate genbank files for these constructs'''
    logger.info('Generating genbank files for sequences.')
    #with open(join(folder, 'sequence_order.tsv'), 'w') as f:
    #    pass




#
# initial_info()
# setup_folders()
# make_fasta()
# filter_fasta()
# build_blast_db()
# blast_fasta_file()
# blast_to_abc()
# build_network()
# cluster()
# cluster_measures()
# make_cluster_fasta_files()
# align_clusters()
# compute_mi()
make_order_file()
make_gb_files()


if __name__ == 'main':
