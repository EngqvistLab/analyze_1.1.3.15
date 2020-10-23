#!/usr/bin/env python3
"""
Get data from BRENDA

Copyright (C) 2020  Martin Engqvist Lab
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

import os
from os.path import join, dirname, basename, exists, isdir
from dotenv import load_dotenv, find_dotenv

### Load environmental variables from the project root directory ###
# find .env automagically by walking up directories until it's found
dotenv_path = find_dotenv()

# load up the entries as environment variables
load_dotenv(dotenv_path)


# set up directory paths
CURRENT_DIR = os.getcwd()
PROJ = dirname(dotenv_path) # project root directory

DATA = join(PROJ, 'data') #data directory
RAW_EXTERNAL = join(DATA, 'raw_external') # external data raw directory
RAW_INTERNAL = join(DATA, 'raw_internal') # internal data raw directory
INTERMEDIATE = join(DATA, 'intermediate') # intermediate data directory
FINAL = join(DATA, 'final') # final data directory

RESULTS = join(PROJ, 'results') # output directory
FIGURES = join(RESULTS, 'figures') # figure output directory
PICTURES = join(RESULTS, 'pictures') # picture output directory


# make folders specific for certain data
folder_name = 'brenda_data_2019_1'
if folder_name != '':
	#make folders if they don't exist
	if not exists(join(RAW_EXTERNAL, folder_name)):
		os.makedirs(join(RAW_EXTERNAL, folder_name))

	if not exists(join(INTERMEDIATE, folder_name)):
		os.makedirs(join(INTERMEDIATE, folder_name))

	if not exists(join(FINAL, folder_name)):
		os.makedirs(join(FINAL, folder_name))


#### Your code here ####

import multiprocessing


def run_cd_hit(infile, outfile, cutoff, memory):
    '''
    '''

    # get the right word size for the cutoff
    if cutoff < 0.5:
        word = 2
    elif cutoff < 0.6:
        word = 3
    elif cutoff < 0.7:
        word = 4
    else:
        word = 5


    my_cmd = './cdhit-master/cd-hit -i %s -o %s -c %s -n %s -T 1 -M %s' % (infile, outfile, cutoff, word, memory)
    os.system(my_cmd)



def cd_hit_ec(start_filepath, filepath, ec):
    '''
    run cd-hit on an ec number
    '''
    mem = 2000 # memory in mb

    cutoff = 1.0
    outfile = join(filepath, '%s_clustered_sequences_%s.fasta' % (ec, str(int(cutoff*100))))
    if not exists(outfile):
        run_cd_hit(infile=join(start_filepath, '%s.fasta' % ec), outfile=outfile, cutoff=cutoff, memory=mem)

    cutoff = 0.9
    outfile = join(filepath, '%s_clustered_sequences_%s.fasta' % (ec, str(int(cutoff*100))))
    if not exists(outfile):
        run_cd_hit(infile=join(filepath, '%s_clustered_sequences_100.fasta' % ec), outfile=outfile, cutoff=cutoff, memory=mem)

    cutoff = 0.7
    outfile = join(filepath, '%s_clustered_sequences_%s.fasta' % (ec, str(int(cutoff*100))))
    if not exists(outfile):
        run_cd_hit(infile=join(filepath, '%s_clustered_sequences_90.fasta' % ec), outfile=outfile, cutoff=cutoff, memory=mem)

    cutoff = 0.5
    outfile = join(filepath, '%s_clustered_sequences_%s.fasta' % (ec, str(int(cutoff*100))))
    if not exists(outfile):
        run_cd_hit(infile=join(filepath, '%s_clustered_sequences_70.fasta' % ec), outfile=outfile, cutoff=cutoff, memory=mem)

    cutoff = 0.4
    outfile = join(filepath, '%s_clustered_sequences_%s.fasta' % (ec, str(int(cutoff*100))))
    if not exists(outfile):
        run_cd_hit(infile=join(filepath, '%s_clustered_sequences_50.fasta' % ec), outfile=outfile, cutoff=cutoff, memory=mem)



def cd_hit_worker(args):
    '''
    A worker function to enable multiprocessing
    '''
    start_filepath, filepath, ec = args
    cd_hit_ec(start_filepath=start_filepath, filepath=filepath, ec=ec)


def run_all():
    # get all ec numbers
    start_filepath = join(RAW_EXTERNAL, folder_name, 'sequence_data')
    all_ec = [s.replace('.fasta', '') for s in os.listdir(start_filepath)]

    # get all parameters together
    data = []
    for ec in all_ec:
        data.append( (start_filepath, join(INTERMEDIATE, folder_name, 'individual_ec_clustering'), ec) )

    #make a process pool that can be used to distribute work on all cores
    num_cores = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(num_cores)
    pool.map(cd_hit_worker, data)

run_all()
