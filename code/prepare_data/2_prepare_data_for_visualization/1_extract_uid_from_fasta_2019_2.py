#!/usr/bin/env python3
"""
Go through the BRENDA fasta file and extract all the uniprot identifiers, their EC number, as well as the organism names.

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

import os
from dotenv import load_dotenv, find_dotenv # do 'pip install python-dotenv'
from os.path import join, dirname, basename, exists, isdir

### Load environmental variables from the project root directory ###
# find .env automagically by walking up directories until it's found
dotenv_path = find_dotenv()

# load up the entries as environment variables
load_dotenv(dotenv_path)

# now you can get the variables using their names

# Check whether a network drive has been specified
DATABASE = os.environ.get("NETWORK_URL")
if DATABASE == 'None':
	pass
else:
	pass
	#mount network drive here

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
folder_name = 'brenda_2019_2'
if folder_name != '':
	#make folders if they don't exist
	if not exists(join(RAW_EXTERNAL, folder_name)):
		os.makedirs(join(RAW_EXTERNAL, folder_name))

	if not exists(join(INTERMEDIATE, folder_name)):
		os.makedirs(join(INTERMEDIATE, folder_name))

	if not exists(join(FINAL, folder_name)):
		os.makedirs(join(FINAL, folder_name))


#### Your code here ####
from orgtools import topfunctions
from Bio import SeqIO

def make_fasta():
	'''
	Go through and collect sequnces from 1.1.3.n and
	assemble into a new file.
	'''
	all_data = []
	for fi in sorted(os.listdir(join(RAW_EXTERNAL, folder_name, 'sequence_data'))):
		if fi.startswith('1.1.3.'):
			with open(join(RAW_EXTERNAL, folder_name, 'sequence_data', fi), 'r') as f:
				for line in f:
					all_data.append(line)

	with open(join(INTERMEDIATE, 'brenda_2019_2', '1_1_3__BRENDA_sequences_2019_2.fasta'), 'w') as f:
		f.write('\n'.join(all_data))



def filter_fasta(min_len = 200, max_len = 580):
    '''Remove sequences that are too long or too short, also sequenes that have X in them'''
    infile = join(INTERMEDIATE, 'brenda_2019_2', '1_1_3__BRENDA_sequences_2019_2.fasta')
    outfile = join(INTERMEDIATE, 'brenda_2019_2', '1_1_3__BRENDA_sequences_filtered_2019_2.fasta')

    retained = 0
    removed_len = 0
    removed_x = 0
    removed_m = 0
    with open(outfile, 'w') as f:
        for record in SeqIO.parse(infile, 'fasta'):
            header = record.description
            seq = str(record.seq).lower()

            if min_len <= len(seq) <= max_len and 'x' not in seq and seq[0] == 'm':
                retained += 1
                f.write('>%s\n%s\n' % (header, seq))
            elif 'x' in seq:
                removed_x += 1
            elif seq[0] != 'm':
                removed_m += 1
            else:
                removed_len += 1



def get_uid_from_fasta():
	'''
	Go through fasta file from BRENDA and extract the uid, organism and ec number.
	'''
	data = {}

	# open file and go through it
	with open(join(INTERMEDIATE, 'brenda_2019_2', '1_1_3__BRENDA_sequences_filtered_2019_2.fasta'), 'r') as f:
		for line in f:

			# only look at header lines
			if line.startswith('>'):
				line_data = line.lstrip('>').rstrip().split(';')
				uid = line_data[0]
				ec = line_data[2]
				org = ' '.join(line_data[3].split()[:2])

				# add ec key if not present
				if data.get(ec) is None:
					data[ec] = {}

				if data[ec].get(org) is None:
					data[ec][org] = []

				# add uid to data structure
				data[ec][org].append(uid)

	return data



def write_outfile(data):
	'''
	Write the parsed data to an outfile.
	'''
	with open(join(FINAL, folder_name, 'ec_uid_org_from_fasta_2019_2.tsv'), 'w') as f:
		f.write('uid\tec\n')

		for ec in sorted(data.keys()):
			for org in sorted(data[ec].keys()):
				for entry in data[ec][org]:
					f.write('%s\t%s\n' % (entry, ec))



def write_info_outfile(data):
	'''
	Write the parsed data to an outfile.
	'''
	# first collect all organism names
	all_ids = []
	for ec in sorted(data.keys()):
		for org in sorted(data[ec].keys()):
			all_ids.extend(data[ec][org])


	# get taxid and lineage data
	properties_object = topfunctions.Properties(all_ids)
	properties_object.flatfile(join(FINAL, folder_name, '1-1-3-n_identifier_info_2019_2.tsv'))


make_fasta()

filter_fasta()


data = get_uid_from_fasta()

write_outfile(data)
write_info_outfile(data)
