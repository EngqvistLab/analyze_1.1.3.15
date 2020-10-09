#!/usr/bin/env python3
"""
The script calculates all pairwise identities (based on MUSCLE alignment) and returns as a pandas data frame.

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
folder_name = 'brenda_2017_1'
if folder_name != '':
	#make folders if they don't exist
	if not exists(join(RAW_EXTERNAL, folder_name)):
		os.makedirs(join(RAW_EXTERNAL, folder_name))

	if not exists(join(INTERMEDIATE, folder_name)):
		os.makedirs(join(INTERMEDIATE, folder_name))

	if not exists(join(FINAL, folder_name)):
		os.makedirs(join(FINAL, folder_name))


#### Your code here ####

from dnapy.resources import fasta
import itertools
import pandas as pd

import io
from Bio.Align.Applications import MuscleCommandline

# Get a list of all selected
def list_ordered_sequences():
	'''
	Get a list of ordered sequences.
	'''
	with open(join(FINAL, 'brenda_2017_1', 'selected_ids.tsv'), 'r') as f:
		all_uid = f.read().split('\n')

	return all_uid


# Get the sequences for the selected uids
def get_ec_to_uid_to_seq(uid_list):
	'''
	Return nested dicitonary with EC number primary keys and uid values.
	In the second level each uid holds its protein sequence.
	Only return for uniprot identifiers in the list.
	'''
	data = {}

	for record in fasta.parse_file(join(INTERMEDIATE, 'brenda_2017_1', '1_1_3__BRENDA_sequences_filtered_brenda_2017_1.fasta')):
		header, sequence = record
		header_data = header.lstrip('>').split(';')
		uid = header_data[0]
		ec = header_data[2]

		if uid in uid_list:
			if data.get(ec) is None:
				data[ec] = {}

			data[ec][uid] = sequence

	return data


# Do all pairwise alignments
def create_combinations(uid_list):
	'''
	Generate all pairwise combinations of uids.
	'''
	return itertools.combinations_with_replacement(uid_list, 2)


def pairwise_align(id1, id2, seq1, seq2):
	'''
	Perform pairwise alignment.
	'''
	records = '>%s\n%s\n>%s\n%s' % (id1, seq1, id2, seq2) #prepare 'virtual' FASTA file

	records_handle = io.StringIO(records) #turn string into a handle
	tempdata = records_handle.getvalue()
	muscle_cline = MuscleCommandline()
	stdout, stderr = muscle_cline(stdin=tempdata)
	aln = fasta.parse_string(stdout)
	output = []
	for entry in aln:
		output.append(entry[0])
		output.append(entry[1])
	return output


# Get percent identity of the aligned sequences
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





def calculate_identities_for_all(seq_data):
	'''
	'''
	# assemble the complete data set
	all_seq_data = {}
	for ec in seq_data.keys():
		all_seq_data.update(seq_data[ec])


	identity_data = {}

	# make all pairs of combinations of the uids
	all_uid = list(all_seq_data.keys())
	uid_pairs = create_combinations(all_uid)

	for pair in uid_pairs:
		id1 = '>%s' % pair[0]
		id2 = '>%s' % pair[1]
		seq1 = all_seq_data[pair[0]]
		seq2 = all_seq_data[pair[1]]
		output = pairwise_align(id1=id1, id2=id2, seq1=seq1, seq2=seq2)

		aln_seq1 = output[1]
		aln_seq2 = output[3]

		identity = pair_ident(Seq1=aln_seq1, Seq2=aln_seq2, single_gaps=True)

		if identity_data.get(pair[0]) is None:
			identity_data[pair[0]] = {}

		if identity_data.get(pair[1]) is None:
			identity_data[pair[1]] = {}

		identity_data[pair[0]][pair[1]] = identity
		identity_data[pair[1]][pair[0]] = identity

	return identity_data



def save_identity_matrix_all(data):
	'''
	Save the identity data as a flatfile matrix.
	'''
	# Convert identity values to data frame and save to flatfile
	df = pd.DataFrame(data)
	df = df.fillna(0)
	df.to_csv(path_or_buf=join(FINAL, 'brenda_2017_1', 'all_identity_matrix.tsv'), sep='\t')


def save_network_all(data):
	'''
	Write a network file for cytoscape.
	'''
	used_pairs = set([])
	out_data = ['source\ttarget\tinteraction\tdirected\tvalue']

	for uid1 in data.keys():
		for uid2 in data[uid1].keys():
			if uid1 == uid2:
				continue

			pair = str(sorted([uid1, uid2]))
			if pair in used_pairs:
				continue

			used_pairs.add(pair)

			out_data.append('%s\t%s\t%s\t%s\t%s' % (uid1, uid2, 'identity', 'false', data[uid1][uid2]))

	with open(join(FINAL, 'visualization_data', 'all_identity_network.tsv'), 'w') as f:
		f.write('\n'.join(out_data))


def calculate_identities_within_ec(seq_data):
	'''
	'''
	identity_data = {}

	#for each ec
	for ec in seq_data.keys():
		identity_data[ec] = {}

		# make all pairs of combinations of the uids
		uid_in_ec = list(seq_data[ec].keys())
		uid_pairs = create_combinations(uid_in_ec)

		for pair in uid_pairs:
			id1 = '>%s' % pair[0]
			id2 = '>%s' % pair[1]
			seq1 = seq_data[ec][pair[0]]
			seq2 = seq_data[ec][pair[1]]
			output = pairwise_align(id1=id1, id2=id2, seq1=seq1, seq2=seq2)

			aln_seq1 = output[1]
			aln_seq2 = output[3]

			identity = pair_ident(Seq1=aln_seq1, Seq2=aln_seq2, single_gaps=True)

			if identity_data[ec].get(pair[0]) is None:
				identity_data[ec][pair[0]] = {}

			if identity_data[ec].get(pair[1]) is None:
				identity_data[ec][pair[1]] = {}


			identity_data[ec][pair[0]][pair[1]] = identity
			identity_data[ec][pair[1]][pair[0]] = identity

	return identity_data


def save_identity_matrix_ec(data):
	'''
	Save the identity data as a flatfile matrix.
	'''
	# Convert identity values to data frame and save to flatfile
	for ec in data.keys():
		df = pd.DataFrame(data[ec])
		df = df.fillna(0)
		df.to_csv(path_or_buf=join(FINAL, 'visualization_data', '%s_identity_matrix.tsv' % ec), sep='\t')


def save_network_ec(data):
	'''
	Write a network file for cytoscape.
	'''
	for ec in data.keys():
		used_pairs = set([])
		out_data = ['source\ttarget\tinteraction\tdirected\tvalue']

		for uid1 in data[ec].keys():
			for uid2 in data[ec][uid1].keys():
				if uid1 == uid2:
					continue

				pair = str(sorted([uid1, uid2]))
				if pair in used_pairs:
					continue

				used_pairs.add(pair)

				out_data.append('%s\t%s\t%s\t%s\t%s' % (uid1, uid2, 'identity', 'false', data[ec][uid1][uid2]))

		with open(join(FINAL, 'visualization_data', '%s_identity_network.tsv' % ec), 'w') as f:
			f.write('\n'.join(out_data))


# get all ordered uids
ordered_uid = list_ordered_sequences()

# get ec mapping and the sequences
ec_uid_seq = get_ec_to_uid_to_seq(ordered_uid)




# calculate the pairwise identities
all_ident_data = calculate_identities_for_all(ec_uid_seq)

# save the whole identity matrix
save_identity_matrix_all(all_ident_data)

# save the data in network format
save_network_all(all_ident_data)



# calculate the pairwise identities
ec_ident_data = calculate_identities_within_ec(ec_uid_seq)

# save the whole identity matrix
save_identity_matrix_ec(ec_ident_data)

# save the data in network format
save_network_ec(ec_ident_data)
