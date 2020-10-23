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
import fasta
import itertools
import io
from Bio.Align.Applications import MuscleCommandline
from os.path import join
import pandas as pd


# Get a list of all selected
def get_sequences():
	'''
	Get all sequencs as a dictionary
	'''
	fasta_data = fasta.parse_file('1_1_3__BRENDA_sequences_filtered_2017_1.fasta')

	seq_data = {}
	for record in fasta_data:
		header, seq = record
		seq_data[header] = seq

	return seq_data


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



def save_identity_matrix_all(data):
	'''
	Save the identity data as a flatfile matrix.
	'''
	# Convert identity values to data frame and save to flatfile
	df = pd.DataFrame(data)
	df = df.fillna(0)
	df.to_csv(path_or_buf='all_identity_matrix_ident.tsv', sep='\t')


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

	with open('all_identity_network.tsv', 'w') as f:
		f.write('\n'.join(out_data))





# get all sequences
seq_data = get_sequences()

# make all sequence combinations
combos = create_combinations(seq_data.keys())

# generate all pairs and align them
counter = 0
identity_data = {}
for pair in combos:
	counter += 1
	if counter % 10000 == 0:
		print('%s done' % counter)
	id1, id2 = pair
	seq1 = seq_data[id1]
	seq2 = seq_data[id2]

	# keep only the uniprot identifiers
	id1 = id1.split(';')[0]
	id2 = id2.split(';')[0]


	# perform the alignment
	output = pairwise_align(id1, id2, seq1, seq2)

	aln_seq1 = output[1]
	aln_seq2 = output[3]

	identity = pair_ident(Seq1=aln_seq1, Seq2=aln_seq2, single_gaps=True)

	if identity_data.get(id1) is None:
		identity_data[id1] = {}

	if identity_data.get(id2) is None:
		identity_data[id2] = {}

	identity_data[id1][id2] = identity
	identity_data[id2][id1] = identity



save_identity_matrix_all(identity_data)
save_network_all(identity_data)
