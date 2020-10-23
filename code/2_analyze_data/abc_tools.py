#!/usr/bin/env python3
"""
abc_tools.py extracts bitscore data from a blast result.

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

import pandas as pd

# read in the abc file
def find_mutual_best_blast_hits(filepath, bitscore_cutoff=0):
	'''
	Go through abc file to find the best blast score between each pair of sequences.
	'''

	data = {}
	all_uid = set([])

	with open(filepath, 'r') as f:
		for line in f:

			# parse the line
			source, target, bitscore = line.strip().split()
			source = source.split(';')[0]
			target = target.split(';')[0]
			bitscore = int(round(float(bitscore)))

			# skip all of the really bad bitscores
			if bitscore < bitscore_cutoff:
				continue

			# make score for same seq 0
			if source == target:
				bitscore = 0

			# keep track of all the uids
			all_uid.add(source)
			all_uid.add(target)

			# add the first reference sequence
			if data.get(source) is None:
				data[source] = {}

			if data.get(target) is None:
				data[target] = {}

			# now add the bitscore
			if data[source].get(target) is None:
				data[source][target] = bitscore
				data[target][source] = bitscore

			elif data[source][target] < bitscore:
				data[source][target] = bitscore
				data[target][source] = bitscore

			else:
				continue
	return data, all_uid


def fill_missing(data, uid_list, missing_value=0):
    '''
    Some sequence pairs can be missing due to the low complexity filter of BLAST.
    Fill these in with values.
    '''
    for source in sorted(uid_list):
        for target in sorted(uid_list):

            if data.get(source) is None:
                data[source] = {}

            if data[source].get(target) is None:
                data[source][target] = missing_value

    return data


def filter_data(data, uid_list):
	'''
	Filter data to retain only desired uniprot identifiers.
	'''
	data_subset = {}
	uid_set = set(uid_list)

	# use uid set to filter the data
	for uid_one in sorted(data.keys()):
		for uid_two in sorted(data[uid_one].keys()):
			if uid_one in uid_set and uid_two in uid_set:

				if data_subset.get(uid_one) is None:
					data_subset[uid_one] = {}

				if data_subset.get(uid_two) is None:
					data_subset[uid_two] = {}

				data_subset[uid_one][uid_two] = data[uid_one][uid_two]
				data_subset[uid_two][uid_one] = data[uid_two][uid_one]

	return data_subset



def to_matrix(data):
    '''
    Convert input data to matrix format.
    '''
    return pd.DataFrame(data)
