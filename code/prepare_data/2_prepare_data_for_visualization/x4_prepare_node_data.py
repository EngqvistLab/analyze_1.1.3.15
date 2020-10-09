#!/usr/bin/env python3
"""
Collect various node data and make cytoscape-compatible data files that can be imported.

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
folder_name = 'visualization_data'
if folder_name != '':
	#make folders if they don't exist
	if not exists(join(RAW_EXTERNAL, folder_name)):
		os.makedirs(join(RAW_EXTERNAL, folder_name))

	if not exists(join(INTERMEDIATE, folder_name)):
		os.makedirs(join(INTERMEDIATE, folder_name))

	if not exists(join(FINAL, folder_name)):
		os.makedirs(join(FINAL, folder_name))


#### Your code here ####


# get ec for each
def map_ec_to_uid():
	'''
	Make an attribute file mapping uids to ec numbers
	'''
	with open(join(FINAL, folder_name, 'uid_ec.attr'), 'w') as fo:
		fo.write('node\tec\n')

		with open(join(INTERMEDIATE, folder_name, 'ec_uid_org_from_fasta.tsv'), 'r') as fi:
			fi.readline()

			for line in fi:
				ec, uid, org = line.strip().split('\t')
				fo.write('%s\t%s\n' % (uid, ec))



# get whether they had been characterized
def map_brenda_to_uid():
	'''
	Make an attribute file mapping uids to whether they have been characterized.
	Based on BRENDA html files.
	'''
	with open(join(FINAL, folder_name, 'uid_characterized.attr'), 'w') as fo:
		fo.write('node\tcharacterized_in_brenda\n')

		with open(join(INTERMEDIATE, folder_name, 'ec_uid_characterized_from_html.tsv'), 'r') as fi:
			fi.readline()

			for line in fi:
				ec, uid = line.strip().split('\t')
				fo.write('%s\t%s\n' % (uid, 'True'))


# get which cluster they were in
def mcl_clusters():
	'''
	Get data regarding which cluster the sequences were in.
	'''
	pass


# get which were selected
def ordered_uids():
	'''
	Make an attribute file mapping uids to whether they were selected for synthesis.
	'''
	with open(join(FINAL, folder_name, 'uid_ordered.attr'), 'w') as fo:
		fo.write('node\tordered\n')

		with open(join(FINAL, 'BRENDA', 'selected_ids.tsv'), 'r') as fi:

			for line in fi:
				uid = line.strip()
				fo.write('%s\t%s\n' % (uid, 'True'))


# get which were synthesized, expressed, and active
def experimental_data():
	'''
	Get info from own experiments.
	'''
	pass



map_ec_to_uid()
map_brenda_to_uid()
mcl_clusters()
ordered_uids()
experimental_data()
