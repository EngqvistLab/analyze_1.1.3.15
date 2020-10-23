#!/usr/bin/env python3
"""
1_get_properties.py leverages a range of other scripts to assemble information for uniprot identifiers.

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
folder_name = 'sequence_properties'
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


def list_ordered_sequences():
    '''
    Get a list of the sequences I ordered.
    '''
    print('Getting a list of selected UniProt ids ...')
    data = []
    with open(join(FINAL, 'brenda_2017_1', 'selected_ids.tsv'), 'r') as f:
        for line in f:
            data.append(line.strip())
    print('Done')

    return data


# get the identifiers
identifier_data = list_ordered_sequences()

# get all the associated data
properties_object = topfunctions.Properties(identifier_data)

# save data to flatfile
properties_object.flatfile(join(FINAL, 'brenda_2017_1', 'uid_data.tsv'))
