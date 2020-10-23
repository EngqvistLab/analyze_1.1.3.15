#!/usr/bin/env python3
"""
Use Regex to extract all the uids that occur in BRENDA pages (indicating that they are characterized)

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

	if not exists(join(FINAL, 'sequence_properties')):
		os.makedirs(join(FINAL, 'sequence_properties'))

#### Your code here ####


import re




def get_all_uniprot_id():
    '''Use regex to get all the uniprot identifiers. Intended as an alternate method that does not depend on parsing the html.'''

	# first check whether the zipped file exists in the scratch folder

	# if it does, compare md5sum and change date

	# only copy the file if these are not the same

	# unzip the file



    #open the list of EC numbers and find all
    filepath = join(RAW_EXTERNAL, folder_name, 'all_enzymes.php.html')
    with open(filepath, 'r') as f:
        data = f.read()
    all_ec = set(re.findall('[0-9]+\.[0-9]+\.[0-9]+\.[0-9a-zA-Z]+', data))

    total = len(list(all_ec))
    print('Number of EC: %s' % total)

    #process each of these
    data = {}
    counter = 0
    for ec in sorted(list(all_ec)):
        print('processing', ec)
        html_doc = join(RAW_EXTERNAL, folder_name , 'html', 'html_data', '%s.html' % ec)

        #read the html page
        with open(html_doc, 'rb') as f:
            document = f.read().decode('utf-8', 'ignore')

        #http://www.uniprot.org/help/accession_numbers
        m = re.findall('[OPQ][0-9](?:[A-Z0-9]){3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}', document)
        data[ec] = set(m)

    #count how many
    with open(join(FINAL, 'sequence_properties', 'ec_uid_characterized_from_html.tsv'), 'w') as f:
        f.write('ec\tuid\n')

        for ec in sorted(data.keys()):
            if data[ec] == set([]):
                f.write('%s\t%s\n' % (ec, 'NA'))
            else:
                for uid in data[ec]:
                    f.write('%s\t%s\n' % (ec, uid))



get_all_uniprot_id()
