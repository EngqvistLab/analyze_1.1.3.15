#!/usr/bin/env python3
"""
Modify this line to briefly discribe the functionality of ./brenda_2017_1/2_second_selection.py

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
folder_name = 'BRENDA_second_selection_swissprot'
if folder_name != '':
	#make folders if they don't exist
	if not exists(join(FINAL, folder_name)):
		os.makedirs(join(FINAL, folder_name))


#### Your code here ####

from seqsample import seqsample
import pandas as pd

df = pd.read_csv(join(FINAL, 'experiments', 'uid_experiments.tsv'), sep='\t')



# get a list of all the sequences that worked
worked = list(df.loc[df['active']==True]['node'])

# add a list of SwissProt identifiers that should be "confirmed"
swisprot = ['Q7X2H8',
'P56216',
'Q9UJM8',
'Q9NYQ3',
'E4QP00',
'Q07523',
'Q0PCD7',
'Q9LRR9',
'Q10CE4',
'P0CS93',
'Q6QWR1',
'Q8J2V8',
'Q6UG02',
'Q9NYQ2',
'Q9WU19',
'P12676',
'P05414',
'Q5G234',
'Q53U15',
'P79076',
'O93852',
'P93762',
'O49506',
'P54783',
'P22637',
'Q8J136',
'P04842',
'B5WWZ8',
'Q24JJ8',
'P58710',
'P10867',
'Q8H3I4',
'Q00922',
'P37339',
'Q7FAS1',
'P13006',
'Q9LW56',
'P9WMV9',
'Q7XPR4',
'D7UQ40',
'B8B8K5',
'Q01KC2',
'P59097',
'F2QY27',
'Q6YT73',
'P04841',
'Q9LJH5',
'B8AKX6',
'B8B7C5',
'B8AUI3',
'Q9LRS0',
'Q9AJD6',
'O86963',
'Q8EEB0',
'F2R038',
'Q6NQ66',
'P81156',
'O52792',
'O81030',
'Q9LYD8',
'I1S2N3',
'Q3ZC33',
'Q3ZBW2',
'Q8HXW0',
'Q90YK3',
'C4R702',
'Q54E41',
'P35596',
'P9WMV8',
'Q9KX73',
'Q5B2E9',
'O65709',
'Q9ZWB9',
'Q9ZBU1',
'Q94BP3',
'Q9C614',
'Q75ZP8',
'Q2MF66',
'Q9HDX8',
'O81032',
'Q9FM82',
'Q6FS20',
'Q6CSY3',
'Q7SGY1',
'Q752Y3',
'Q6BZA0',
'Q9CG65',
'Q92452',
'Q5XAK0',
'Q8NZX0',
'Q9FM84',
'Q99YI8',
'B5WWZ9',
'P0DB20',
'P0DB21',
'Q6CG88',
'P21800',
'P16101']

worked.extend(swisprot)


# get a list of all the sequences that did not work
unwanted = list(df.loc[df['active']!=True]['node'])

# # add in other sequences that I don't want to pick (the iron sulfer cluster proteins)
unwanted.extend(['A2CCU9',
'A0A0B0H8K0',
'A0A0J7JD82',
'A5CXB1',
'D4YN11',
'D7FFM6',
'F5SJF4',
'G4CXE1',
'G6FBF9',
'H0RSU3',
'H0SLS5',
'H0TLA1',
'I9Q9Y3',
'I9QC67',
'I9QLS9',
'I9QQ35',
'I9QXN7',
'I9SA22',
'I9SAM5',
'I9SN73',
'I9SW61',
'I9SY86',
'I9TFH2',
'I9THB2',
'I9U1D8',
'I9UMF5',
'I9VIR9',
'I9VJ91',
'I9VVE0',
'I9WXT2',
'I9X7P1',
'I9XDF2',
'I9YBY4',
'I9YPI2',
'I9YRU8',
'I9Z6T9',
'I9ZFV5',
'I9ZW66',
'J0A799',
'J0AF82',
'J0AIG3',
'J0AKV7',
'J0AM02',
'J0B9A4',
'J0BT75',
'J0C7Q0',
'J0CJ25',
'J0FCM7',
'J0HUJ5',
'J0HVI6',
'J0I449',
'J0I969',
'J0IGI0',
'J0INP2',
'J0IR62',
'J0JI36',
'J0JPM0',
'J0K0H0',
'J0K1I9',
'J0KRB8',
'J0LV76',
'J0MA07',
'J0MQ57',
'J0MTC9',
'J0N6S8',
'J0N9A2',
'J0NGH8',
'J0NMR6',
'J0NQ35',
'J0NY54',
'J0PDY3',
'J0Q1Y2',
'J0QF80',
'J0QQU2',
'J0T3Q0',
'J0T775',
'J0UJF0',
'K9NCU5',
'Q1CV18',
'Q3JHC3',
'Q3JQM9',
'U2EQX2',
'S6BUT8',
'Q96Y72'])


# setup logger
seqsample.initial_info()
#
# # make sure I have the filepaths
# seqsample.define_paths()
#
# # compute new selection
# seqsample.mi(data_path=join(FINAL, 'BRENDA'), output_path=join(FINAL, 'BRENDA_second_selection_swissprot'), preselected=worked, unwanted=unwanted)
#
# print('Mi done')

# make a sequence selection and generate order files
seqsample.make_order_file(data_path=join(FINAL, 'BRENDA_second_selection_swissprot'), output_path=join(FINAL, 'BRENDA_second_selection_swissprot'), preselected=worked, unwanted=unwanted, include_preselected=False, mi_cutoff = 90)

print('All done')
