#!/usr/bin/env python3
"""
A script for combining the KOMODO data with that from DSMZ

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
from os.path import join, dirname, basename, exists, isdir, isfile

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
folder_name = 'komodo_data'
if folder_name != '':
	#make folders if they don't exist
	if not exists(join(RAW_EXTERNAL, folder_name)):
		os.makedirs(join(RAW_EXTERNAL, folder_name))

	if not exists(join(INTERMEDIATE, folder_name)):
		os.makedirs(join(INTERMEDIATE, folder_name))

	if not exists(join(FINAL, folder_name)):
		os.makedirs(join(FINAL, folder_name))


#### Your code here ####


# Download KOMOD page with all the info, save to .tsv file
# http://komodo.modelseed.org/servlet/KomodoTomcatServerSideUtilitiesModelSeed?MediaList


# iterate through the tables using Python library


# extract pH and medium

# write to flatfile

import pandas as pd
import numpy as np


#################################################################################

#setting up variables for flatfiles with databases the code operates on
KOMODO_DB = join(RAW_EXTERNAL, folder_name, 'komodo.txt')
DSMZ_DB = join(RAW_EXTERNAL, folder_name, 'DSM_parsed_taxid_cleaned.txt')
ORGANISM_PH = join(INTERMEDIATE, folder_name, 'organism_pH_averaged.tsv')
#################################################################################
'''
the code aims to merge KOMODO database (media - pH information) with
DSMS database (organism - media information), and extract microorganism - pH pairs.
In case of given pH range, spliting the values and taking a mean value is done
'''

'''
To do/ questions
- clean the code, putting in functions - which steps together?
- in the future: include additional databases and join with KOMODO?
'''


# Check that the data files are there
if not isfile(KOMODO_DB):
	print('The file "%s" does not exist. Please obtain the requisite data before running the script.' % KOMODO_DB)

if not isfile(DSMZ_DB):
	print('The file "%s" does not exist. Please obtain the requisite data before running the script.' % DSMZ_DB)


def _filter_files():
    
  
    '''
    Go through the input files and generate an outut flatfile that will be used in data retrieval.
    '''


    df_komodo = pd.read_table(KOMODO_DB)
    df_DSMZ = pd.read_table(DSMZ_DB, error_bad_lines=False)

    #extracting media number from a link given in a Media column of DSMZ
    df_DSMZ['Medium'] = df_DSMZ['Medium'].map(lambda x: x.lstrip('http://www.dsmz.de/microorganisms/medium/pdf/DSMZ_Medium').rstrip('.pdf'))

    #merging the two databases by Medium number
    df_DSMS_komodo = pd.merge(df_DSMZ, df_komodo, left_on = "Medium", right_on = "ID")

    #creating series with microorganism name and pH, with null values dropped
    df_organism_pH = df_DSMS_komodo[['Name_x', 'PH']].copy()
    df_organism_pH['PH'].replace('', np.nan, inplace=True)
    df_organism_pH.dropna(subset=['PH'], inplace=True)
    #df_organism_pH = df_organism_pH[df_organism_pH.PH != '']
    df_organism_pH = df_organism_pH[df_organism_pH.PH != 'null']
    



    #separation of pH ranges using hyphen into two columns, updating the dataframe with separated columns
    pH_separated = df_organism_pH['PH'].str.split("-").apply(pd.Series)
    df_organism_pH = pd.concat([df_organism_pH.drop('PH', axis = 1), pH_separated], axis = 1)


    #turning separated pH values into numeric- perhaps an easier way to do it?
    #creating a column with average of the two pH value
    #updating dataframe so it only has two columns: name and average pH
    #grouping by microorganism's name and averaging the result, so there is only one pH value for one microorganism
    df_organism_pH[0] = pd.to_numeric(df_organism_pH[0], errors = 'coerce')
    df_organism_pH[1] = pd.to_numeric(df_organism_pH[1], errors = 'coerce')

    df_organism_pH['avg_pH'] = df_organism_pH.mean(axis = 1)

    df_organism_pH = pd.concat([df_organism_pH.Name_x, df_organism_pH.avg_pH], axis = 1)

    by_name = df_organism_pH.groupby(by = 'Name_x')
    df_organism_pH = by_name.mean().round(2)

    df_organism_pH


    #saving a file with microorganism - pH pairs
    df_organism_pH.to_csv(ORGANISM_PH, sep = '\t')






# Check whether the file has been filtered
if not isfile(ORGANISM_PH):
	_filter_files()



def _normalize_name(organism):
    '''
    Normalize a single organism name.
    Should be: "Escherichia coli"
    '''
    if len(organism.split()) < len(organism.split('_')):
        organism = ' '.join(organism.split('_'))
    return organism.lower().capitalize()


def _normalize_org_names(organism_list):
    '''
    Takes a list of organism names and normalizes how they are written.
    Reurns list of normalized organism names
    '''
    return [_normalize_name(x) for x in organism_list]

def get_ph(organism_list):
    '''
    Takes a list of organism names and returns their growth pH values
    '''

    organism_set = set(_normalize_org_names(organism_list))
    out_data = {key:None for key in organism_set}

    with open(ORGANISM_PH, 'r') as file:
        for line in file:
            fields = line.split('\t')
            Name = fields[0].strip()
            pH = fields[1].rstrip()

            if Name in organism_set:
                out_data[Name] = pH

    return out_data
