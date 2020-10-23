#!/usr/bin/env python3
"""
Modify this line to briefly discribe the functionality of prepare_visualization/prepare_graph.py

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





def get_all_ec():
	'''
	Get a list of all ec numbers in data file.
	'''
	ec_set = set([])
	with open(join(INTERMEDIATE, folder_name, 'ec_uid_org_from_fasta.tsv'), 'r') as f:
		f.readline()

		for line in f:
			ec, uid, org = line.split('\t')
			ec_set.add(ec)
	return ec_set



def write_network_ouput_file(data):
	'''
	Write an network file that can be used for cytoscape visualization.
	'''

	with open(join(FINAL, folder_name, 'network.tsv'), 'w') as f:
		out_list = ['source\ttype\ttarget\tvalue']
		for uid_one in sorted(data.keys()):

			for uid_two in sorted(data[uid_one].keys()):

				out_list.append('%s\tbitscore\t%s\t%s' % (uid_one, uid_two, data[uid_one][uid_two]))

		#write it
		f.write('\n'.join(out_list))


#
# def write_json_outfile(data, uids, ID):
# 	'''
# 	For visualization with D3js
# 	Write the json data file
# 	'''
#
# 	uids = sorted(uids)
#
# 	filename = '%s.json' % ID
# 	with open(join(FINAL, folder_name, filename), 'w') as f:
# 		group = '"#006d2c"'
# 		#group=uids[uid]
#
# 		#write the nodes
# 		f.write('{\n  "nodes":[\n')
#
# 		out_list = []
# 		for uid in uids:
# 			out_list.append('    {"name":"%s","group":%s}' % (uid, group))
#
# 		#write it
# 		f.write(',\n'.join(out_list))
# 		f.write('\n  ],\n')
#
#
#
# 		#write the links
# 		f.write('  "links":[\n')
#
# 		out_list = []
# 		for uid_one in sorted(data.keys()):
# 			uid_one_pos = uids.index(uid_one)
#
# 			for uid_two in sorted(data[uid_one].keys()):
# 				uid_two_pos = uids.index(uid_two)
#
# 				out_list.append('    {"source":%s,"target":%s,"value":%s}' % (uid_one_pos, uid_two_pos, data[uid_one][uid_two]))
#
# 		#write it
# 		f.write(',\n'.join(out_list))
# 		f.write('\n  ]\n}')
#
#
#
# def write_html_file(ID):
# 	'''
# 	For visualization with D3js
# 	Write the html file that goes with the json data.
# 	'''
# 	filename = '%s.html' % ID
# 	with open(join(FINAL, folder_name, filename), 'w') as f:
# 		html_string = '''
# <!DOCTYPE html>
# <meta charset="utf-8">
# <style>
#
# .node {
#   stroke: #fff;
#   stroke-width: 1.5px;
# }
#
# .link {
#   stroke: #999;
#   stroke-opacity: .6;
# }
#
# </style>
# <body>
# <script src="http://d3js.org/d3.v3.min.js"></script>
# <script>
#
# var width = 900,
#     height = 900;
#
# var color = d3.scale.category20();
#
# var force = d3.layout.force()
#     .charge(-90)
#     .linkDistance(50)
#     .size([width, height]);
#
# var svg = d3.select("body").append("svg")
#     .attr("width", width)
#     .attr("height", height);
#
# d3.json("%s.json", function(error, graph) {
#   force
#       .nodes(graph.nodes)
#       .links(graph.links)
#       .start();
#
#   var link = svg.selectAll(".link")
#       .data(graph.links)
#     .enter().append("line")
#       .attr("class", "link")
#       .style("stroke-width", function(d) { return Math.sqrt(d.value); });
#
#   var node = svg.selectAll(".node")
#       .data(graph.nodes)
#     .enter().append("circle")
#       .attr("class", "node")
#       .attr("r", 5)
#       .style("fill", function(d) { return color(d.group); })
#       .call(force.drag);
#
#   node.append("title")
#       .text(function(d) { return d.name; });
#
#   force.on("tick", function() {
#     link.attr("x1", function(d) { return d.source.x; })
#         .attr("y1", function(d) { return d.source.y; })
#         .attr("x2", function(d) { return d.target.x; })
#         .attr("y2", function(d) { return d.target.y; });
#
#     node.attr("cx", function(d) { return d.x; })
#         .attr("cy", function(d) { return d.y; });
#   });
# });
#
# </script>
# ''' % ID
# 		f.write(html_string)
#


# get all the bitscore data
bitscore_data, bitscore_uids = abc_tools.find_mutual_best_blast_hits(filepath=join(INTERMEDIATE, 'BRENDA', 'all-vs-all_bitscore.abc'), bitscore_cutoff=20)

# write output files for cytoscape
write_network_ouput_file(bitscore_data)
