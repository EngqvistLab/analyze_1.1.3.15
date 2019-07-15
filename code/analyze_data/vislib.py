#!/usr/bin/env python3
"""
This library draws visualizations from enzyme data.

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


#### Your code here ####

import wsvg
from colcol import scales
import math
import os




class _DrawingBaseClass(object):
	'''
	'''
	def __init__(self, data, filepath, main, group_dict, group_colors, text_color):
		'''
		'''
		assert type(data) is dict, 'Error, the variable "data" must be a dictionary with identifier keys and sequence values.'
		assert type(filepath) is str, 'Error, the variable "filepath" must be a string.'

		self.data = data
		self.filepath = filepath
		self.main = main
		self.group_dict = group_dict

		self.id_list = sorted(list(self.data.keys()))
		self.num_sequences = len(self.id_list)

		self.fontsize = 7 #font size
		self.x_spacing = 80 # space to the left
		self.y_spacing = 40 # space on the top
		self.text_color = text_color

		if group_dict is not None:
			assert type(group_dict) is dict, 'Error, the variable "group_dict" must be a dictionary with identifier keys and string values to identify group belonging.'
			self._setup_group_variable()
			self.usegroups = True

			# assign group colors
			if group_colors is None:
				self.group_colors = scales.rainbow(self.number_of_groups)
			else:
				self.group_colors = group_colors

			assert self.number_of_groups == len(self.group_colors), 'Error, the number of colors supplied in "group_colors" must match the number of groups.'

		else:
			self.number_of_groups = 0
			self.usegroups = False


	def _setup_group_variable(self):
		'''
		'''
		self.all_groups = sorted(list(set(self.group_dict.values())))
		self.number_of_groups = len(self.all_groups)

		uid_order_list = []
		for group in self.all_groups:
			for uid in self.id_list:
				if self.group_dict[uid] == group:
					uid_order_list.append(uid)

		self.id_list = uid_order_list


	def _start_plotting(self):
		'''
		Start an SVG scene to which plotting elements can be added.
		'''
		self.scene = wsvg.Scene(name=self.filepath, size=(self.page_x_size, self.page_y_size))


	def _title(self):
		'''
		'''
		self.scene.add(wsvg.Text(text=self.main, origin=(self.page_x_size/2, self.y_spacing/1.5), angle=0, size=self.fontsize*1.5, weight="bold", color=self.text_color, anchor="center", align='center'))


	def _stop_plotting(self):
		'''
		'''
		self.scene.write_svg()


	def _convert_to_pdf(self):
		'''
		'''
		# for one file
		# inkscape t.svg --export-pdf=t.pdf

		# for many files
		#$ (echo t.svg --export-pdf=t.pdf; echo u.svg --export-pdf=u.pdf) | DISPLAY= inkscape --shell

		mycmd = 'inkscape %s --export-pdf=%s' % (self.filepath, self.filepath.replace('.svg', '.pdf'))
		os.system(mycmd)



class network(_DrawingBaseClass):
	'''
	'''
	def __init__(self, data, filepath, main=None, group_dict=None, group_colors=None, text_color='#424242'):
		'''
		Takes ...... input data .....
		filepath specifies the output file
		main is optional and is a string that specifies the plot title
		group_dict is optional, it is a dictionary with identifier keys and string values that indicate group belonging
		'''
		_DrawingBaseClass.__init__(self, data=data, filepath=filepath, main=main, group_dict=group_dict, group_colors=group_colors, text_color=text_color)

		# setup page size
		self.page_x_size = 744.09
		self.page_y_size = 1052.36

		# setup the svg scene
		self._start_plotting()

		# finish document
		self._stop_plotting()

		# convert to pdf
		self._convert_to_pdf()



class operons(_DrawingBaseClass):
	'''
	'''
	def __init__(self, data, filepath, main=None, group_dict=None, group_colors=None, text_color='#424242'):
		'''
		Takes ...... input data .....
		filepath specifies the output file
		main is optional and is a string that specifies the plot title
		group_dict is optional, it is a dictionary with identifier keys and string values that indicate group belonging
		'''
		_DrawingBaseClass.__init__(self, data=data, filepath=filepath, main=main, group_dict=group_dict, group_colors=group_colors, text_color=text_color)

		# setup page size
		self.page_x_size = 744.09
		self.page_y_size = 1052.36

		# setup the svg scene
		self._start_plotting()

		# finish document
		self._stop_plotting()

		# convert to pdf
		self._convert_to_pdf()



class domains(_DrawingBaseClass):
	'''
	'''
	def __init__(self, data, filepath, main=None, group_dict=None, group_colors=None, text_color='#424242'):
		'''
		Takes ...... input data .....
		filepath specifies the output file
		main is optional and is a string that specifies the plot title
		group_dict is optional, it is a dictionary with identifier keys and string values that indicate group belonging
		'''
		_DrawingBaseClass.__init__(self, data=data, filepath=filepath, main=main, group_dict=group_dict, group_colors=group_colors, text_color=text_color)

		# setup page size
		self.page_x_size = 744.09
		self.page_y_size = 1052.36

		# setup the svg scene
		self._start_plotting()

		# finish document
		self._stop_plotting()

		# convert to pdf
		self._convert_to_pdf()



class substrate(_DrawingBaseClass):
	'''
	'''
	def __init__(self, data, dend_data, filepath, subst_data=None, subst_data_replicate=None, main=None, group_dict=None, group_colors=None, property_dict=None, property_colors=None, text_color='#424242'):
		'''
		Takes ...... input data .....
		filepath specifies the output file
		main is optional and is a string that specifies the plot title
		group_dict is optional, it is a dictionary with identifier keys and string values that indicate group belonging
		'''
		_DrawingBaseClass.__init__(self, data=data, filepath=filepath, main=main, group_dict=group_dict, group_colors=group_colors, text_color=text_color)

		self.dend_data = dend_data
		self.subst_data = subst_data
		self.subst_data_replicate = subst_data_replicate
		self.group_dict = group_dict
		self.property_dict = property_dict
		self.property_colors = property_colors

		# calculate dendrogram size
		self._calculate_dendrogram_size()

		# calculate the rectangle size for the substrate plot
		self._calculate_column_width()

		# setup page size
		self.page_x_size = 1052.36
		self.page_y_size = 744.09

		# space between dendrogram and actiity data
		self.space_between_plots = 90

		# dendrogram offset
		self.x_offset = self.page_x_size / 2 - self.dend_x_size / 2
		self.y_offset = self.page_y_size - self.dend_y_size - self.y_spacing

		# adjust dendrogram data coordinates
		self._adjust_data_offset()

		# the data is displayed in rows, calculate indexes for where these should go
		self._calculate_row_index_positions()

		# setup the svg scene
		self._start_plotting()

		# draw any category information (sequence info)
		self.indexes_used = []
		self._draw_seq_info()

		# draw the substrate activity component
		self.all_subst = []
		if self.subst_data is not None:
			self._draw_substrate_activity_info()

		# draw the vertical grey support lines
		self._draw_vertical_guides()

		# draw the howizontal grey support lines
		#self._draw_horizontal_guides()

		# draw the dendrogram component
		self._draw_dendrogram()

		# set title if present
		if main is not None:
			assert type(main) is str, 'Error, the variable "main" must be a string.'
			self._title()

		# finish document
		self._stop_plotting()

		# convert to pdf
		self._convert_to_pdf()


	def _calculate_dendrogram_size(self):
		'''
		Calculate how big the dendrogam is on the x- and y-axis
		'''
		all_x = []
		all_y = []
		for polygon in self.dend_data:
			point_data = []
			for point in polygon:
				x, y = point
				all_x.append(x)
				all_y.append(y)

		self.dend_min_x = min(all_x)
		self.dend_max_x = max(all_x)
		self.dend_min_y = min(all_y)
		self.dend_max_y = max(all_y)
		self.dend_x_size = self.dend_max_x - self.dend_min_x
		self.dend_y_size = self.dend_max_y - self.dend_min_y


	def _calculate_column_width(self):
		'''
		Calculate column width used for plotting
		'''
		self.rect_size = self.dend_x_size / len(self.id_list)


	def _calculate_row_index_positions(self):
		'''
		Calculate where data rows should go.
		'''
		self.index_y_coordinates = {}
		max_indexes = math.floor((self.page_y_size - self.y_spacing*2 - self.space_between_plots - self.dend_y_size)/self.rect_size)

		y = self.dend_min_y
		for i in range(0, max_indexes):
			self.index_y_coordinates[i] = y - self.space_between_plots - i*self.rect_size


	def _adjust_data_offset(self):
		'''
		Modify dendrogram coordinates by offset
		'''
		# adjust dendrogram coordinates with the offset
		self.dend_data_adjusted = []
		for polygon in self.dend_data:
			point_data = []
			for point in polygon:
				x, y = point
				point_data.append((x+self.x_offset, y+self.y_offset))
			self.dend_data_adjusted.append(point_data)

		# adjust label coordinates with the offset
		self.data_adjusted = {}
		for key in sorted(self.data.keys()):
			x, y = self.data[key]
			self.data_adjusted[key] = (x+self.x_offset, y+self.y_offset)

		self.dend_min_x += self.x_offset
		self.dend_max_x += self.x_offset
		self.dend_min_y += self.y_offset
		self.dend_max_y += self.y_offset


	def _draw_vertical_guides(self):
		'''
		Draw the grey vertical guide lines.
		'''
		guide_col = '#DDDDDD'

		if self.subst_data is None:
			y_start = self.page_y_size - self.y_spacing - self.dend_y_size - self.space_between_plots - self.rect_size*len(self.all_subst) - (len(self.indexes_used)-1)*self.rect_size
		else:
			y_start = self.page_y_size - self.y_spacing - self.dend_y_size - self.space_between_plots - self.rect_size*len(self.all_subst) - (len(self.indexes_used)+1)*self.rect_size

		y_end = self.page_y_size - self.y_spacing - self.dend_y_size

		for uid in self.data_adjusted.keys():
			label_x, label_y = (self.data_adjusted[uid])
			start = (label_x-self.rect_size/2, y_start)
			end = (label_x-self.rect_size/2, y_end)
			self.scene.add(wsvg.Line(start=start, end=end, line_color=guide_col, line_width=0.75))

		self.scene.add(wsvg.Line(start=(self.dend_max_x+self.rect_size/2, y_start), end=(self.dend_max_x+self.rect_size/2, y_end), line_color=guide_col, line_width=0.75))


	def _draw_horizontal_guides(self):
		'''
		Draw grey horizontal guide lines
		'''
		guide_col = '#DDDDDD'

		if self.property_dict is not None:
			self.scene.add(wsvg.Line(start=(self.dend_min_x-self.rect_size/2, self.index_y_coordinates[0]+self.rect_size), end=(self.dend_max_x+self.rect_size/2, self.index_y_coordinates[0]+self.rect_size), line_color=guide_col, line_width=0.75))
			self.scene.add(wsvg.Line(start=(self.dend_min_x-self.rect_size/2, self.index_y_coordinates[len(self.indexes_used)]+self.rect_size), end=(self.dend_max_x+self.rect_size/2, self.index_y_coordinates[len(self.indexes_used)]+self.rect_size), line_color=guide_col, line_width=0.75))

		self.scene.add(wsvg.Line(start=(self.dend_min_x-self.rect_size/2, self.index_y_coordinates[min(self.subst_index.values())]+self.rect_size), end=(self.dend_max_x+self.rect_size/2, self.index_y_coordinates[min(self.subst_index.values())]+self.rect_size), line_color=guide_col, line_width=0.75))
		self.scene.add(wsvg.Line(start=(self.dend_min_x-self.rect_size/2, self.index_y_coordinates[max(self.subst_index.values())]), end=(self.dend_max_x+self.rect_size/2, self.index_y_coordinates[max(self.subst_index.values())]), line_color=guide_col, line_width=0.75))


	def _draw_dendrogram(self):
		'''
		Draw dendrogram with the labels
		'''
		for polygon in self.dend_data_adjusted:
			self.scene.add(wsvg.PolyLine(points=polygon, fill_color=None, line_color='#111111', line_width=1))

		for key in sorted(self.data_adjusted.keys()):
			x, y = self.data_adjusted[key]
			self.scene.add(wsvg.Text(text=key, origin=(x+self.fontsize*0.30, y-self.fontsize), angle=-90, size=self.fontsize, weight="normal", color=self.text_color, anchor="start", align='start'))


	def _draw_datapoint(self, index, identifier, plot_type='square', color='#FF5500'):
		'''
		Draw a datapoint at a specific inxex (rows) for a specific identifier (columns)
		'''
		if identifier not in self.data_adjusted.keys():
			return

		x, junk = self.data_adjusted[identifier]
		y = self.index_y_coordinates[index]

		if plot_type == 'square':
			self.scene.add(wsvg.Rectangle(origin=(x-self.rect_size/2, y), height=self.rect_size, width=self.rect_size, fill_color=color, line_color=color, line_width=0.5))

		elif plot_type == 'bottom':
			self.scene.add(wsvg.PolyLine(points=[(x+self.rect_size/2, y+self.rect_size),(x-self.rect_size/2, y+self.rect_size), (x-self.rect_size/2, y)], closed=True, fill_color=color, line_color=color, line_width=0.5))

		elif plot_type == 'top':
			self.scene.add(wsvg.PolyLine(points=[(x-self.rect_size/2, y), (x+self.rect_size/2, y), (x+self.rect_size/2, y+self.rect_size)], closed=True, fill_color=color, line_color=color, line_width=0.5))

		else:
			raise ValueError


	def _draw_row_label(self, index, text, color):
		'''
		Add a label at a specific row index
		'''
		x = self.dend_min_x
		y = self.index_y_coordinates[index]
		self.scene.add(wsvg.Text(text=text, origin=(x-self.rect_size, y+self.fontsize), angle=0, size=self.fontsize, weight="normal", color=color, anchor="end", align='end'))


	def _assign_property_colors(self, prop):
		'''
		Automatically get property colors if needed
		'''
		# collect all unique values for this property
		all_vals = set([])
		for uid in self.property_dict.keys():
			all_vals.add(self.property_dict[uid][prop])

		# count their number and get the equivalent number of colors
		all_vals = sorted(all_vals)
		num_vals = len(all_vals)
		colors = scales.rainbow(num_vals, start_col='#ADFF2F', clockwise=True, full=True)[0::2]
		colors.extend(scales.rainbow(num_vals, start_col='#6EA31E', clockwise=True, full=True)[1::2])
		new_col_dict = {k:v for k, v in zip(all_vals, colors)}

		# add to the color dictionary
		self.property_colors[prop] = new_col_dict


	def _draw_seq_info(self):
		'''
		Draw sequence information that is in addition to the substrate activity data.
		'''
		color = '#111111'
		uids = list(self.property_dict.keys())
		property_types = sorted(self.property_dict[uids[0]].keys(), reverse=True)
		#self.all_properties = sorted(self.property_dict[uids[0]].keys())


		# add the property row labels
		i = 0
		for prop in property_types:
			# do I have pre-assigned colors? If not, automatically assign.
			if self.property_colors.get(prop) is None:
				self._assign_property_colors(prop)

			self._draw_row_label(index=i, text=prop, color=color)
			self.indexes_used.append(i)
			i += 1

		# now draw the data
		i = 0
		for prop in property_types:
			for uid in sorted(self.property_dict.keys()):
				value = self.property_dict[uid][prop]
				color = self.property_colors[prop][value]
				self._draw_datapoint(index=i, identifier=uid, color=color, plot_type='square')
			i += 1


	def _draw_substrate_activity_info(self):
		'''
		Add activity information as a kind of heatmap.
		'''
		# get a list of all substrates
		plt_type = 'square'
		self.all_subst = sorted(self.subst_data[list(self.subst_data.keys())[0]].keys())

		# assign an index for these
		self.subst_index = {k:v for k,v in zip(self.all_subst, range(len(self.indexes_used)+2, len(self.all_subst)+len(self.indexes_used)+2))}

		# assign a color for each substrate
		subst_color = {k:v for k,v in zip(self.all_subst, scales.rainbow(len(self.all_subst)))}

		# add the substrate row labels
		for subst in self.all_subst:
			self._draw_row_label(index=self.subst_index[subst], text=subst, color=subst_color[subst])

		# is there a repeat that should be plotted as well?
		if self.subst_data_replicate is not None:
			for uid in sorted(self.subst_data_replicate.keys()):
				for subst in sorted(self.subst_data_replicate[uid]):
					if self.subst_data_replicate[uid][subst] is True:
						self._draw_datapoint(index=self.subst_index[subst], identifier=uid, color=subst_color[subst], plot_type='top')
			plt_type = 'bottom'

		# now draw the data
		for uid in sorted(self.subst_data.keys()):
			for subst in sorted(self.subst_data[uid]):
				if self.subst_data[uid][subst] is True:
					self._draw_datapoint(index=self.subst_index[subst], identifier=uid, color=subst_color[subst], plot_type=plt_type)



class alignment(_DrawingBaseClass):
	'''
	'''
	def __init__(self, data, filepath, main=None, group_dict=None, group_colors=None, text_color='#424242'):
		'''
		Take a dictionary with identifier keys and protein sequence values.
		The protein sequenes must have been aligned before.
		filepath specifies the output file
		main is optional and is a string that specifies the plot title
		group_dict is optional, it is a dictionary with identifier keys and string values that indicate group belonging
		'''
		_DrawingBaseClass.__init__(self, data=data, filepath=filepath, main=main, group_dict=group_dict, group_colors=group_colors, text_color=text_color)

		self.colors = {'R':'#8694fa', 'K':'#baaafc',
				'E':'#f93333', 'D':'#fb7979',
				'I':'#ffff4f', 'L':'#ffff79', 'V':'#fefeab', 'A':'#ffffc9',
				'C':'#e2f9ad', 'H':'#d5f6fb', 'M':'#c3ed27',
				'N':'#ee72a7', 'Q':'#f8c3e3',
				'F':'#c7c88a', 'Y':'#7dafb9', 'W':'#85b0cd',
				'S':'#ca9ec8', 'T':'#f0e4ef',
				'G':'#c0c0c0',
				'P':'#f1f2f3',
				'-':'#ffffff'} # as in UGENE}

		self.seq_len = len(self.data[self.id_list[0]])
		self.rect_size = 8 #rectangle size
		self.residues_per_line = 70
		self.number_of_blocks = math.ceil(self.seq_len / self.residues_per_line)
		self.space_between_blocks = self.rect_size * 3

		# setup x offset for the numbering
		self.left_text_offset = self.rect_size/2.0 - 3
		self.right_text_offset = self.rect_size/2.0 + self.residues_per_line*self.rect_size + 7

		# setup group vairable
		if self.usegroups is True:
			self.group_spacer = self.rect_size/2
		else:
			self.group_spacer = 0

		# setup page size
		self.page_x_size = 744.09
		self.page_y_size = self.y_spacing*2 + self.number_of_blocks*self.rect_size*self.num_sequences + self.number_of_blocks*self.space_between_blocks + self.number_of_blocks*self.group_spacer*self.number_of_groups

		# setup the svg scene
		self._start_plotting()

		# set title if present
		if main is not None:
			assert type(main) is str, 'Error, the variable "main" must be a string.'
			self._title()

	 	# go through each identifier and draw the sequence
		uid_number = 0
		for uid in self.id_list:
			uid_number += 1
			self._draw_single_sequence(uid, uid_number)

		# finish document
		self._stop_plotting()

		# convert to pdf
		self._convert_to_pdf()


	def _draw_single_sequence(self, uid, uid_number):
		'''
		'''
		block_number = -1
		residue_number = 0 # keep track of the actual number of amino acid residues (skip spaces)

		# which group in the order? First, second, third  etc...
		current_group = self.group_dict[uid]
		current_group_num = self.all_groups.index(current_group)

		# get the sequence
		seq = self.data[uid]

		for i in range(0, self.seq_len):

			# setup values for drawing
			aa = seq[i]
			rect_color = self.colors[aa]
			text = aa

			# count the amino acid residues
			if aa != '-':
				residue_number += 1

			# when n number of residues are reached, skip down one block
			if i % self.residues_per_line == 0:
				block_number += 1
				x = self.x_spacing
				y = self.y_spacing + (block_number * self.rect_size * self.num_sequences) + (self.rect_size * uid_number) + (self.space_between_blocks * block_number) + (self.group_spacer * self.number_of_groups * block_number) + (self.group_spacer * current_group_num)

				#draw the group, if applicable
				if self.usegroups is True:
					self.scene.add(wsvg.Rectangle(origin=(x+self.right_text_offset+20, y), height=self.rect_size, width=self.rect_size/2, fill_color=self.group_colors[current_group_num], line_color=self.group_colors[current_group_num], line_width=0.5))
					self.scene.add(wsvg.Text(text=current_group, origin=(self.x_spacing + self.right_text_offset + 30, y+self.fontsize*0.80), angle=0, size=self.fontsize, weight="normal", color=self.text_color, anchor="start", align='start'))

				# add the identifier
				self.scene.add(wsvg.Text(text=uid, origin=(x+self.rect_size/2.0 - 30, y+self.fontsize*0.80), angle=0, size=self.fontsize, weight="normal", color=self.text_color, anchor="end", align='end'))

				# add left numbering
				if i == 0:
					num = 1
				else:
					num = residue_number
				self.scene.add(wsvg.Text(text=num, origin=(x+self.left_text_offset, y+self.fontsize*0.80), angle=0, size=self.fontsize, weight="normal", color=self.text_color, anchor="end", align='end'))

			# add right numbering
			if (i+1) % self.residues_per_line == 0:
				self.scene.add(wsvg.Text(text=residue_number, origin=(self.x_spacing + self.right_text_offset, y+self.fontsize*0.80), angle=0, size=self.fontsize, weight="normal", color=self.text_color, anchor="start", align='start'))

			# increment x-position by 1
			x += self.rect_size

			# actually draw the box for the individual amino acid
			self.scene.add(wsvg.Rectangle(origin=(x, y), height=self.rect_size, width=self.rect_size, fill_color=rect_color, line_color=rect_color, line_width=0.5))
			self.scene.add(wsvg.Text(text=text, origin=(x+self.rect_size/2.0, y+self.fontsize*0.80), angle=0, size=self.fontsize, weight="bold", color=self.text_color, anchor="middle", align='center'))

		# add final right number
		self.scene.add(wsvg.Text(text=residue_number, origin=(x + 10, y+self.fontsize*0.80), angle=0, size=self.fontsize, weight="normal", color=self.text_color, anchor="start", align='start'))


#
# x = LoadAlignments()
# data = x.alignments()
# ec = '1.1.3.21'
# group = 0
# groups = {'U4KM46':'hit', 'Q2SSQ7':'pop', 'Q6MTY6':'pop', 'F9EP79':'tart', 'D5RC44':'tart', 'Q7NC64':'tart', 'Q70GP9':'tart'}
# alignment(data=data[ec][group], filepath='test.svg', main='%s %s' % (ec, group), group_dict=groups)
