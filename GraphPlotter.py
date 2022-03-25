import os
import datetime
import math
import csv
import numpy as np

"""
	Plots high-dimensional plots from a list of energiy values of a network as
	graphs (i.e. as vertices connected by edges). The vertices are tagged with
	their coordinates (represented with filled sectors for 1 and empty sectors
	for 0) and coloured by their values. The edges are directed and weighted
	by the direction and magnitude of the slope of the polynomial along them.

	There are three ways to add an extra dimension to the graph; a combination
	is used to achieve somewhat readable graphs up to 7 dimensions.

		* The "natural" way: Just draw line in s new cardinal direction. Only
		works for the first two dimensions, since the screen is two-dimensional.

		* The "flip" way: draw another graph concentric with the firstm
		flipped inside out so that the connecting edges become as simpla as
		possible. A flip step increases the the size of the graph a lot, but
		not its density, so it is well suited when more dimensions are needed
		afterwards. Flip steps are used for dimension 3 and 4 and for 5 if we
		go all the way to 7.

		* The "shift" way: draw another graph on top of the other, rotated
		radially shifted. A shift step gives a very dense graph, and it is
		therefore suited for the final dimensions. Shift steps are used for
		the last two dimensions (more than so are practically impossible),
		starting at dimension 5.
		
	The output is a .tex file generating a Tikz figure.
"""

""" 
	Values to indicate overexpression, kockdown, or lack thereof.
	Also present medium (M) are absent medium (A). 
"""
OE, KD, FREE, M, A, ON, OFF = +1, -1, 0, +2, -2, +3, -3


MAX_DIMENSION = 7 # More dimensions results in VERY messy and large plots.

class GraphPlotter:

	"""
		Creates a graph plotter with a fixed colour range
		* PATH: PATH to where to save the printed files.
		* mini: the minimum value.
		* maxi: the maximum value.
		(min and max are standard functions and, hence, avoided.)
		If mini or maxi is not given, the min or max value of the given list
		of values to plot is used instead.
	"""
	def __init__(self, PATH, mini=None, maxi=None):
		self.PATH = PATH
		if(mini >= maxi):
			print("Minimum must be strictly less than maximum")	# Probably better with exception
		if(mini is not None):
			self.mini = mini
		if(maxi is not None):
			self.maxi = maxi
		if((mini is None) and (maxi is None)): # mini and maxi is in this case set later.
			self.adaptive = True
		else:
			self.adaptive = False
			self.norm = 0xff /(self.maxi - self.mini)

		self.max_r = 1	# used to know the maximum radius so the figure can be sized correctly
		
		""" Since Python cannot handle subclasses we must wrap the subclassses
			into a container which we access here.
		"""
		subclasses = self._subclass_container()
		self.Node = subclasses["Node"]
		del subclasses


	def setup(self, dimension):
		self.dim = dimension
		if(self.dim > MAX_DIMENSION):
			print("Only graphs of %s or fewer dimensions can be produced." % MAX_DIMENSION)

		# Sets up all nodes
		self.nodes = [0 for x in range(1 << self.dim)]
		self.nodes[0] = self.Node(0, 0)
		
		

	""" Normalises the colour scale and yields a colour for a given values. """
	def colour(self, val):
		c = int((val - self.mini) * self.norm)	# Normalise value
		c = max(0, min(0xff,c))	# Limits value to 0 - 255
		c = self.mapColour(c)	# applies mapping
		return c

	"""
		Imitates the default scale of pm3d in gnuplot. Similar to
		RGBColourScale, but goes through violet rather than blue.
	"""
	def mapColour(self, val):
		rgb = 0

		# RED
		if(val < 0x80):
			rgb |= (val * 2) << 16
		else:
			rgb |= 0xff << 16

		# GREEN
		if(val >= 0x80):
			rgb |= ((val - 0x80) * 2) << 8

		# BLUE
		if(val < 0x40):
			rgb |= val * 4
		elif(val < 0x80):
			rgb |= 0xff - (val - 0x40) * 4
		elif(val > 0xc0):
			rgb |= (val - 0xc0) * 4

		return rgb



	"""
		Creates a graph from a list of energy values (ordered in the correctly).
		* values: lis of energy values of the different states of a network.
		* filename: the graph will be drawn in TikZ code to a .tex file of
			this name.
		* info: some information that will be printed in the file header. All
			lines will automatically be commented out to avoid interfering
			with the graph.
		* geneStatus: Optional argument. A list that describes if a gene is
			overexpressed, knocked down or free. If none is supplied it is
			assumed that all genes are free.
		* mash: Optional arguments. Dictionary of genes to mash together. Each 
			key is the composition name and the corresponding value is a list of 
			the names of the genes to mash.
		* hyperplane: Optional argument. A list with a subset of specific genes 
			to plot.
	"""
	def generateGraph(self, values, filename, info, labels, geneStatus=[], mash={}, hyperplane=[], attractor_colour={}, basins={}):
		# If a list of gene statuses is provided, the unallowed states (i.e.
		# those where KD genes are turned on or those where OE genes are
		# turned off) are first marked and then removed. They are not to be
		# included in the plot.

		
		self.attractor_colour = attractor_colour
		self.basins = basins
		
		orgLabels = labels.copy()
		orgGeneStatus = geneStatus.copy()
		orgLabel_index = {orgLabels[i]: i for i in range(len(orgLabels))}
		print('Labels_ini', labels)
		
		""" Deal with OE/KD and media. """
		markedValues = [0 for x in range(len(values))]
		print('Number of marked before anything:', markedValues.count(1))
		if(len(geneStatus) > 0):
			for i in range(len(geneStatus)):
				if(geneStatus[i] == KD):
					for mask in range(len(values)):
						if((1 << i) & mask != 0):
							markedValues[mask] = 1
				elif(geneStatus[i] == OE):
					for mask in range(len(values)):
						if((1 << i) & mask == 0):
							markedValues[mask] = 1
				elif(geneStatus[i] == M):
					for mask in range(len(values)):
						if((1 << i) & mask == 0):
							markedValues[mask] = 1
				elif(geneStatus[i] == A):
					for mask in range(len(values)):
						if((1 << i) & mask != 0):
							markedValues[mask] = 1
				if(geneStatus[i] == OFF):
					for mask in range(len(values)):
						if((1 << i) & mask != 0):
							markedValues[mask] = 1
				elif(geneStatus[i] == ON):
					for mask in range(len(values)):
						if((1 << i) & mask == 0):
							markedValues[mask] = 1
		print('Number of marked after geneStatus:', markedValues.count(1))
		
		""" Deal with mashing together genes. """
		markedValues_2 = [0 for x in range(len(values))]
		if(len(mash) > 0):
			print('mash', mash)
			neg_gene = []
			for m in mash:
				for i, g in enumerate(mash[m]):
					if(g.startswith('neg')):
						mash[m][i] = g[3:]
						neg_gene.append(g[3:])
				markedValues_3 = [0 for x in range(len(values))]
				toBeMashed = []
				new_label = m
				for i in range(len(labels)):
					if(labels[i] in mash[m]):
						toBeMashed.append(i)
				c = 0
				tmp_lab = np.array(labels)
				print('tbm pre reordering', np.array(labels)[np.array(toBeMashed)])
				while(labels[toBeMashed[0]] in neg_gene):
					toBeMashed.append(toBeMashed.pop(0))
					c += 1
					if(c == len(toBeMashed)):
						print('You should not assign all genes as opposite, that the same as none is opposite...')
						break
				print('tBM', np.array(labels)[np.array(toBeMashed)], 'neg_genes', neg_gene)
				for t in range(len(values)):
					good = True
					foam = (1 << toBeMashed[0]) & t
					for i in toBeMashed[1:]:
						if(labels[i] in neg_gene):
							if((foam == 0 and (1 << i) & t == 0) or (foam != 0 and (1 << i) & t != 0)):
								good = False
								break
						else:
							if((foam == 0 and (1 << i) & t != 0) or (foam != 0 and (1 << i) & t == 0)):
								good = False
								break
					if(not good):
						markedValues_2[t] = 1
						markedValues_3[t] = 1
				print('Number of marked after mashing '+m+':', markedValues_3.count(1))
			for m in mash:
				first = True
				labelsToDelete = []
				statusesToDelete = []
				for i in range(len(labels)):
					if(labels[i] in mash[m]):
						if(first):
							labels[i] = m
							if(len(geneStatus) > 0):
								geneStatus[i] = 0	# Just don't mash genes you KD/OE, it's not wise!
							first = False
						else:
							labelsToDelete.append(i)
							statusesToDelete.append(i)
				labelsToDelete.sort(reverse=True)
				statusesToDelete.sort(reverse=True)
				for i in labelsToDelete:
					del(labels[i])
				if(len(geneStatus) > 0):
					for i in statusesToDelete:
						del(geneStatus[i])
			label_index = {labels[i]: i for i in range(len(labels))}
				

		if(len(geneStatus) > 0 or len(mash) > 0):
			markedValues = [markedValues[x] + markedValues_2[x] for x in range(len(values))]
			print('Number of marked totally:', markedValues.count(1) + markedValues.count(2))
			reducedValues = [0 for x in range(len(values) - (markedValues.count(1) + markedValues.count(2)))]
			red_labels = []
			unmarkedValues =  np.where(np.array(markedValues) == 0)[0]
			for x in range(len(labels)):
				if(geneStatus[x] == FREE):
					red_labels.append(labels[x])
			red_label_index = {red_labels[i]: i for i in range(len(red_labels))}
			const_val = 0
			
			map_red_to_org = np.arange(len(values))[unmarkedValues]
			map_org_to_red = {map_red_to_org[i]: i for i in range(len(map_red_to_org))}
			print(map_red_to_org)
			print(map_org_to_red)
			if(len(self.basins) > 0):
				reducedBasins = {}
				for new, old in enumerate(map_red_to_org):
					print(self.basins[old])
					reducedBasins[new] = [map_org_to_red[x] for x in self.basins[old] if x in map_red_to_org]
				self.basins = reducedBasins
			reducedValues = list(np.array(values)[unmarkedValues])	
			values = reducedValues
		
		s = bin(len(values))[2:]
		dimension = len(s) - len(s.rstrip('0'))
		
		
		


		""" Deal with hyperplanes. """
		non_disp_stat = {}
		if(len(hyperplane) > 0):
			print('Hyperplane', hyperplane)
			non_disp_stat = hyperplane.pop()
			print(non_disp_stat)
			markedValues_of = [0 for x in range(len(values))]
			HP_num = []
			HP_lab = []
			reduced_labels = []
			for i in range(len(labels)):
				if(labels[i] in hyperplane):
					reduced_labels.append(labels[i])
				else:
					HP_num.append(i)
			print('HP_num', HP_num)
			print('red_lab', reduced_labels)
			for i in HP_num:
				if(non_disp_stat[labels[i]] == 0):
					for mask in range(len(values)):
						if((1 << i) & mask != 0):
							markedValues[mask] = 1
				elif(non_disp_stat[labels[i]] == 1):
					for mask in range(len(values)):
						if((1 << i) & mask == 0):
							markedValues[mask] = 1
			reducedValues = []
			for i in range(len(values)):
				if(markedValues[i] == 0):
					reducedValues.append(values[i])
			values = reducedValues
			labels = reduced_labels
			
			# Change colour scale according to new landscape
			self.mini = min(values)
			self.maxi = max(values)
			self.norm = 0xff /(self.maxi - self.mini)
			
			s = bin(len(values))[2:]
			dimension = len(s) - len(s.rstrip('0'))
			
		print('dim', dimension)
		print('val', values)


		if(len(values) != 1 << dimension):
			print("The size of the value array must be a power of two!")
		self.setup(dimension)
		print('Dim', dimension)

		for i in range(len(self.nodes)):
			self.nodes[i].value = values[i]
		if(self.adaptive):
			self.mini = min(values)
			self.maxi = max(values)
			self.norm = 0xff/(self.maxi - self.mini)

		print('labels', labels)
		print(values[0b11], values[0b01], values[0b00], values[0b10])
		self.generateGraph_part2(filename, info, labels, geneStatus, non_disp_stat)

	""" Gnerates the TikZ-document. """
	def generateGraph_part2(self, filename, info, labels, geneStatus=[], non_disp_stat={}):

		self.sector = 360 / self.dim

		# Prints graphic
		out = open(self.PATH+filename, 'w')

		self.printHeader(filename, out)

		if(info != 0):
			for s in info.split(sep='\\R'):
				out.write('%%% ' + s)
			out.write('\n\n')

		# LaTeX document (Only for testing, better to write just TikZ-code so
		# that the file can be imported as a figure in a LaTeX-document.)
		out.write("\\documentclass[class=minimal,tikz=true]{standalone}\n")
		out.write("\\usepackage[utf8]{inputenc}\n")
		out.write("\\tikzset{>=latex}\n\n")

		out.write("\\begin{document}\n")
		out.write("\\begin{tikzpicture}\n\n")
		

		# Marks node coordinates
		
		# Made a change so I could clip the figure before writing out all the nodes.
		nodes_to_print = ''
		for n in self.nodes:
#			n.markNode(out)
			nodes_to_print += n.markNode(out)
		
		print('max_r', self.max_r)
		out.write("\\clip (0:0) circle ({:.1f});\n".format(self.max_r*1.05))
		out.write(nodes_to_print)
		
		# Print edges
		self.nodes[0].printEdges(out)
		# Print nodes on top of edges
		for n in self.nodes:
			n.printNode(out)
		for attr in self.attractor_colour:
			if(len(self.basins) > 0):
				n_att = len(self.basins[attr])
				for i, a in enumerate(self.basins[attr]):
					start = 0 + i  * (360 / n_att)
					stop = 0 + (i + 1)  * (360 / n_att)
					out.write('\n\\draw [%s, line width=0.9] (%s) +(%f:.155) arc (%f:%f:.155);' % (
									self.attractor_colour[a], ('{:0>' + str(self.dim) + '}').format(bin(attr)[2:]),
									start, start, stop))
			else:
				out.write('\n\\draw [%s, line width=0.9] (%s) circle [radius = .155];' % (
									self.attractor_colour[attr], ('{:0>' + str(self.dim) + '}').format(bin(attr)[2:])))
		
		out.write("\n\\end{tikzpicture}\n")
		out.write("\\end{document}")
		out.close()
		
		print('labels', labels)
		self.generateLegend(filename, labels, geneStatus, non_disp_stat)


	"""
		Generate the cake-styled legend. Each piece of cake represents a gene.
		If a gene is OE/KD, it is written out seperately below the cake-diagrame
	"""
	def generateLegend(self, filename, labels, geneStatus=[], non_disp_stat={}):
		print('labels', labels)
		print('non_disp_stat', non_disp_stat	)
		cakeLabels = []
		OEKD_labels = []
		if(len(geneStatus) > 0 and len(non_disp_stat) >= 0):
			for i in range(len(geneStatus)):
				if(labels[i] not in non_disp_stat):
					if(geneStatus[i] == FREE):
						cakeLabels.append(labels[i])
					elif(geneStatus[i] == KD):
						OEKD_labels.append(labels[i] + '-KD')
					elif(geneStatus[i] == OE):
						OEKD_labels.append(labels[i] + '-OE')
					elif(geneStatus[i] == M):
						OEKD_labels.append(labels[i] + '-Medium')
					elif(geneStatus[i] == A):
						OEKD_labels.append(labels[i] + '-Absent')
					elif(geneStatus[i] == ON):
						OEKD_labels.append(labels[i] + '-On')
					elif(geneStatus[i] == OFF):
						OEKD_labels.append(labels[i] + '-Off')
				else:
					OEKD_labels.append(labels[i] + '=' + non_disp_stat[labels[i]])
		else:
			cakeLabels = labels
			if(len(non_disp_stat) > 0):
				for l in non_disp_stat:
					OEKD_labels.append(l + '=' + str(non_disp_stat[l]))
			
		print('cake', cakeLabels)
		print('OEKD', OEKD_labels)

		if(len(cakeLabels) != self.dim):
			print("Wrong number of labels. Got %d dimensions and %d labels" % (self.dim,  len(cakeLabels)))
			print('The labels are', cakeLabels)
			print(geneStatus)

		out = open(self.PATH + filename.replace('.tex', '-legend.tex'), 'w')

		self.printHeader(filename.replace('.tex', '-legend.tex'), out)
		out.write("%%% This serves as a legend for a graph,\n")
		out.write("%%% which most likely is in a file with the same\n")
		out.write("%%% name as this one, save the \"-legend\".\n")
		if(len(OEKD_labels) > 0):
			out.write("%%% The following genes are over expressed (OE) / knocked downed (KD) / present medium / absent medium: \n")
			out.write("%%% " + str(OEKD_labels) + " \n")
		out.write("\n\n")

		# LaTeX document
		out.write("\\documentclass[class=minimal,tikz=true, border=1pt]{standalone}\n")
		out.write("\\usepackage[utf8]{inputenc}\n")
		out.write("\\tikzset{>=latex}\n\n")

		out.write("\\begin{document}\n")
		out.write("\\begin{tikzpicture}\n\n")
		
		node_r = str(float(self.NODE_R) * 2)
		if(self.dim <= 1):
			out.write("\\draw [black, fill=black] (0,0) circle [radius=%s];\n"
						% node_r)
			out.write("\\draw (.5,0) node {%s};\n" % cakeLabels[0])
		else:
			angle = 90 + (self.sector / 2)
			fill = True
			scale_cut = 0.7
			max_colour = self.colour(self.maxi * scale_cut)
			min_colour = self.colour(self.mini * scale_cut)
			
			
			for i in range(self.dim):
				# Draws the "cake" coordinate representation
				fill = not fill
				out.write("\\definecolor{colour}{HTML}{%06X};\n" % (min_colour if fill else max_colour))
				out.write(("\\draw [black, fill=%s] (0,0)"
							"-- (%.4f:%s) arc (%.4f:%.4f:%s);\n") %
							(("colour"),
							angle, node_r, angle,
							angle - self.sector, node_r))
				out.write(("\\draw [black] (%.4s:%s) -- (%.4f:0.75) node "
							"[at end, sloped, %s] {\\makebox{%s}};\n") %
							(angle - self.sector/2, node_r,
							angle - self.sector/2,
							('right' if ((angle - self.sector/2) >= -90 and
							(angle - self.sector/2) <= 90 )else 'left'),
							cakeLabels[i]))
				angle -= self.sector
			# Surrounds it with a black circle to make it easier to see
			out.write("\\draw [black] (0,0) circle [radius = %s];\n" % node_r)

		# The labels of OE/KD genes are separated from the cake labels.
		if(len(OEKD_labels) > 0):
			l = 0
			if(len(geneStatus) > 0):
				for i in range(len(labels)):
					if(geneStatus[i] == OE):
						out.write("\\draw [black] (0,%s) node []{\\makebox{%s - OE}};\n" % (-1.5 - l, labels[i]))
						l += 0.5
					elif(geneStatus[i] == KD):
						out.write("\\draw [black] (0,%s) node []{\\makebox{%s - KD}};\n" % (-1.5 - l, labels[i]))
						l += 0.5
					elif(geneStatus[i] == M):
						out.write("\\draw [black] (0,%s) node []{\\makebox{%s - Medium}};\n" % (-1.5 - l, labels[i]))
						l += 0.5
					elif(geneStatus[i] == ON):
						out.write("\\draw [black] (0,%s) node []{\\makebox{%s - On}};\n" % (-1.5 - l, labels[i]))
						l += 0.5
					elif(geneStatus[i] == OFF):
						out.write("\\draw [black] (0,%s) node []{\\makebox{%s - Off}};\n" % (-1.5 - l, labels[i]))
						l += 0.5
			if(len(non_disp_stat) > 0):
				for i in range(len(non_disp_stat)):
					out.write("\\draw [black] (0,%s) node []{\\makebox{%s}};\n" % (-1.5 - l, OEKD_labels[i]))
					l += 0.5
					

		out.write("\n\\end{tikzpicture}\n")


		# Generate colour bar.
		out.write("\\begin{tikzpicture}\n\n")

		w = 0.5 # Width of colourbar
		H = 3 # Height of colourbar
		c = 100 # Number of compartments
		h = H / c # Height of compartments

		COLOURSCALE = 0.9
		K = (self.maxi - self.mini) / H
		for i in range(c):
			col = self.colour((i * h * K +self.mini) * COLOURSCALE + (1 - COLOURSCALE) * self.mini)
			out.write("\\definecolor{colour}{HTML}{%06X};\n" % col)
			out.write("\\draw [colour, fill](0,%4f) rectangle (%4f,%4f);\n"
			 % (i * h, w, (i + 1) * h))


		# Writes out the scale of the colour bar. Only uses integers as ticks.
		# Dynamically uses 3 to 7 ticks depending on how the range can be divided.
		mm = self.maxi - self.mini
		if(mm <= 6):
			k = mm
			r = 1
		else:
			r = 0
			for i in range(6,2,-1):
				if(mm % i == 0):
					k = i
					r = mm / i
					break
			if(r == 0):
				if(mm % 2 == 0):
					k = 2
					r = mm / 2
				else:
					k = 2
					r = mm / 2


		for i in range(k, - 1, -1):
			out.write("\\draw [black] (%4f,%4f) node [] {\\makebox{$%s$}};\n"
					%(2* w, H / k * i, int(self.mini + i * r)))
		out.write("\\draw [black, anchor=center] (%4f,%4f) node [rotate=90] {\\makebox{Energy}};\n"
					%(-w / 2, H / 2))


		out.write("\n\\end{tikzpicture}\n\n")


		out.write("\\end{document}")
		out.close()

	""" Prints out the file header. """
	def printHeader(self, filename, out):
		out.write("%%% In file \"" + self.PATH + filename + "\"\n")
		out.write("%%%\n")
		out.write("%%% This file was generated by a graph-drawing module of CELLoGeNe" +
						" written and developed by Mattias Sjö and Emil Andersson.\n")
		out.write("%%% Graph created " + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + "\n")
		out.write("")


	# Constants used in the code to come.
	NODE_R = ".14"
	def NODE_R_add(self, add=0):
		return str(float(self.NODE_R) + add)
	THICKNESS_FACTOR = 0.6
	MAX_THICKNESS = 2

	KEY_1D = 0b0000001
	KEY_2D = 0b0000010
	KEY_3D = 0b0000100
	KEY_4D = 0b0001000
	KEY_5D = 0b0010000
	KEY_6D = 0b0100000
	KEY_7D = 0b1000000

	UPHILL = "<-"
	DOWNHILL = "->"
	FLAT = "dashed"


	""" Since Python cannot handle subclasses we must wrap the subclassses
		into a container.
	"""
	def _subclass_container(self):
		_parent_class = self

		"""
			Internal representation of graph nodes. Nodes form the "node tree"
			where each node map to those whose keys are subsets of its own,
			restricted in a systematic fashion so that no node has more than
			one parent. Thus, the connections of the tree allows the graph to
			be completely and non-redundantly traversed when drawing edges.
		"""
		class Node:
			"""
				Recursive constructor. Call on the root (zero key) to build
				entire node tree representing the graph.
				* key: the key that encodes the point this node represents.
				* index: the index of the bit that was edited in the key of
					the caller to create the key of this node.
			"""
			def __init__(self, key, index):
				# The GraphPlottes self-variables.
				self._GP = _parent_class

				# Key coordinates at this node.
				self.key = key

				# Indicess in the node-tree rooted at this node.
				self.sub = [0 for x in range(self._GP.dim - index)]

				# Functions value at this node.
				self.value = 0

				# TikZ tag for this node, representes key in binary.
				self.tag = ("(0)" if self._GP.dim == 0
							else ('({:0>' + str(self._GP.dim) + '})').format(str(bin(key))[2:]))

				for i in range(index, self._GP.dim):
					k = self.key ^ (1 << i)
					self.sub[i - index] = k
					self._GP.nodes[k] = self._GP.Node(k,i+1)

				# TikZ polar coordinates of this node.
				self.location = self.mapLocation()


			"""
				Takes a value, adds it to the node value, and passes it on to
				other nodes that should also take the value. This is used to
				convert polynomial coefficients to values: the nodes that
				accept a coefficient are exactly the nodes at whose
				coodinates the coefficient contributes to the polynomial.
				When all coefficients have been offered to the root of the
				node tree.
				* k: the key of the coefficient.
				* v: the value of the coefficient.
			"""
			# Never used!
			def accept(self, k, v):
				# The all-1's node should take all values, and the root should only
				# take one of them. In fact, value passing uses the entirely
				# opposite structure to all other graph traversals, so we simply
				# apply all changes to the node with the complemental key, and
				# everything turns out nicely.

				compl = (( 1 << self._GP.dim) - 1) ^ self.key
				self._GP.nodes[compl].value += v

				for s in self.sub:
					if((k & s) == k):
						self._GP.nodes[s].accept(k,v)





			"""
				Writes the node as a TikZ node specification. The visuals are
				created by calling printNode(out)
				* out: a printer to the TikZ file.
			"""
			def markNode(self, out):
				return "\\node %s at %s {};\n" % (self.tag, self.location)


			"""
				Prints the node to TikZ. The node is represented as a circle
				divided into dim sectors, where each sector is filled only if
				the corresponding bit in key is high (gene is turned on). The
				bit is clockwise around the circle, starting at noon. The fill
				colour is determined from the value.
				* out: a printer to the TikZ file.
			"""
			def printNode(self, out):
				COLOURSCALE = 0.9
				# Finds colour based on value
				# The colour scale is adapted so that the minimum value is
				# black, but even though the original colour scales maximum
				# value is white, the scale is altered such that white is
				# never used due to visual reasons.
				self.col = self._GP.colour(self.value * COLOURSCALE + (1 - COLOURSCALE) * self._GP.mini)

				# Defines the proper colour
				out.write("\n\\definecolor{colour}{HTML}{%06X};\n" % self.col)
				
				if(len(self._GP.basins) > 0):
					n_att = len(self._GP.basins[self.key])
					for i, a in enumerate(self._GP.basins[self.key]):
						start = 0 + i  * (360 / n_att)
						stop = 0 + (i + 1)  * (360 / n_att)
						out.write('\\draw [%s, line width=0.9] (%s) +(%.4f:.155) arc (%.4f:%.4f:.155);\n' % (
										self._GP.attractor_colour[a], ('{:0>' + str(self._GP.dim) + '}').format(bin(self.key)[2:]),
										start, start, stop))
				
				# Draws the node
				if(self._GP.dim <= 1):
					out.write("\\draw [black, fill=%s] %s circle [radius=%s];\n" %
							(("colour" if ((self.key & 1) == 1 or self._GP.dim == 0) else "white"),
							self.location, self._GP.NODE_R))
				else:
					angle = 90 + self._GP.sector / 2
					m = 1
					
					white_pie = [False if (self.key & i) == i else True for i in [1<<j for j in range(self._GP.dim)]]
					
					for i in range(self._GP.dim):
						# Draws the "cake" coordinate representation
						if(self.key != 0 ):
							out.write(("\\draw [%s, fill=%s, line width=%f] %s"
								+ "-- +(%.4f:%s) arc (%.4f:%.4f:%s);\n") %
								("black", 
								"colour" if (self.key & m) == m else "white",
								0.4,
								self.location, angle, self._GP.NODE_R, angle,
								angle - self._GP.sector, self._GP.NODE_R))
						else:
							out.write(("\\draw [%s, fill=%s, line width=%.4f] %s"
								+ "-- +(%.4f:%s) -- +(%.4f:%s);\n") %
								("colour", # Black cake markers for all the nodes except the 0-node which has the corresponding colour value.
								"colour" if (self.key & m) == m else "white",
								1.0, # The 0-node gets thicker cake markers for higher visual clarity.
								self.location, angle, self._GP.NODE_R, angle,
								self._GP.NODE_R))
						if((self.value - self._GP.mini) / (self._GP.maxi - self._GP.mini) <= 0.15):
							out.write(("\\draw [%s, line width=%f] %s"
								+ "-- +(%.4f:%s) arc (%.4f:%.4f:%s);\n") %
								("black" if (white_pie[i] and white_pie[i-1]) else "white",
								0.4,
								self.location, angle, self._GP.NODE_R, angle,
								angle - self._GP.sector, self._GP.NODE_R))
						m <<= 1
						angle -= self._GP.sector

					# Surrounds it with a black circle to make it easier to see.
					out.write("\\draw [black, line width=%f] %s circle [radius = %s];\n" %
								(0.4 if (self.key != 0) else 0.4, # The 0-node gets thicker cake markers for higher visual clarity.
								self.location, 
								self._GP.NODE_R
								))


			"""
				Prints the edged originating at this node. Uses the node to
				tree to avoid redundant edges.
				* out: a printer to the TikZ file.
			"""
			def printEdges(self, out):

				# Iterates over all keys that are supersets with Hamming distance 1
				m = 1
				for i in range(self._GP.dim):
					if((self.key & m) == 0):
						s = self.key ^ m
						if(self.key in self._GP.attractor_colour):
							attr = ', ' + self._GP.attractor_colour[self.key]
						elif(s in self._GP.attractor_colour):
							attr = ', ' + self._GP.attractor_colour[s]
						else:
							attr = ''
							
						out.write(("\\draw[%s%s] " + self.mapEdge(s) + ";\n") %
									(self.mapLine(self._GP.nodes[s].value - self.value), attr))
					m <<= 1

				# Recurse on all unique such keys
				for s in self.sub:
					self._GP.nodes[s].printEdges(out)


			"""
				Calculates the absolute lacoation of a node.
				* return: a TikZ polar coordinate specification for the node,
					on the form angle:radius
			"""
			def mapLocation(self):
				return "(" + self.mapAngle() + ":" + self.mapRadius() + ")"

			"""
				Utility for mapLocation().
				* return: the polar angle of a node's lacation.'
			"""
			def mapAngle(self):
				a = 0
				# Handles the natural steps
				if(self.key & 0b11 == 0b00):
					a = 225
				elif(self.key & 0b11 == 0b01):
					a = 315
				elif(self.key & 0b11 == 0b10):
					a = 135
				elif(self.key & 0b11 == 0b11):
					a = 45

				# Handles the shift steps
				if(self._GP.dim < 5):
					return str(a)
				elif(self._GP.dim < 7):
					if((self.key & self._GP.KEY_5D) != 0):
						a += 45
					if((self.key & self._GP.KEY_6D) != 0):
						a += 22.5
				else:
					if((self.key & self._GP.KEY_6D) != 0):
						a += 45
					if((self.key & self._GP.KEY_7D) != 0):
						a += 22.5
				return str(a)


			"""
				Utility for mapLocation().
				* return: the polar radius of a node's lacation.
			"""
			def mapRadius(self):
				r = 1

				# Handles the flip steps
				if((self.key & self._GP.KEY_3D) != 0):
					r = 3 - r
				if((self.key & self._GP.KEY_4D) != 0):
					r = 5 - r

				# Handles the shift steps
				if(self._GP.dim < 7):
					if((self.key & self._GP.KEY_5D) != 0):
						r += 0.5
					if((self.key & self._GP.KEY_6D) != 0):
						r += 0.25
				else:
					# Extra flip step
					if((self.key & self._GP.KEY_5D) != 0):
						r  = 9 - r

					if((self.key & self._GP.KEY_6D) != 0):
						r += 0.5
					if((self.key & self._GP.KEY_7D) != 0):
						r += 0.25

					# For 7-cubes, we widen the central hole and compress the rest.
					# r = 0.5 * (r + 1)
				
				if(r > self._GP.max_r):
					self._GP.max_r=r
				
				return str(r)


			"""
				Maps an energy difference to a line decoration.
				* diff: the energy difference.
				* return: a TikZ parth argument, to be used like \path[<argument>]...
			"""
			def mapLine(self, diff):
				if(diff == 0):
					return self._GP.FLAT

				if(self._GP.dim < 7):
					return("%s, line width = %.2fpt" %
						((self._GP.UPHILL if diff > 0 else self._GP.DOWNHILL),
						abs(diff)*self._GP.THICKNESS_FACTOR))
				# For the 7th dimension, the line thickness are mapped
				# through a tanh function in order to make the plot less dense.
				else:
					return("%s, line width = %.2fpt" %
						((self._GP.UPHILL if diff > 0 else self._GP.DOWNHILL),
						math.tanh(abs(diff)/self._GP.maxi)*self._GP.MAX_THICKNESS))

			"""
				Calculates the parameters for an edge connecting two nodes.
				* tkey: the key of the target node.
				return: a format string for a TikZ path specification, with
				two %s slots for the source an destination nodes.
			"""
			def mapEdge(self, tkey):
				if(self._GP.dim <= 1):
					return "%s to %s" % (self.tag, self._GP.nodes[tkey].tag)

				if(self.key ^ tkey == self._GP.KEY_1D):	# 1st dimension: draw around circle
					return(("%s to[bend "
							+ ("right" if ((self.key & (self._GP.KEY_1D|self._GP.KEY_2D)) == 0) else "left")
							+ " = 43] %s") % (self.tag, self._GP.nodes[tkey].tag))

				elif(self.key ^ tkey == self._GP.KEY_2D):	# 2nd dimension: draw around circle
					return(("%s to[bend "
							+ ("left" if ((self.key & (self._GP.KEY_1D|self._GP.KEY_2D)) == 0) else "right")
							+ " = 43] %s") % (self.tag, self._GP.nodes[tkey].tag))

				elif(self.key ^ tkey == self._GP.KEY_3D):	# 3rd dimension: draw radially
					return("%s to %s" % (self.tag, self._GP.nodes[tkey].tag))

				elif(self.key ^ tkey == self._GP.KEY_4D):	# 4th dimension: draw bent line if necessaty
					if((self.key & self._GP.KEY_3D) == 0):
						return(("%s to[bend "
								+ ("right" if ((self.key & 0b10000) != 0) else "left")
								+ " = "
								+ ("45" if(self._GP.dim < 7) else "20")
								+ " ] %s") % (self.tag, self._GP.nodes[tkey].tag))
					else:
						return "%s to %s" % (self.tag, self._GP.nodes[tkey].tag)

				elif(self.key ^ tkey == self._GP.KEY_5D):	# 5th dimension: draw between circles
					if(self._GP.dim < 7):
						return(("%s to[bend "
								+ ("right" if self._GP.dim < 6 else "left")
								+ " = 20] %s") % (self.tag, self._GP.nodes[tkey].tag))
					# Special version for 7-cube (form "rainbow")
					else:
						# if high 4th and low 3rd, draw straicht, otherwise bend
						if((self.key & self._GP.KEY_4D) != 0 and (self.key & self._GP.KEY_3D) == 0):
							return "%s to %s" % (self.tag, self._GP.nodes[tkey].tag)
						else:
							return "%s to [bend right = 15] %s" % (self.tag, self._GP.nodes[tkey].tag)

				elif(self.key ^ tkey == self._GP.KEY_6D):	# 6th dimension: draw between circles
					if(self._GP.dim < 7):
						return "%s to [bend right = 10] %s" % (self.tag, self._GP.nodes[tkey].tag)
					else:
						# To avoid collisisons in 7D, some lines need adjustment
						return(("%s to[bend left = "
								+ ("20" if ((self.key & (self._GP.KEY_4D|self._GP.KEY_5D)) != (self._GP.KEY_4D|self._GP.KEY_5D))
										else "10")
								+ "] %s") % (self.tag, self._GP.nodes[tkey].tag))

				elif(self.key ^ tkey == self._GP.KEY_7D):	# 7th dimension: draw between circles
					return "%s to [bend right = 10] %s" % (self.tag, self._GP.nodes[tkey].tag)

				else:
					print("Unvalid")




		# End of subclass container.
		return {"Node": Node}


if(__name__ == "__main__"):

	land = [0, 0, 0, 2, 2, 2, 0, 2, 1, -1, -3, -1, 3, 1, -1, 1, 1, 0, -1, 2, -1, -2, -4, -1, 1, -3, -2, 2, 1, -3, -1, 3]
	labels = ['A', 'B', 'C', 'D', 'E']

	
	# If the landscape is stored in a csv file.
	'''
	landscape_file = ''
	with open(landscape_file, 'r') as f:
		reader = csv.reader(f, skipinitialspace=True)
		land = list(reader)
	land = [int(x) for x in land[0]]
	#'''


	geneStatus = []	# List with genestatuses, i.e. if genes are OE/KD, medium, abscent medium or free.
	hyperplane = []	# Hyperplane of the hypercube to plot.
	mash = {}		# If several genes should be represented by the same piece of cake. Recommended to only do this to genes that behaves similarly for the states of interest.
	
	PATH = os.getcwd()	# PATH where generated .tex file will be stored.
	out_file = '/LandscapePlot_demo.tex'
	
	gp = GraphPlotter(PATH, min(land), max(land))
	gp.generateGraph(land, out_file,
		 '', labels.copy(), geneStatus=geneStatus.copy(), hyperplane=hyperplane, mash=mash)
