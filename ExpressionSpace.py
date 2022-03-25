import numpy as np
import copy
import time

"""
	A class which contains all the points of the space and keep track on their
	values, if they are attractors, and the attractors basins.
"""
class ExpressionSpace:
	def __init__(self, keys):
		# The dimensions of the space.
		self.dim = len(keys)
		self.keys = keys.copy()

		# Masks all the bits used in Point.attr and Point.coords. Consists of
		# dim consecutive 1-bits.
		self.mask = (1 << self.dim) - 1

		""" Since Python cannot handle subclasses we must wrap the subclassses
			into a container which we access here.
		"""
		# Contains all the points in the space, indexed på Point.coords.
		subclasses = self._subclass_container()
		self.Point = subclasses["Point"]
		del subclasses

		self.points = [0 for x in range(1 << self.dim)]
		self.points[0] = self.Point(0,0)
		for p in self.points:
			p.setup()

		# The attractors present in this space. Each element contains one
		# dictionary containing information about an attractor. The elements
		# the dictionary are: 'points' (the attrator's Point objects),
		# 'states' (the attractor's coordinates), 'degeneracy' (how many
		# states the attractor consists of), 'inclusive' (the states inclueded
		# in the inclusive (more often called soft) basin), and 'exclusive' (the states inclueded in
		# the exclusive (more often called strics) basin).
		self.attractors = []
		self.possibleFates = {i:[] for i in range(1 << self.dim)} # The attractors possible to reach from each state.
		


	def	update(self, oldVal, newVal, key):
		self.points[0].update(oldVal, newVal, key, 0, 1)

	def update_vals(self, vals):
		self.points[0].update_vals(vals)


	"""
		Checks if a state is an attractor or degenerate attractor and adds to
		the counting of for how many configurations this state has been an
		attractor.
	"""
	def check(self):
		self.attractors.clear()
		pot_attr = []
		pot_attr_dic = {}
		pot_attr_dic_map = {}
		checked = [False for x in range(1 << self.dim)]
		for p in self.points:
			if(p.attr == 0):
				pot_attr.append([p])
				pot_attr_dic[p.coords] = p
				pot_attr_dic_map[p.coords] = len(pot_attr) - 1
		
		
		for i in range(len(pot_attr)):
			if(not pot_attr[i]):
				continue
			p = pot_attr[i][0]
			c = p.coords
			checked[c] = True
			if(p.flat == 0):
				p.attrCount += 1
				self.attractors.append({'points': copy.copy([p])})
				self.attractors[-1]['states'] = [self.attractors[-1]['points'][0].coords]
				self.attractors[-1]['degeneracy'] = 1
				continue
			
			plateau = False
			tmp = [p]
			new = 1
			while(True):
				newnew = 0
				for t in tmp[-new:]:
					for m in range(self.dim):
						if(t.flat & (1 << m) != 0):
							neigh = t.coords ^ (1 << m)
							if(not checked[neigh]):
								if(neigh in pot_attr_dic):
									tmp.append(pot_attr_dic[neigh])
									pot_attr[pot_attr_dic_map[neigh]].pop()
								else:
									tmp.append(self.points[neigh])
									plateau = True
								newnew += 1
								checked[neigh] = True
				new = newnew
				if(new == 0):
					break
			if(plateau):
				pot_attr[i] = []
			else:
				pot_attr[i] = tmp
				self.attractors.append({'points':copy.copy(tmp)})
				self.attractors[-1]['states'] = [self.attractors[-1]['points'][x].coords for x in range(len(tmp))]
				self.attractors[-1]['degeneracy'] = len(tmp)
				for t in tmp:
					t.attrCount += 1
			
		pot_attr = []	
		checked = [False for x in range(1 << self.dim)]


	def binom(self, n, k):
		if(k == n or k == 0):
			return 1
		return self.binom(n-1, k) + self.binom(n-1, k-1)


	"""
		Summarises the attractor counts for points in this space.
		* attrs: a list of size dim, which this method will fill with the
			attractors counts. Each entry is found at the index equal to the
			point's coordinates.
		* return: the maximum attractor count of any point.
	"""
	def countAttractors(self, attrs):
		if(len(attrs) != len(self.points)):
			print("Attractor count array of wrong size " + str(len(attrs))
					+ " (should be " + str(len(self.points)) + ")")
		# Ensures that the attractors created by latest update are counted
		self.check()

		maxCount = 0
		for p in self.points:
			if(p.attrCount > maxCount):
				maxCount = p.attrCount
			attrs[p.coords] = p.attrCount
		return maxCount

	"""
		Ranks all the attractors after their prevalence.
		* display: The number of attractors to display. If the number given is
			higher than the total number of states, that is used instead.
		* out: optional filehandler argument given if it is desired to brint
			the ranking to an open file.
		* return: the number of times the most prevalent attractor was an
			attractor.
	"""
	def rankAttractors(self, display=49, out=None):
		print("\nRanking attractors...\n")

		display = min(display, 1 << self.dim)

		ranked = self.points.copy()
		ranked = sorted(ranked, key=lambda a: a.attrCount, reverse=True)


		count = -1

		f1 = "%" + str(len(str(ranked[0].attrCount))) + "d : %s\n"
		f2 = "%" + str(len(str(ranked[0].attrCount))) + "s : %s\n"

		print("The %d most prevalent attractors:\n" % display)
		if(out is not None):
			out.write("\n\nThe %d most prevalent attractors:\n" % display)

		for i in range(display):
			if(ranked[i].attrCount != count):
				count = ranked[i].attrCount
				print(f1 % (count, self.ptos(ranked[i].coords)))
				if(out is not None):
					out.write(f1 % (count, self.ptos(ranked[i].coords)))
			else:
				print(f2 % ("", self.ptos(ranked[i].coords)))
				if(out is not None):
					out.write(f2 % ("", self.ptos(ranked[i].coords)))
		return ranked[0].attrCount


	def makeAttractorHistogram(self, out):
		self.check()

		maxCount = 0
		for p in self.points:
			if(p.attrCoun > maxCount):
				maxCount = p.attrCount

		format = "%s : %" + str(len(maxCount)) + "d | s\n"

		for p in self.points:
			out.write(format % (self.ptos(p.coords), p.attrCount, self.histoBar(p.arrtrCount, maxCount)))

		return maxCount


	def histoBar(self, val, max):
		height = 64 * val / max
		bar = [0 for x in range(65)]
		for i in range(height):
			bar[i] = 'X'
		dec = int(10 * height - int(height))
		if(dec != 0):
			bar[int(height)] = chr('0' + dec)
		return str(bar)


	"""
		Finds the inclusive basins for each attractor. An inclusive basin is
		defined as all the states from where it is possible to reach the
		attractor The attractors are also added to the list of possible fates
		for each state, which are stored as a dictionary.
	"""
	def findInclusiveBasins(self):
		""" Recursive sub method for finding basins """
		
		start = time.time()
		
		self.possibleFates = {i:[] for i in range(1 << self.dim)}
		for a in self.attractors:
			incl = []
			checked = [False for x in range(1 << self.dim)]
			a_coord = a['states'][0]
			incl.append(a_coord)
			checked[a_coord] = True

			self.possibleFates[a_coord].append(a_coord)
			new = 1
			while(True):
				newnew = 0
				for s in incl[-new:]:
					attr = self.points[s].attr
					flat = self.points[s].flat
					for k in range(self.dim):
						if(attr & (1 << k) == 0 or flat & (1 << k) != 0):
							s_next = s ^ (1 << k)
							if(not checked[s_next]):
								incl.append(s_next)
								self.possibleFates[s_next].append(a_coord)
								newnew += 1
								checked[s_next] = True
				new = newnew
				if(new == 0):
					break
			a['inclusive'] = incl
		
		return

	"""
		Finds the exclusive basins for each attractor. An exclusive basin is
		defined as all the states that can only reach the attractor. When
		the inclusive basins are known, the exclusive basin for attractor A
		is the set difference of inclusive basin A and the other attractors'
		inclusive basins.
	"""
	def findExclusiveBasins(self):
		mapping = {}
		ex = {}
		for a in range(len(self.attractors)):
			ex[self.attractors[a]['states'][0]] = []
			mapping[self.attractors[a]['states'][0]] = a
		for s in self.possibleFates:
			if(len(self.possibleFates[s]) == 1):
				ex[self.possibleFates[s][0]].append(s)
		for e in ex:
			self.attractors[mapping[e]]['exclusive'] = ex[e]

		return


	""" The container of subclasses. """
	def _subclass_container(self):
		_parent_class = self

		"""
			Represents the points in the space in a very particular fashion.
			While they are sequentially accessible from points, they are also
			interconnected in a directed graph whose traversal makes recursive
			manipulation of the space very efficient.

			A point are linked directly only to the points at Hamming distance
			1 from it (i.e. those coords are this points coords with exactly
			one bit flipped), but the graph is directed so that each node only
			links to the points whose coords contain more 1-bits than its own.
			Thus, all possible links exist in exactly one direction.

			Links where differing 1-bits is at a higher position than any
			1-bit in its coords are called strong links, while the rest are
			called weak links. Traversing the graph using only strong links
			ensures that each node is visited exactly once. A method typically
			uses this by comparing with weakly linked nodes an recursing on
			strongly linked ones, ensuring total coverage of the nodes and
			edges of the graph with minimal redundancy.

			The graph has one further advantage over iteration through points:
			if one uses the links in the order they are present in links, a
			point is guaranteed to be updated through strong-links recursion
			before it is touched through a weak link. This is often a critical
			feature that makes up for the recursion overhead.

			If the points are ordered numerically by coords, the ith point
			that has j 1-bits will have dim-j links, of which dim-j-i are
			strong. Thus, the point with coords==0 has dim links, all strong,
			while the point with coords==mask has no links at all.
		"""
		class Point:
			""" Builds a point, and, recursively, all points strong-linked by it. """
			def __init__(self, coords, index):
				self._ES = _parent_class # The ExpressionSpace's self-variables.

				# A pattern of bits representing the coordinates of this point.
				self.coords = coords

				# The current value of the function at this point.
				self.value = 0

				# Each bit in this integer is set to 0 if this point is stable
				# in the corresponding direction, and to 1 if it is not.
				# The point is an attractor if attr==0.
				self.attr = self._ES.mask

				# Like attr, but instead marks directions where the energy
				# around the point is flat.
				self.flat = self._ES.mask

				# Counts the number of configurations at which this point is
				# an attractor.
				self.attrCount = 0

				# Gradient in energy function
				self.grad = [0 for x in range(self._ES.dim)]
				


				""" Builds link list and finds the strong index. """
				# The coords of all points linked to by this one.
				self.links = [0 for x in range(bin(self._ES.mask ^ coords).count('1'))]
				s = len(self.links)
				j, m = 0, 1
				for i in range(self._ES.dim):
					if((m & coords) == 0):
						self.links[j] = coords | m

						if(i == index):
							s = j
						j += 1
					m <<= 1

				# Entries in links with at least this index are strong links.
				self.strong = s

				# Recursion!
				for i in range(self.strong, len(self.links)):
					index += 1
					self._ES.points[self.links[i]] = Point(self.links[i], index)

				self.indices = [0 for x in range(len(self._ES.keys))]


			"""
				Applies a new set of values of a function to the space, and
				performs the necessary comparisons.

				* oldVal: the old values supplied by the function.
				* newVal: the new values supplied by the function.
				* key: defines the function's subspace by masking only bits
					corresponding to dimensions contained in it.
				* idx: the current index used to extract values.
				* bit: the next bit to set in idx. Its use gets a little
					complicated in order to efficiently sync the function's
					subspace with global space.
			"""
			def update(self, oldVal, newVal, key, idx, bit):
				self.value += newVal[idx] - oldVal[idx] - 1 # Must take - 1 since the values in the operators are shifted by 1 (due to program technical reasons)
				m = 0
				for i in range(len(self.links)):
					m = self.coords ^ self.links[i]

					# If the bit is present in the function's subspace...
					if((m & key) !=0 ):

						# If we have a strong link, recurse with the next value.
						if(i >= self.strong):
							nextBit = bit << 1
							self._ES.points[self.coords | m].update(oldVal,
							  newVal, key, idx | bit, nextBit)
							bit = nextBit
						self.compare(self._ES.points[self.coords | m])
					else:
						# If we have a strong link, recurse with the same value.
						# This creates an exact copy of the function's subspace
						# in all dimensions outside it.
						if(i >= self.strong):
							self._ES.points[self.coords | m].update(oldVal,
							  newVal, key, idx, bit)

			def setup(self):
				for i in range(len(self._ES.keys)):
					k, m = 1, 1
					while(m < self._ES.keys[i]):
						if((m & self._ES.keys[i]) != 0):
							if((m & self.coords) != 0):
								self.indices[i] |= k
							k <<= 1
						m <<= 1


			def update_vals(self, vals):
				self.value = 0
				for i in range(len(self._ES.keys)):
					self.value += (vals[i][self.indices[i]] - 1) # Must take - 1 since the values in the operators are shifted by 1 (due to program technical reasons)

				for i in range(self.strong):
					self.compare(self._ES.points[self.links[i]])
				for i in range(self.strong, len(self.links)):
					self._ES.points[self.links[i]].update_vals(vals)
					self.compare(self._ES.points[self.links[i]])


			"""
				Compares the valuess and two neighbouring points and updates
				their stability information based on the direction of the
				gradient between them.
			"""
			def compare(self, other):
				diff = self.coords ^ other.coords
				compl = self._ES.mask ^ diff
				grad = self.value - other.value

				if(grad < 0):#self.value < other.value):
					# Self is stable, other is not
					other.attr |= diff
					self.attr &= compl

					self.flat &= compl
					other.flat &= compl
				elif(grad > 0):#other.value < self.value):
					# Other is stable, self is not
					self.attr |= diff
					other.attr &= compl

					self.flat &= compl
					other.flat &= compl
				else:
					# Equal energy, both are treated as stable
					self.attr &= compl
					other.attr &= compl
#					other.attr |= diff

					self.flat |= diff
					other.flat |= diff
				
				
				# The length of the binary string representation minus one is
				# the directional index.
				index = len(str(bin(diff)[2:])) - 1
				self.grad[index] = grad
				other.grad[index] = -grad

			def addTo(self, graph, ctr):
				level = self._ES.dim - len(self.links)
				if(graph[level] == None):
					graph[level] = ['' for x in range(self._ES.binom(self._ES.dim, level))]
				ctr[level] += 1
				graph[level][ctr[level]] = "%s(%d)%s" % (self._ES.ptos(self.coords),
													 self.value,
													 self._ES.ptos(self.attr))

				for i in range(self.strong, len(self.links)):
					self._ES.points[self.links[i]].addTo(graph, ctr)

	# Class point is done
	# End of subclass container
		return {"Point": Point}


	def toString (self):
		graph = [[] for x in range(self.dim + 1)]
		ctr = [0 for x in range(self.dim)]

		self.points[0].addTo(graph, ctr)

		sb = ''
		for G in graph:
			for g in G:
				sb += g
			sb += '\n\n'
		return sb

	def pointToString(self, point, length):
		s = ['' for x in range(length)]
		m = 1
		for i in reversed(range(length)):
			s[i] = '1' if (point & m) != 0 else '0'
			m <<= 1
		return s

	def ptos(self, point):
		return self.pointToString(point, self.dim)
















