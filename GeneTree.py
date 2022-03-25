from Operators import Operators as op
import sys
import random
import math


""" The values are defined in the Operator module. """
REP, NON, ACT = op.REP, op.NON, op.ACT

""" Values to indicate overexpression, kockdown, or lack thereof. """
OE, KD, FREE, M, A = +1, -1, 0, +2, -2


"""
	A gene structure that retains the details of how its genes are combined,
	allowing for edits to its configuration, either randomly through
	makeRandomChange/undoRandomChange, or efficiently and sequentially through
	nextConfig.
"""
class GeneTree:
	"""
		When a GeneTree is initiated, it must be provided with a list of
		operator names and a list of geneStatuses. In addition, it must also
		be provided with either the pair target and genes, or stringList
		and names. Further, a boolean flag fixed.
		* operators: a list of the names of the operators used to form logical values.
		* geneStatus: a list wich indicates if a gene is overexpressed, kocked down or free.
		* target: the index of the target genem for which this tree serves as input.
		* genes: the indices of genes acting on the target. Each gene willl be represented with a leaf in a tree.
		* stringList: a list where the first entry is the name of the target gene and the second element is the string representation of the trees configuration.
		* names: a list of all the genes in the network where the names are ordered according to the gene indices.
		* fixed: boolean value which flags if the tree should be able to change configuration or not.
	"""
	def __init__(self, operators, geneStatus, target=None, genes=None, stringList=None, names=None, fixed=False):

		self.ops = operators	# List of strings with names of operators
		self.p = len(operators) # The number of operators used by this tree.
		self.fixed = fixed # Boolean value which states if the tree has a fixed configuration or not.
		self.geneStatus = geneStatus # List where it is indicates if the genes are overexpressed, knocked down or free.


		""" Since Python cannot handle subclasses we must wrap the subclassses
			into a container which we access here.
		"""
		subclasses = self._subclass_container()
		self.Node = subclasses['Node']
		del subclasses


		""" Initiate a tree from a target gene and a list of input genes. """
		if((target is not None) and (genes is not None)):
			self.treeOrder = len(genes) # Number of genes present in the tree
			self.target = target # The index of the target tree.
			if(self.treeOrder > 0):
				self.root = self.Node(genes, 0) if self.treeOrder > 1 else self.Node(genes[0]) # The node which the rest of the tree is built.

				self.geneKey = self.root.key | (1 << self.target) # A mask with bits set in the posiitons corresponding to gene indices in the tree.

				# step = 2 ** (number of gene indices lower than target)
				self.step = 1 << (bin(self.root.key & ((1 << self.target) - 1)).count('1')) # Used for achieving the proper ordering of energy values. It is two to the power of the number of gene indices in the tree that are lower than target.

				self.newEnergy = [0 for x in range(1 << (self.treeOrder + 1))] # The latest energy values of the tree, numerically ordered by coordinates in the activation subspace covered by this tree.
				self.oldEnergy = [0 for x in range(1 << (self.treeOrder + 1))] # Contain the previous contents of newEnergy, updated every time it is generated with updateEnergy().
			else:	# No-input tree
				self.root = None

				self.geneKey = 1 << target
				self.step = 0

				self.newEnergy = [NON, NON]
				self.oldEnergy = [NON, NON]
			self.geneList = [g for g in genes] # The gene indices included in this tree.

		elif((stringList is not None) and (names is not None)):
			split = stringList
			self.target = self.linearSearch(names, split[0])
			
			self.geneList = []

			if(split[1] == 'NONE' or split[1] == 'Medium' or split[1] == 'Absent medium'): # If the target has no input.
				self.root = None

				self.geneKey = (1 << self.target)
				self.treeOrder = 0
				self.step = 0

				self.newEnergy = [NON, NON]
				self.oldEnergy = [NON, NON]
			else:
				self.root = self.rebuildTree(split[1], names)
				self.geneKey = self.root.key | (1 << self.target)
				self.treeOrder = self.root.order
				
				# step = 2 ** (number of gene indices lower than target)
				self.step = 1 << (bin(self.root.key & ((1 << self.target) - 1)).count('1'))

				self.newEnergy = [0 for x in range(1 << (self.treeOrder + 1))]
				self.oldEnergy = [0 for x in range(1 << (self.treeOrder + 1))]
			if(names is not None):
				self.names = names
		self.rand = random.random() # Used for scrambling the tree.
			


	"""
		Rebuilds a node and its children from their string representation.
		Also fills geneList with the encountered gene indices.
	"""
	def rebuildTree(self, string, names):
		idx = [0]
		left = self.extract(string, idx)

		# Builds a leaf
		if(idx[0] >= len(string)):
			if(left[0] == '~'):
				ls = self.linearSearch(names, left[1:])
				if(ls == 0):		# Special treatment needed to get a negative zero!
					gene = math.copysign(ls, -1)
				else:
					gene = -self.linearSearch(names, left[1:])
				self.geneList.append(-gene)
			else:
				gene = self.linearSearch(names, left)
				self.geneList.append(gene)
			return self.Node(gene)


		# Builds a non-leaf
		o = self.extract(string, idx)
		right = self.extract(string, idx)

		n = self.Node()

		n.opr = self.linearSearch(self.ops, o)

		n.L = self.rebuildTree(left, names)
		n.R = self.rebuildTree(right, names)

		n.restructureAndApply(o)
		return n


	"""
		Calculates the number of unique configurations of this tere, which
		should be the number of times nextConfig() can be applied before
		returning True.
	"""
	def getNumberOfConfigs(self):
		if(self.fixed):
			return 1
		if(self.treeOrder < 2):
			return 1
		if(self.treeOrder < 3):
			return self.p
		return self.p * (3 * self.p - 2)**(self.treeOrder - 2)


	"""
		Tries the next congiguration of the tree. Returns True if it was the
		last configuration, otherwise False.
	"""
	def nextConfig(self):
		done = self.root.change()
		return done


	"""
		Updates and returns the current energy values given by this tree
		acting on its target.
	"""
	def updateEnergy(self):
		self.oldEnergy = self.newEnergy.copy()

		k = 0
		for i in range(0, len(self.newEnergy), 2 * self.step):
			for j in range(self.step):
				self.newEnergy[i + j] = self.root.values[k]
				self.newEnergy[i + j + self.step] = 2 - self.root.values[k]
				k += 1


	"""
		Makes a random edit to the structure of the tree. It will not be truly
		random in the sense that all states are equally possible; the probability
		of a new state is instead very roughly the same as its distance from the
		current one in terms of applications of {@link #nextConfig()}.

		The edit will happen at or below depth n with probability 2^(-n), and
		the nature of the edit will be a permutation with probability 1/4
		(only if permutations are possible for the node in question) and an
		operator change otherwise.

	     * random: a random integer that will be used to generate the change.
	"""
	def makeRandomChange(self, random):
		self.root.changeRandom(random)

	""" Undoues the effect of the last call to makeRandomChange(). """
	def undoRandomChange(self):
		self.root.undoRandom()


	"""
		Rebuilds the tree with a completely random structure. Unlike
		makeRandomChange(), the resulting configuration is unrelated to
		the current one, and unlike nextConfig(), the result may be a
		right-heavy tree (although commutative operators ensures that all
		trees are equivalent to their mirror image).

		The uniformly random structure is created using the balanced
		parenthesis-generation algorithm of Arnold & Sleep (1980), but directly
		translated to a tree along the lines of Knuth TAOCP 7.2.1.6, adapted and
		optimised for this use.

		The genes are randomly placed by Fisher-Yates shuffling the gene index
		list and assigning them to leaves as they are put on the tree. Operator
		indices are simply randomised.
	"""
	def shuffle(self):
		random.shuffle(self.geneList)

		# The root is implicitly a left child, so we start at r = 1
		k = self.root.buildRandomTree(0, 2 * self.treeOrder - 2, iter(self.geneList))
		if(k != 0):
			sys.exit(1)

	"""
		Auxiliary method to buildRandomTree. Used to determine if the next
		parenthesis should be left or right.

		Notable special cases are prob(0, k) == 0 since there is no left
		parenthesis to close, and prob(k, k) == 1 since we must put only right
		parentheses to close the left ones before we run out.

		* r: the number of unbalanced left parentheses.
		* k: the number of parenthes left to insert.

		return the probability that the next parenthesis is right.
	"""
	def prob(self, r, k):
		return r*(r + k + 2)/(2 * k * (r+1))


	""" Writes the logical function represented by this tree as a string. """
	def toString(self, names):
		if(self.geneStatus[self.target] == OE):
			return names[self.target] + " <- OE"
		elif(self.geneStatus[self.target] == KD):
			return names[self.target] + " <- KD"
		elif(self.geneStatus[self.target] == M):
			return names[self.target] + " <- Medium"
		elif(self.geneStatus[self.target] == A):
			return names[self.target] + " <- Absent medium"
		elif(self.treeOrder == 0):
			return names[self.target] + " <- NONE"
		return names[self.target] + " <- " + self.root.toString(names)



	""" The container of subclasses. """
	def _subclass_container(self):
		_parent_class = self


		"""	The nodes that build up the tree. """
		class Node:
			"""
				A node-object can be initiated without arguments, with only
				genes or with genes and idx.
			"""
			def __init__(self, genes=None, idx=None):
				self._GT = _parent_class # The GeneTree's self-variables

				# Creates an empty node.

				# Default values
				self.latestRandom = -2 # Encodes the latest random change made to the node.
				self.per = 0 # Permutation counter.


				# Creates a leaf containing a single gene.
				# genes is a list of input genes, alternatively only an int (or float)
				# if there only is one input gene.
				if(genes is not None):

					# If there only is 1 input gene.
					if(((type(genes) is int) or (type(genes) is float)) and (idx is None)):
						gene = genes # Only 1 gene (hence no plural).
						self.order = 1 # The number of genes in the subtree rooted at this node.
						self.opr = -1 # Operator index.

						self.L = None # The left child of the node.
						self.R = None # The right child of the node.


						# If the gene index is zero, it must be dealt with seperately.
						# Then it is supposed to be a float and the zero is signed.
						if(gene == 0):
							self.indices = [int(gene), sys.maxsize] # An ordered list of all gene indices present in the subtree rooted at this node.
							if(self._GT.geneStatus[0] == FREE):
								self.key = 1 << int(gene) # A mask with bits set in the positions corresponding to gene indices in indices.
							else:
								self.key = 0
							val = ACT if(math.copysign(1, gene) > 0) else REP
						elif(gene < 0):
							self.indices = [-gene, sys.maxsize]
							if(self._GT.geneStatus[-gene] == FREE):
								self.key = 1 << -gene
							else:
								self.key = 0
							val = REP
						elif(gene > 0):
							self.indices = [gene, sys.maxsize]
							if(self._GT.geneStatus[gene] == FREE):
								self.key = 1 << gene
							else:
								self.key = 0
							val = ACT
						self.values = [NON, val] # Contains the energy values of subtree rooted at this node.



					# Recursively  constructs a maximally left-heavy tree. Adds a
					# leaf as the right child and calls itself for the left child,
					# until it reaches a leaf there as well.
					#
					# * genes: a list of gene indices from which the tree is built.
					# * idx: the current index in the array.
					elif((len(genes) > 1) and (idx is not None)):
						self.opr = 0

						self.L = ( self._GT.Node(genes[idx + 1])
									if idx >= len(genes) - 2
									else self._GT.Node(genes, idx + 1) )
						self.R = self._GT.Node(genes[idx])

						self.restructureAndApply(self._GT.ops[0])




			"""
				Randomly builds a tree. Applies the Arnold & Sleep
				algorithm (1980) to generate a random sequence of balanced
				parentheses, and immediately translates it to the equivalent tree.
				It handles energy calculations as it goes, leaving the correct energy
				in the root.

				The translation is as follows: if you get a left parenthesis,
				place the left child and recurse on it. If you get a right
				parenthesis, stop recursing, back up through the tree until the
				latest ancestor with no right child, place that child, and recurse
				on it. Which parenthesis to place next is based on prob(r, k)}.

				For optimisation, we do a kind of look-ahead recursion: before
				placing a child, we check which parenthesis it places; if it is a
				right, we place a leaf instead of recursing. Thus, we can
				conveniently delegate to the Node(int) leaf constructor.

				This method is guaranteed to uniformly generate a random sequence
				of balanced parentheses, and therefore a uniformly random tree
				structure. It works in O(order) time. (not sure thats correct any more)

				* r: the current number of unclosed left parentheses.
				* k: the current number of parentheses left to place.
				* genes: an iterator to a randomly shuffled list of gene indices
					to place in the leaves of the tree.

				return the number of parentheses left after the method finishes
			"""
			def buildRandomTree(self, r, k, genes):
				if(random.random() < self._GT.prob(r+1, k-1)):
					# L closes immediately: add () = left leaf
					k -= 2
					self.L = self._GT.Node(next(genes))
				else:
					# L does not close: add (
					self.L = self._GT.Node()

					# Recurse on L until it closes
					k = self.L.buildRandomTree(r+1, k-1, genes)

				if(k == 0):
					# Special treatment for the final leaf, lest we divide by zero
					self.R = self._GT.Node(next(genes))
				elif(random.random() < self._GT.prob(r,k)):
					# R closes immediately: add ) = right leaf
					k -= 1
					self.R = self._GT.Node(next(genes))
				else:
					# R does not close: no paren; we added the ) when L closed
					self.R = self._GT.Node()

					# Recurse on R until it closes
					k = self.R.buildRandomTree(r, k, genes)

				# When all children are done, we randomise and apply the operator
				self.opr = random.randint(0,self._GT.p-1)
				self.restructureAndApply(self._GT.ops[self.opr])

				return k


			"""
				Applies the given operator on between the node's two children
				and maps the the the indices between the different lists correctly.
			"""
			def restructureAndApply(self, o):
				self.key = self.L.key | self.R.key
				self.order = self.L.order + self.R.order
				self.values = [0 for x in range(1 << self.order)]
				contL = self.contract(self.L.key, self.key)
				contR = self.contract(self.R.key, self.key)

				for i in range(1 << self.order):
					# We exploit the identity that
					# contract(expand(i,A),B) = contract(i,contract(A,B))
					# for computational efficency.
					self.values[i] = op.operators[o][self.L.values[self.contract(i,contL)]][self.R.values[self.contract(i,contR)]]


			def contract(self, A, B):
				C = 0
				i, j = 1, 1
				while(i <= B and j <= A):
					if(B & i):
						if(A & i):
							C |= j
						j <<= 1
					i <<= 1
				return C

			def expand(self, A, B):
				C = 0
				i, j = 1, 1
				while(i <= B and j <= A):
					if(B & i):
						if(A & j):
							C |= i
						j <<= 1
					i <<= 1
				return C


			"""
				Moves the tree efficiently through all its possible structures,
				operator combinations, and leaf permutations.

				return: a signal of completion, only used internally between
				layers of recursion. When one invocation returns True, all
				layers halt and return True as well. An external call will
				never return False, so the value is of no external use.
			"""
			def change(self):
				# When leaf is reached, signal compleation
				if(self.opr < 0):
					return True
				# If the fixed flag is True, then do not change the configuration.
				if(self._GT.fixed):
					return True

				# Use next operator
				# Skip operator if it is the same as the child's and we have already
				# made a permutation (then we have already visisted this state).
				while True:
					self.opr += 1
					if(not (self.per > 0 and self.opr == self.L.opr)):
						break

				# When last operator is reached, reset operator and permute subtrees
				if(self.opr >= self._GT.p):
					self.opr = 0

					# If child is leaf, skip permutation and recurse
					if(self.L.opr < 0):
						self.restructureAndApply(self._GT.ops[self.opr])
						return True
					else:
						# Permute subtrees
						temp = self.R
						self.R = self.L.R
						self.L.R = self.L.L
						self.L.L = temp

						# Permutation requires restructuring before apply works
						self.L.restructureAndApply(self._GT.ops[self.L.opr])
						self.restructureAndApply(self._GT.ops[self.opr])

						# Permutation if all permutations are done
						self.per += 1
						if(self.per > 2):
							self.per = 0
							if(self.L.change()):
								# If child is done, we are done.
								self.restructureAndApply(self._GT.ops[self.opr])
								return True
						# Otherwise, just ensure we don't have the same operators
						elif(self.L.opr == 0):
							self.opr += 1

							# Update L only if L.change() didn't already do it
							self.L.restructureAndApply(self._GT.ops[self.L.opr])

				# Update function and signal that the method is not done
				self.restructureAndApply(self._GT.ops[self.opr])
				return False




			def changeRandom(self, random):

				# Recurse if leftmost bit is 1 and recursion is possible
				if(self.opr >= 0 and self.L.opr >= 0 and (random & 1) != 0):
					self.L.changeRandom(random >> 1)
					self.restructureAndApply(self._GT.ops[self.opr])
					self.latestRandom = -1
					return

				# Next two bits determine change
				adv = random & 0b11

				# If 0: permute subtrees if possible
				if(adv == 0):
					if(self.opr >= 0 and self.L.opr >= 0 and self.L.opr != self.opr):
						# Permute subtrees
						temp = self.R
						self.R = self.L.R
						self.L.R = self.L.L
						self.L.L = temp

						# Permutation requires restructuring before apply works
						self.L.restructureAndApply(self._GT.ops[self.L.opr])
						self.restructureAndApply(self._GT.ops[self.opr])
					# If Permutation is impossible/unnecessary, advance 4 operators
					else:
						adv = 4 if (4 % self._GT.p != 0) else 5
						self.opr = (self.opr + adv) % self._GT.p
				# Otherwise: advance 1, 2 or 3 operators
				else:
					# Ensures that operator is actually changed
					if(adv % self._GT.p == 0):
						adv += 1

						# Advance operator
						self.opr = (self.opr + adv) % self._GT.p
					self.restructureAndApply(self._GT.ops[self.opr])

					self.latestRandom = adv


			### Has not been tested! ###
			def undoRandom(self):

				# The value of latest random determines what kind of action to undo
				if(self.latestRandom == -2):
					# Hos not changed.
					print("No random changes to undo")
					# Might be better with throwe some exception
				elif(self.latestRandom == -1):	# Recurse on child
					self.L.undoRandom()
				elif(self.latestRandom == 0):	# Permute subtree
					temp = self.L.L
					self.L.L = self.L.R
					self.L.R = self.R
					self.R = temp

					# Permutation requires restructuring before apply works
					self.L.restructureAndApply(self._GT.ops[self.L.opr])
					self.restructureAndApply(self._GT.ops[self.opr])

				else:	# Changed operators
					self.opr = (self.opr - self.latestRandom) % self._GT.p

				self.latestRandom = -2
				self.restructureAndApply(self._GT.ops[self.opr])


			"""
				Recursively composes a string showing the logical function
				represented by the tree rooted at this node. If it is a leaf,
				the string will simply be the name of its gene, with a tilde
				marking repression.

				* names: the gene names. List indices must match gene indices.
				* return: gene names combined with operators, structured to
					reflect the form of the tree.
			"""
			def toString(self, names):
				if(self.opr == -1):
					return ("~" if self.values[1] == REP else "") + names[self.indices[0]]

				left = self.L.toString(names)
				if(self.L.opr >= 0):
					left = '(' + left + ')'

				right = self.R.toString(names)
				if(self.R.opr >= 0):
					right = '(' + right + ')'

				return left + ' ' + self._GT.ops[self.opr] + ' ' + right


		# End of subclass container
		return {"Node": Node}




	"""
		Extracts a token from a string. It is either everything inside a set
		of matching parentheses (excluding the parentheses themselves), or
		if no parentheses are found, everything up until the next space.

		* string: the string from which the token is to be extracted.
		* ptr: a single-element array, used like a c-style pointer.
			Its element is initially the index of the first character of the
			token (or the left parenthesis), and is left as the index of
			the first character after the space terminating the token (a
			space is assumed to follow the right parenthesis).

		* return: the extracted token.
	"""
	def extract(self, string, ptr):
		begin = ptr[0]
		end = begin + 1
		if(string[begin] == '('):
			par = 1
			while(end < len(string) and par > 0):
				if(string[end] == '('):
					par += 1
				elif(string[end] == ')'):
					par -= 1
				end += 1
			ptr[0] = end + 1
			# Now, begin and end mark the parentheses enclosing the LHS.
			# We extract everything inside them.
			return string[begin + 1: end - 1]
		else:
			while(end < len(string) and string[end] != ' '):
				end += 1
			ptr[0] = end + 1
			# No parentheses: extract everything until the next space
			return string[begin:end]


	def linearSearch(self, array, key):
		i = 0
		for a in array:
			if(key == a):
				return i
			i += 1
		return -1
