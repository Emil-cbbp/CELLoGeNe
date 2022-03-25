""" Imports """
import GeneTree as GT
from Operators import Operators
from RunData import RunData
import numpy as np
import random
import time
import math
import copy
import datetime
import sys


sys.setrecursionlimit(15000)


""" Values used to represent how genes affect each other and their status. """
REP, NON, ACT = Operators.REP, Operators.NON, Operators.ACT
OE, KD, FREE, M, A = GT.OE, GT.KD, GT.FREE, GT.M, GT.A



"""
	The Network class sets up and contains the gene network. Takes dictionaries
	geneData	and commands as input which are created by the lgnReader.
"""
class Network:
	def __init__(self, geneData, commands):
		self.geneData = geneData
		self.commands = commands

		self.N = self.geneData['N'] # Number of genes in network.
		self.genes = self.geneData['genes'] # List of genes and their input.
		self.ops = self.geneData['ops'] # List of operators.
		self.specInput = self.geneData['specInput'] # A specified configuration as input.
		self.constraints = self.geneData['constraints'] # List of states that are expected to be attractors.
		self.names = self.geneData['names'] # List of gene names.
		self.inputs = self.geneData['inputs'] # List of gene inputs

		self.indices = {self.names[i]:i for i in range(len(self.names))}	# Dictionary of names corresponding do indices

		self.geneStatus = []	#If the genes are OE/KD/FREE. Is filled when adjecency matrix is formed.
		self.adj = self.constructAdjMatrix()	# Adjecency matrix of network
		self.network = self.constructConAdjMatrix()	# Condensed adjecency matrix

		# If no operators are given, then a randomised priority list is used instead.
		## Not implemented yet though.
		if(len(self.ops) == 0):
			self.priorityList = self.randomisePriorityList()


		self.currentRun = RunData(self)

		# The gene trees that make up the logical structure of the network
		# interactions. These are changed when various cofigurations are tried.
		self.geneForest = [0 for x in range(self.N)]
		self.keys = [0 for x in range(self.N)]

		# Flags if one or some of the trees' configuration is specifiec.
		self.fixedFlag = self.commands['fixedFlag']

		# Filename used to save results.
		self.filename = self.commands['filename']

		# Path where result file is saved.
		self.PATH = self.commands['PATH']
		
		# Sets up the initial gene forest.
		self.formNetwork()

		# Determines which trees are trivial and which are not.
		ntr = 0
		while(ntr < self.N and self.geneForest[ntr].treeOrder < 2):
			ntr += 1
		self.nontriv = ntr # The index of the first nontrivial (tree order > 1) tree in the forest.

		# Storage lists for the configurations with the best consistency score.
		self.bestTau = []	# Stores the energy.
		self.bestConfigurations = [] # Stores the configuration strings.
		self.bestSpaces = [] # Stores the expression spaces.
		self.bestConsistencyScore = -math.inf # The so far best achieved consistency score.

		return


	"""
		Builds the network (or rather the trees list (forest)) by, for each
		gene, building a tree containing its input genes. This method is only
		called the first time this is done.
	"""
	def formNetwork(self):
		print("\nForming network...")

		tmp = self.geneForest.copy()

		# Creates a tree for each gene's input.
		for i in range(self.N):
			# If the configuration for this tree is given.
			if(self.names[i] in self.specInput):
				tmp[i] = GT.GeneTree(self.ops, self.geneStatus, stringList=[self.names[i], self.specInput[self.names[i]]], names=self.names, fixed=self.fixedFlag)
			# If the configuration for this tree is not given.
			else:
				tmp[i] = GT.GeneTree(self.ops, self.geneStatus, target=i, genes=self.network[i])

		# Sorts the forest for more efficient iteration.
		self.geneForest = sorted(tmp, key=lambda a: a.treeOrder)


		format = ("(%" + str(len(str(self.geneForest[-1].getNumberOfConfigs())))
					+ "d) %s\n")

		# Prints out the number of configurations for each tree as well as
		# its starting configuration.
		for i in range(len(self.geneForest)):
			self.keys[i] = self.geneForest[i].geneKey
			print(format % (self.geneForest[i].getNumberOfConfigs(),
						   self.geneForest[i].toString(self.names)))


	"""
		Builds the network (or rather the trees list (forest)) by, for each
		gene, building a tree containing its input genes. This resets all
		changes previously made to the configurations of the network.
	"""
	def reformNetwork(self):
		print("\nReforming network...")
		for i in range(self.nontriv, self.N):
			tgt = self.geneForest[i].target
			self.geneForest[i] = GT.GeneTree(tgt, self.network[tgt], self.ops)

	""" Calculates the total number of configurations in the forest. """
	def getNumberOfConfigs(self):
		conf = 1
		for tree in self.geneForest:
			conf *= tree.getNumberOfConfigs()
		return conf


	"""
		Iterates through all the configurations of the forest, calculating
		the energy for each and saves the configurations with the best
		consistency score.
	"""
	def runExhaustive(self):
		# Open a file to write the results.
		out = open(self.PATH+self.filename, 'w')

		# Writes out a preamble.
		out.write("### In file \"" + self.PATH + self.filename + "\"\n")
		out.write("###\n")
		out.write("### This file was generated by CELLoGeNe, " +
					"written and developed by Emil Andersson and Mattias Sjö.\n")
		out.write("### File created " + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + "\n")
		out.write("\n")

		out.write('#Using genes: ' + str(self.names))
		out.write('\n#Network connections: ' + str(self.inputs))
		out.write('\n#Using operators: ' + str(self.ops))
		out.write('\n#Constraints: ' + str([('{:0>'+str(self.N)+'}').format(bin(c)[2:]) for c in self.constraints]))
		out.write('\n#There is : ' + str(self.getNumberOfConfigs()) + ' configurations')
		out.write('\n#Exhaustive search through configuration space.')
		out.write('\n')



		count = 1 # Keeps count on how many configurations that have been iterated over.
		self.currentRun.initialise()

		# If only 1 operator is used, all thre trees are trivial, or all the
		# trees configurations are specified, then there is only 1 configuration.
		if(self.geneForest[0].p < 2 or self.nontriv >= self.N or
				 (self.fixedFlag and len(self.specInput) == self.N)):
			print("\nNetwork has only 1 configuration.")
			out.write('\n#Network has only 1 configuration')
			out.write('\n')
			self.bestTau = [[self.currentRun.space.points[i].value for i in range(1<<self.N)]]
			self.bestConfigurations = [self.recordConfig()]
			self.bestConsistencyScore = self.consistencyScore()
			

			self.currentRun.report(display=1<<self.N ,out=out)
			
			if(self.bestConsistencyScore < 0):
				self.currentRun.space.findInclusiveBasins()
				self.currentRun.space.findExclusiveBasins()
			self.bestSpaces = [copy.deepcopy(self.currentRun.space)]

			
			print('Best consistency score: %.3f\tNumber of best configurations: %d'
				% (self.bestConsistencyScore, len(self.bestConfigurations)))
			out.write('\n\n')
			out.write('#Best consistency score: %.3f\n#Number of best configurations: %d\n'
				%(self.bestConsistencyScore, len(self.bestConfigurations)))
			out.write('\n')
			if(len(np.where(np.array(self.geneStatus) == FREE)[0]) <= 20): 
				print(self.bestTau[0])
				out.write('\n#conf nr\tconfiguration\ttau\tInclusive\tExclusive\n')
				out.write('%d\t%s\t%s\t%s\t%s\n'
							 % (1, str(self.bestConfigurations[0]),
								str(self.bestTau[0]),
								str({str([('{:0>'+str(self.N)+'}').format(bin(s)[2:])
									for s in a['states']]): len(a['inclusive'])
									for a in self.bestSpaces[0].attractors}),
								str({str([('{:0>'+str(self.N)+'}').format(bin(s)[2:])
									for s in a['states']]): len(a['exclusive'])
									for a in self.bestSpaces[0].attractors})))
			else:
				print('OR HERE?')
				out.write('\n#conf nr\tconfiguration\tInclusive\tExclusive\n')
				out.write('%d\t%s\t%s\t%s\n'
							 % (1, str(self.bestConfigurations[0]),
								str({str([('{:0>'+str(self.N)+'}').format(bin(s)[2:])
									for s in a['states']]): len(a['inclusive'])
									for a in self.bestSpaces[0].attractors}),
								str({str([('{:0>'+str(self.N)+'}').format(bin(s)[2:])
									for s in a['states']]): len(a['exclusive'])
									for a in self.bestSpaces[0].attractors})))

			print('Best consistency score: %.3f\tNumber of best configurations: %d'
				% (self.bestConsistencyScore, len(self.bestConfigurations)))


			return


		# If more than 1 configuration.

		vals = [self.geneForest[i].newEnergy for i in range(self.N)]

		total = self.getNumberOfConfigs()
		format = "%" + str(len(str(total))) + "d configurations done (%3d %%)\n"
		percent = total / 100
		p = percent
		out.write('\n#There is ' + str(total) + ' configurations')

		print("\nNetwork predicted to have " + str(total) + " configurations.")
		print("\nTrying all configurations of the network...")

		# Calculate and save configuration, energy and consistency score
		# for the first configuration
		self.bestConfigurations.append(self.recordConfig())
		self.bestTau.append([self.currentRun.space.points[i].value for i in range(1<<self.N)])
		self.bestConsistencyScore = self.consistencyScore()
		self.bestSpaces.append(copy.deepcopy(self.currentRun.space))


		start = time.time()

		done = [False for i in range(self.N)]
		i = self.nontriv

		# The master loop which iterates through all configurations.
		while(i < self.N):
			# Move to the next configuration
			done[i] = self.geneForest[i].nextConfig()
			self.geneForest[i].updateEnergy()
			vals[i] = self.geneForest[i].newEnergy

			# If all configurations of tree i has been iterated over,
			# move on to next tree.
			if(done[i]):
				i += 1
			else:
				self.currentRun.space.update_vals(vals)
				self.currentRun.space.check()

				i = self.nontriv

				# Compile the energy landscape tau for this configuration
				tau = [self.currentRun.space.points[x].value for x in range(1<<self.N)]


				# Check consistency score. If same as previously highest,
				# append to lists of bests, if higher than earlier, reseset
				# lists, if lower, move on.
				C = self.consistencyScore()
				if(C >= self.bestConsistencyScore):
					if(C > self.bestConsistencyScore):
						self.bestConfigurations = []
						self.bestTau = []
						self.bestSpaces = []
						self.bestConsistencyScore = C
					self.bestConfigurations.append(self.recordConfig())
					self.bestTau.append(tau)
					self.bestSpaces.append(copy.deepcopy(self.currentRun.space))



				count += 1
				# Prints out progress to the terminal
				if(count >= p):
					while True:
						p += percent
						if(p > count):
							break
					self.time_elapsed(start)
					print(format % (count, int(p / percent) - 1))

		print("Tried " + str(count) + " configurations.")
		# Reports (both to terminal and file) how many time each state has been an attractor
		self.currentRun.report(display=1<<self.N ,out=out)

		print('Best consistency score: %.3f\tNumber of best configurations: %d'
			% (self.bestConsistencyScore, len(self.bestConfigurations)))
		out.write('\n\n')
		out.write('#Best consistency score: %.3f\n#Number of best configurations: %d\n'
			%(self.bestConsistencyScore, len(self.bestConfigurations)))
		out.write('\n')
		if(len(np.where(np.array(self.geneStatus) == FREE)[0]) <= 20): 
			out.write('\n#conf nr\tconfiguration\ttau\tInclusive\tExclusive\n')
		else:
			out.write('\n#conf nr\tconfiguration\tInclusive\tExclusive\n')

		# Saves the best configuration, their energy and basins of attractions
		# in file.
		for i in range(len(self.bestConfigurations)):
			print('%s'
					 % (str(self.bestConfigurations[i])))
#					   str(self.bestTau[i])))
		for i in range(len(self.bestConfigurations)):
			self.bestSpaces[i].findInclusiveBasins()
			self.bestSpaces[i].findExclusiveBasins()
			if(len(np.where(np.array(self.geneStatus) == FREE)[0]) <= 20):
				out.write('%d\t%s\t%s\t%s\t%s\n'
						 % (i+1, str(self.bestConfigurations[i]),
							str(self.bestTau[i]),
							str({str([('{:0>'+str(self.N)+'}').format(bin(s)[2:])
								for s in a['states']]): len(a['inclusive'])
								for a in self.bestSpaces[i].attractors}),
							str({str([('{:0>'+str(self.N)+'}').format(bin(s)[2:])
								for s in a['states']]): len(a['exclusive'])
								for a in self.bestSpaces[i].attractors})))
			else:
				out.write('%d\t%s\t%s\t%s\n'
						 % (i+1, str(self.bestConfigurations[i]),
							str({str([('{:0>'+str(self.N)+'}').format(bin(s)[2:])
								for s in a['states']]): len(a['inclusive'])
								for a in self.bestSpaces[i].attractors}),
							str({str([('{:0>'+str(self.N)+'}').format(bin(s)[2:])
								for s in a['states']]): len(a['exclusive'])
								for a in self.bestSpaces[i].attractors})))

			# Print to terminal for check.
			print(self.bestConfigurations[i])
		print('Best consistency score: %.3f\tNumber of best configurations: %d'
			% (self.bestConsistencyScore, len(self.bestConfigurations)))

		out.close()
		return


	"""
		Iterates through a set number or fractions of the configurations of the forest, calculating
		the energy for each and saves the configurations with the best
		consistency score.
	"""
	def runMonteCarlo(self, frac=0, tot=0):
		self.currentRun.initialise()

		# Open a file to write the results.
		out = open(self.PATH+self.filename, 'w')

		# Writes out a preamble.
		out.write("### In file \"" + self.PATH + self.filename + "\"\n")
		out.write("###\n")
		out.write("### This file was generated by CELLoGeNe, " +
					"written and developed by Emil Andersson and Mattias Sjö.\n")
		out.write("### File created " + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + "\n")
		out.write("\n")

		out.write('#Using genes: ' + str(self.names))
		out.write('\n#Network connections: ' + str(self.inputs))
		out.write('\n#Using operators: ' + str(self.ops))
		out.write('\n#Constraints: ' + str([('{:0>'+str(self.N)+'}').format(bin(c)[2:]) for c in self.constraints]))
		out.write('\n#There is : ' + str(self.getNumberOfConfigs()) + ' configurations')

		# If only 1 operator is used, all thre trees are trivial, or all the
		# trees configurations are specified, then there is only 1 configuration.
		if(self.geneForest[0].p < 2 or self.nontriv >= self.N or
				 (self.fixedFlag and len(self.specInput) == self.N)):
			print("\nNetwork has only 1 configuration.")
			out.write('\n#Network has only 1 configuration')
			out.write('\n')
			self.bestTau = [[self.currentRun.space.points[i].value for i in range(1<<self.N)]]
			self.bestConfigurations = [self.recordConfig()]
			self.bestConsistencyScore = self.consistencyScore()

			self.currentRun.report(display=1<<self.N ,out=out)

			# Print to terminal for check.
			print('best consistency score', self.bestConsistencyScore)
			for i in range(len(self.bestConfigurations)):
				print(self.bestConfigurations[i])
				print(self.bestTau[i])

			# Find the basins of the energy landscape.
			self.currentRun.space.findInclusiveBasins()
			self.currentRun.space.findExclusiveBasins()

			# Print to terminal for check.
			self.printBasinData()

			return

		# When there are more than 1 configuration.

		vals = [self.geneForest[i].newEnergy for i in range(self.N)]

		total = self.getNumberOfConfigs()
		print('Network predicted to have {:d} configurations'.format(total))
		out.write('\n#Stochastic search through configuration space.')


		if(frac > 0):
			total *= frac
		else:
			if(type(tot) is int):
				total = tot
		format = '{:' + str(len(str(total))) + 'd} configs done ({:3d} %)\n'

		percent = total / 100
		p = percent

		print('\nTrying approx {:d} random configurations...\n'.format(math.ceil(total)))
		
		# Printout format
		sfrac = str(frac*100)
		dec = len(sfrac[sfrac.find('.')+1:].rstrip('0'))
		if(frac == 0):
			out.write(('\n#Trying {:d} of the configurations.').format(total))
		else:
			out.write(('\n#Trying {:' + ('.'+str(dec)+'f' if dec < 4 else '.1e') + '} % ({:d}) of the configurations.').format(frac*100 if frac != 0 else 100, math.ceil(total)))
		out.write('\n')

		start = time.time()

		count = 0

		while(count < total):
			for i in range(self.nontriv, self.N):
				self.geneForest[i].shuffle()
				self.geneForest[i].updateEnergy()
				vals[i] = self.geneForest[i].newEnergy

			self.currentRun.space.update_vals(vals)
			self.currentRun.space.check()
			tau = [self.currentRun.space.points[x].value for x in range(1<<self.N)]

			# Check consistency score. If same as previously highest,
			# append to lists of bests, if higher than earlier, reseset
			# lists, if lower, move on.
			C = self.consistencyScore()
			if(C > self.bestConsistencyScore):
				if(C > self.bestConsistencyScore):
					self.bestConfigurations = []
					self.bestTau = []
					self.bestSpaces = []
					self.bestConsistencyScore = C
				self.bestConfigurations.append(self.recordConfig())
				self.bestTau.append(tau)
				self.bestSpaces.append(copy.deepcopy(self.currentRun.space))

			count += 1
			if(count >= p):
				p = count / percent

				self.time_elapsed(start)
				print(format.format(count, int(p)))

				p = count + percent
		print('Tried {:d} configurations\n'.format(count))
		self.currentRun.report(display=1<<self.N ,out=out)

		print('Best consistency score: %.3f\tNumber of best configurations: %d'
			% (self.bestConsistencyScore, len(self.bestConfigurations)))
		out.write('\n\n')
		out.write('#Best consistency score: %.3f\n#Number of best configurations: %d\n'
			%(self.bestConsistencyScore, len(self.bestConfigurations)))
		out.write('\n')
		if(len(np.where(np.array(self.geneStatus) == FREE)[0]) <= 20):
			out.write('\n#conf nr\tconfiguration\ttau\tInclusive\tExclusive\n')
		else:
			out.write('\n#conf nr\tconfiguration\tInclusive\tExclusive\n')

		# Saves the best configuration, their energy and basins of attractions
		# in file.
		for i in range(len(self.bestConfigurations)):
			print('%s'
					 % (str(self.bestConfigurations[i])))
		for i in range(len(self.bestConfigurations)):
			self.bestSpaces[i].findInclusiveBasins()
			self.bestSpaces[i].findExclusiveBasins()
			if(len(np.where(np.array(self.geneStatus) == FREE)[0]) <= 20):
				out.write('%d\t%s\t%s\t%s\t%s\n'
						 % (i+1, str(self.bestConfigurations[i]),
							str(self.bestTau[i]),
							str({str([('{:0>'+str(self.N)+'}').format(bin(s)[2:])
								for s in a['states']]): len(a['inclusive'])
								for a in self.bestSpaces[i].attractors}),
							str({str([('{:0>'+str(self.N)+'}').format(bin(s)[2:])
								for s in a['states']]): len(a['exclusive'])
								for a in self.bestSpaces[i].attractors})))
			else:
				out.write('%d\t%s\t%s\t%s\n'
						 % (i+1, str(self.bestConfigurations[i]),
							str({str([('{:0>'+str(self.N)+'}').format(bin(s)[2:])
								for s in a['states']]): len(a['inclusive'])
								for a in self.bestSpaces[i].attractors}),
							str({str([('{:0>'+str(self.N)+'}').format(bin(s)[2:])
								for s in a['states']]): len(a['exclusive'])
								for a in self.bestSpaces[i].attractors})))

			# Print to terminal for check.
			print(self.bestConfigurations[i])
		print('Best consistency score: %.3f\tNumber of best configurations: %d'
			% (self.bestConsistencyScore, len(self.bestConfigurations)))

		out.close()





	""" Records the current configuration"""
	def recordConfig(self):
		return [self.geneForest[x].toString(self.names) for x in range(self.N)]

	""" Loads a configuration stored using recordConfiguration()"""
	def loadConfig(self, config):
		self.geneForest = [GT.GeneTree(self.ops, self.geneStatus, stringList=config[i].split(' <- '),
								 names=self.names) for i in range(self.N)]
		self.geneForest = sorted(self.geneForest, key=lambda a: a.treeOrder)



	
	"""
		Checks the consistency score, and returns it, for either the current
		or a given expression space.
	"""
	def consistencyScore(self, space=None):
		epsilon = 1
		B, R, W = self.checkConstraints(space)
		C = B if B != 0 else -W/(R + epsilon)

		return C
	
	"""
		The consistency score is based on the sizes of the inclusive and 
		exclusive (soft and strict) basins and the degeneracy, if the current 
		attractor is a constraint, otherwise it is based on gradients.
		The value of the score does not say much if it is more or less 
		biologically plausible. 
		
		The most (as of now only(?)) important feature of the consistency score
		is that it is only positive if ALL constraints are fullfilled,
		otherwise it is negative.
	"""
	def checkConstraints(self, space=None):
		if(space is None):
			space = self.currentRun.space
		n_c = len(self.constraints)
		correct = []
		B, R, W = 0, 0, 0
		alpha = 0.75
		epsilon = 1
		for c in self.constraints:
			correct.append(True if (any(c in a['states'] for a in space.attractors)) else False)
		if(all(correct)):
			space.findInclusiveBasins()
			space.findExclusiveBasins()
			for i in range(n_c):
				for a in space.attractors:
					if(self.constraints[i] in a['states']):
						B += (alpha * len(a['exclusive']) + (1 - alpha) * len(a['inclusive'])) / (2 **(a['degeneracy'] - 1)) # Alternative form
		else:
			for i in range(n_c):
				if(correct[i]):
					for g in space.points[self.constraints[i]].grad:
						if(g > 0):
							print('i', i, 'g', g)
						R += -g
				else:
					pos = 0
					neg = 0
					for g in space.points[self.constraints[i]].grad:
						if(g > 0):
							pos += g
						else:
							neg += -g
					W += neg / (pos + epsilon)
		return B, R, W





	"""
		Constructs the adjecency matrix from the gene input data and sets the
		gene status.
	"""
	def constructAdjMatrix(self):
		matrix = []
		for gene in self.inputs:
			if(gene[0] == 'X'):
				self.geneStatus.append(OE)
				gene = '.' * self.N
			elif(gene[0] == 'K'):
				self.geneStatus.append(KD)
				gene = '.' * self.N
			elif(gene[0] == 'M'):
				self.geneStatus.append(M)
				gene = '.' * self.N
			elif(gene[0] == 'A'):
				self.geneStatus.append(A)
				gene = '.' * self.N
			elif((gene.find('[') != 0 or gene.find(']') != self.N+1)
				or len(gene) != self.N+2):
				print("File on wrong format")
			else:
				gene = gene[1:self.N+1]
				self.geneStatus.append(FREE)
			row = []
			for a in gene:
				if(a == '+'):
					row.append(ACT)
				elif(a == '-'):
					row.append(REP)
				elif(a == '.'):
					row.append(NON)
				elif(a == '_'):
					row.append(random.choice([REP, NON, ACT]))
			matrix.append(row)
		return matrix

	"""Constructs the condesed adjency matrix from the adjency matrix. """
	def constructConAdjMatrix(self):
		cond = []
		for row in self.adj:
			inputs = []
			# Since 0 is neither - nor + a special trick with signed zeroes
			# is used for the first element.
			if(row[0] == ACT):
				inputs.append(0.0)
			elif(row[0] == REP):
				inputs.append(-0.0)
			for i in range(1, len(row)):
				if(row[i] == ACT):
					inputs.append(i)
				elif(row[i] == REP):
					inputs.append(-i)
			cond.append(inputs)
		return cond


	""" Randomises the priority list. """
	def randomisePriorityList(self):
		priority = [x for x in range(self.N)]
		random.shuffle(priority)
		return priority

	""" Sets the priority list (not used)"""
	def setPriorityList(self, order):
		self.priorityList = [self.indices[x] for x in order]
		return



	"""
		Converts a configuration list (such as generated by
		recordConfiguration) inta a dictionary. (Not used but usefull)
	"""
	def listToDic(self, conf):
		dic = {}
		for c in conf:
			tmp = c.split(' <- ')
			dic.update({tmp[0].strip(): tmp[1].strip()})
		return dic


	""" Prints basin data to the terminal. """
	def printBasinData(self, space=None):
		if(space is None):
			space = self.currentRun.space
		print('\nResulting basins\n')

		for x in space.attractors:
			print([('{:0>'+ str(self.N) + '}').format(str(bin(y)[2:])) for y in x['states']])

		print()
		for a in space.attractors:
			print('\nAttractor:')
			print([('{:0>'+ str(self.N) + '}').format(str(bin(x)[2:])) for x in a['states']],'\n')
			print('Degeneracy: {}'.format(a['degeneracy']))
			print('Inclusive Basin:')
			print([('{:0>'+ str(self.N) + '}').format(str(bin(s)[2:])) for s in a['inclusive']])
			print('Bin size', len(a['inclusive']))
			print()
			print('Exclusive Basin:')
			print([('{:0>'+ str(self.N) + '}').format(str(bin(s)[2:])) for s in a['exclusive']])
			print('Bin size', len(a['exclusive']))

			print()


	""" Prints the time elapsed since the given start time. """
	def time_elapsed(self, start_time):
		run_time = time.time() - start_time
		hours = run_time // 3600
		minutes = (run_time % 3600) // 60
		seconds = run_time - hours * 3600 - minutes * 60
		print("--- Time elapsed: %d h %d min and %.3f s --- Current time: %s ---"
			% (int(hours), int(minutes), seconds, time.strftime("%H:%M:%S", time.localtime())))


