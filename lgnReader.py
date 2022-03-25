""" This module is used to read .lgn-files and can be used to extract the
geneData and commands"""


""" Imports """
from Operators import Operators
import os


add = 3	#numberOfAdditionalRows: operators, specInput, constraints

""" The reader class. Opens and trims the .lgn-files """
class Reader:
	def __init__(self, file_name):
		self.file_name = file_name
		with open(self.file_name, "r") as f:
			self.lines = f.readlines()
		self.lines = self.lineReader(self.lines)
		self.N = int(self.lines[0].split(':')[1].strip())	# The number of genes
		self.geneData = {}
		self.setGeneData()
		self.commands = {}
		self.setCommands()
		return




	""" Reads each line but skips empty lines an comments"""
	def lineReader(self, strings):

		""" Sub method that is called from findLongComment which removes
		long comments recursively. """
		def removeLongComment(strings, j):
			start = strings[j].find('"""')
			tmp = strings[j][start+3:]
			save = strings[j][:start]
			stop = tmp.find('"""')
			if(stop != -1):
				strings[j] = save + tmp.replace(tmp[:stop+3], '')
				if(strings[j].count('"""') > 0):
							strings, j = removeLongComment(strings, j)
				else:
					return strings, j
				return strings, j
			else:
				strings[j] = strings[j].replace(strings[j][start:], '')
				for i in range(j + 1,len(strings)):
					if(strings[i].count('"""') == 0):
						strings[i] = ''
					else:
						strings[i] = strings[i][strings[i].find('"""')+3:]
						if(strings[i].count('"""') > 0):
							sub, j = removeLongComment(strings[i:], i)
							strings[:i] += sub
						return strings, i + j
				return strings, len(strings)

		""" Sub method used to find where a long comment starts. """
		def findLongComment(strings):
			skip = 0
			if(strings[0].count('"""') == 0):
				return (strings, skip)
			else:
				strings, skip = removeLongComment(strings, 0)
				return strings, skip


		count = 0
		for i in range(len(strings)):
			i = count
			""" Removes #-comments """
			if(strings[i].find('#') != -1):
				strings[i] = strings[i].replace(strings[i][strings[i].find('#'):], '')
			""" Removes long comments """
			save = strings[:i]
			sub, k = findLongComment(strings[i:])
			strings = save + sub
			count += 1 + k
			""" Removes leading and trailing withespaces """
			strings[i] = strings[i].strip()
			if(count == len(strings)):
				break
		""" Removes empty lines and removes leading and trailing withespaces
		from non-empty lines"""
		strings = [x.strip() for x in strings if x.strip()]
		return strings


	""" Returns a dictionary were each entry is a row containing input for one gene. """
	def setGeneData(self):
		self.geneData['N'] = self.N
		self.geneData['genes'] = self.lines[1:self.N+1]
		self.geneData['ops'] = self.extractOperators(self.lines[self.N+1].split(':')[1].strip())
		self.geneData['specInput'] = self.extractDic(self.lines[self.N+2].split('=')[1].strip())
		self.geneData['constraints'] = self.extractConstraints(self.lines[self.N+3].split(':')[1].strip())
		self.geneData['names'] = [x.split(':')[0].strip() for x in self.geneData['genes']]
		self.geneData['mapName'] = {self.geneData['names'][i]: i for i in range(self.N)}
		self.geneData['inputs'] = self.encodeInput([x.split(':')[1].strip() for x in self.geneData['genes']])


		return

	""" Returns a dictionary were each entry contains information about a command"""
	def setCommands(self):
		self.commands['fixedFlag'] = True if self.lines[self.N+1+add].split('=')[1].strip() == 'True' else False
		self.commands['filename'] = self.lines[self.N+2+add].split(':')[1].strip()
		tmp = self.lines[self.N+3+add].split(':')[1].strip()
		self.commands['PATH'] = os.getenv('HOME') + tmp[1:] if tmp[0] == '~' else tmp
		self.extractRunType(self.lines[self.N+4+add])

		return


	""" Encodes the readable input format to the input format used by the Network class. """
	def encodeInput(self, inList):
		encoded = []
		for i in range(self.N):
			K, X, M, A = False, False, False, False
			if(inList[i][0] == 'K'):
				K = True
				inList[i] = inList[i][1:]
			elif(inList[i][0] == 'X'):
				X = True
				inList[i] = inList[i][1:]
			elif(inList[i][0] == 'M'):
				M = True
				inList[i] = inList[i][1:]
			elif(inList[i][0] == 'A'):
				A = True
				inList[i] = inList[i][1:]
			
			split = [s.strip() for s in inList[i][1:-1].split(',')]
			in_string = ''
			for g in self.geneData['names']:
				if(not any(g == s[1:] for s in split)):
					in_string += '.'
				elif(any('+{:s}'.format(g) == s for s in split)):
					in_string += '+'
				elif(any('-{:s}'.format(g) == s for s in split)):
					in_string += '-'
				elif(any('?{:s}'.format(g) == s for s in split)):
					in_string += '_'
			in_string = '[' + in_string + ']'
			if(K):
				in_string = 'K' + in_string
			elif(X):
				in_string = 'X' + in_string
			elif(M):
				in_string = 'M' + in_string
			elif(A):
				in_string = 'A' + in_string
			encoded.append(in_string)
		
		
		return encoded

	""" Extracts the input string of operators into a list of operator names. """
	def extractOperators(self, opString):
		op = []
		if(opString.find('[') != 0 or opString.find(']') != len(opString)-1):
				print("File on wrong format (1)")

		else:
			opString = opString[1:-1]
			oper = opString.strip().split(',')
		allowed = Operators.operators_names # These are: ['AND', 'OBAND', 'UBAND', 'OBOR', 'UBOR', 'MIN', 'MAX']
		for o in oper:
			o = o.strip()
			if(o not in allowed):
				print("%s is not an allowed operator." % o)
			else:
				if(o not in op):
					op.append(o)
				else:
					print("Duplicates of operator %s ." % o)
		print('Using operators:', op)
		return op

	""" Extracts a dictionary from a string-representation of a dictionary. """
	def extractDic(self, dicToBe):
		if(dicToBe.find('{') != 0 or dicToBe.find('}') != len(dicToBe)-1):
			print("File on wrong format (2)")
		dicToBe = dicToBe[1:-1]
		if(dicToBe == ''):
			return {}
		dicToBe = dicToBe.strip().split(',')
		dicToBe = [dicToBe[i].strip().split(':') for i in range(len(dicToBe))]
		dic = {dicToBe[i][0]: dicToBe[i][1].strip() for i in range(len(dicToBe))}
		return dic


	""" Extracts the constraints as a list from a given string. """
	def extractConstraints(self, constr):
		if(constr[0] != '[' or constr[-1] != ']'):
				print("File on wrong format (3)")
		else:
			constr = constr[1:-1]
			constr = constr.strip().split(',')
			constr = [c.strip() for c in constr]
			constraints = [0 for x in range(len(constr))]
			for i in range(len(constr)):
				m = 1
				for c in reversed(constr[i]):
					if(c == '1'):
						constraints[i] |= m
					m <<= 1
		return constraints

	"""
		Converts a configuration list (such as generated by
		recordConfiguration) inta a dictionary.
	"""
	def listToDic(self, conf):
		dic = {}
		for c in conf:
			tmp = c.split(' <- ')
			dic.update({tmp[0].strip(): tmp[1].strip()})
		return dic


	def extractRunType(self, string):
		split = string.split(':')
		runType = split[1].strip()
		if(runType == 'E'):
			self.commands['runType'] = 'E'
			return
		elif(runType == 'M'):
			self.commands['runType'] = 'M'
			if(len(split) > 2):
				for e in split[2:]:
					keyword = e.split('=')[0].strip()
					arg = e.split('=')[1].strip()
					if(keyword == 'frac'):
						self.commands['frac'] = float(arg)
					elif(keyword == 'tot'):
						if(arg == 'tot'):
							self.commands['tot'] = arg
						else:
							self.commands['tot'] = int(float(arg))
					else:
						print('Unvalid keyword.')
			else:
				print('Additional keywors argument is necessary for Monte Carlo')
		else:
			print('Unvalid run type.')
		return
		
