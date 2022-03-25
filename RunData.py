from ExpressionSpace import ExpressionSpace

"""
	Contains the necessary data needed to perform an attractor-finding run of
	the networks.
"""
class RunData:
	def __init__(self, network_self):
		# The self-variables of network.
		self._NS = network_self

	""" Initialises the ExpressionSpace. """
	def initialise(self):
		print('\nInitialising...\n')

		self.space = ExpressionSpace(self._NS.keys)
		
		for i in range(self._NS.N):
			if(self._NS.geneForest[i].treeOrder < 1):
				continue
			self._NS.geneForest[i].updateEnergy()
			self.space.update(self._NS.geneForest[i].oldEnergy,
							 self._NS.geneForest[i].newEnergy,
							 self._NS.geneForest[i].geneKey)
		self.space.check()

	""" Report on the prevalence of attractors. """
	def report(self, display=51, out=None):
		maxCount = self.space.rankAttractors(display, out=out)
		print("Maximum: " + str(maxCount))

