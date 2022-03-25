""" Operators """

"""
	The operators are most easily accessed from other modules by importing
	this module and calling the dictionary 'operators'.
"""

class Operators:

	""" Values used in three-state logic. A more natural choice would be
		REP = -1, NON = 0 and ACT = +1, but the shifted values simplify
		the implementation.
	"""
	REP, NON, ACT = 0, 1, 2


	""" The rather trivial && operator.
		If either input is NON,
		or if the inputs are opposite, returns NON.
	"""
	AND = 	[
			[REP, NON, NON],
			[NON, NON, NON],
			[NON, NON, ACT]
			]
	""" OBAND and UBAND are not monotonous, meaning that they can turn activation
		into repression and vice versa. Therefore, their use it not recommended.
    """

	"""	The over-barred && operator.
		If either input is NON, returns NON.
		Maps opposite inputs to ACT.
	"""
	OBAND = [
			[REP, NON, ACT],
			[NON, NON, NON],
			[ACT, NON, ACT]
			]

	"""	The under-barred && operator.
		If either input is NON, returns NON.
		Maps opposite inputs to REP.
	"""
	UBAND = [
			[REP, NON, REP],
			[NON, NON, NON],
			[REP, NON, ACT]
			]

	"""	The basic || operator.
		In a sense the most natural operator, but unfortunately not associative.
		Therefore, it should only be used in stochastic methods, as the
		exhaustive methods assume associativity.
	"""
	OR = 	[
			[REP, REP, NON],
			[REP, NON, ACT],
			[NON, ACT, ACT]
			]
	"""	The over-barred || operator.
		If either input is not NON, it returns that.
		Maps opposite inputs to ACT.
	"""
	OBOR = 	[
			[REP, REP, ACT],
			[REP, NON, ACT],
			[ACT, ACT, ACT]
			]

	"""	The under-barred || operator.
		If either input is not NON, it returns that.
		Maps opposite inputs to REP.
	"""
	UBOR = 	[
			[REP, REP, REP],
			[REP, NON, ACT],
			[REP, ACT, ACT]
			]

	"""	The minimum operator.
		Returns its smallest output.
	"""
	MIN = 	[
			[REP, REP, REP],
			[REP, NON, NON],
			[REP, NON, ACT]
			]

	"""	The maximum operator.
		Returns its largest output.
	"""
	MAX = 	[
			[REP, NON, ACT],
			[NON, NON, ACT],
			[ACT, ACT, ACT]
			]

	"""
		The operators are most easily accessed from other modules by importing
		this module and calling this dictionary.
	"""
	operators = {	'AND': AND, 'OBAND': OBAND, 'UBAND': UBAND,
					'OBOR': OBOR, 'UBOR': UBOR, 'MIN': MIN, 'MAX': MAX}
	operators_names = ['AND', 'OBAND', 'UBAND', 'OBOR', 'UBOR', 'MIN', 'MAX']

