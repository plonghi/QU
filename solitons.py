class Dash:
	"""
	A path on Sigma, characterized by an ordered sequence of streets, 
	with a choice of orientation for each.
	The orientation is +1 if the dash runs from the first .endpoint
	to the second .endpoint of the street. It is -1 otherwise.
	The attribute .path collects this information in the form
	[[p_1, +1], [p_2, -1], ...]

	A dash can be extended, from either end. The extension is specified 
	by just adding a street. The orientation will be determined 
	automatically by the choice of joint or branch point at which 
	the extension is performed.
	"""
	def __init__(self, label=None):
		self.label = label
		self.path = []
	
	def extend_dash_along_street(self, end_pt=None, street=None):
		if len(self.path) == 0:
			orientation = set_orientation_from_starting_point(street, end_pt)
			self.path = [street, orientation]
		elif end_pt == 'last':
			terminal_street = self.path[-1]



class SolitonPath:
	"""
	The attribute .streets_set just tracks the homology of a soliton.
	But the attribute .path contains full information
	about its path as it winds on the Seiberg-Witten curve.
	"""
	def __init__(self, label=None):
		self.label = label
		self.streets_set = []
		self.path = []
		self.starting_point = None
		self.ending_point = None
		

class ClosedSoliton:
	"""
	The closed path obtained by joining two open paths
	of opposite types and supported on the same street.
	"""
	def __init__(self, label=None, soliton_a=None, soliton_b=None):
		self.label = label
		self.streets_set = []
		self.path = []
		self.homology = None
		pass


def set_orientation_from_starting_point(street, starting_pt):
	if starting_pt == street.endpoints[0]:
		return +1
	elif starting_pt == street.endpoints[1]:
		return -1
	else:
		raise Exception(
			'The point {} is not an endpoint of street {}'
			.format(starting_pt.label, street.label)
		)





