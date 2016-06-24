from mcsn import MCSN

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

	One may impose restrictions on how a dash can grow: forward_only, 
	backward_only, or None (if can grow both ways).
	"""
	def __init__(self, label=None, growth_restriction=None,):
		self.label = label
		self.path = []
		if growth_restriction is not None:
			self.growth_restriction = growth_restriction
		else:
			self.growth_restriction = None
	
	def extend_dash_along_street(self, street=None, end_pt=None, slot=None):
		if len(self.path) == 0:
			if end_pt is None or slot is None:
				raise Exception(
					'For the first street of a dash, '
					'a starting point must be given.'
				)
			else:
				new_street_orientation = set_orientation_from_starting_point(
					street, end_pt, slot
				)
				self.path = [[street, new_street_orientation]]
		
		elif end_pt == 'last' and (
				self.growth_restriction == 'forward_only' or 
				self.growth_restriction is None
			):
			terminal_street, terminal_orientation = self.path[-1]
			if terminal_orientation == +1:
				intermediate_pt = terminal_street.final_point().end_point
			elif terminal_orientation == -1:
				intermediate_pt = terminal_street.initial_point().end_point
			else:
				raise ValueError
			if terminal_street != street:
				new_street_orientation = set_orientation_from_starting_point(
						street, intermediate_pt, slot
					)
				self.path.append([street, new_street_orientation])
			else:
				new_street_orientation = -set_orientation_from_starting_point(
						street, intermediate_pt, slot
					)
				self.path.append([street, new_street_orientation])

		elif end_pt == 'first' and (
				self.growth_restriction == 'backward_only' or 
				self.growth_restriction is None
			):
			terminal_street, terminal_orientation = self.path[0]
			if terminal_orientation == +1:
				intermediate_pt = terminal_street.initial_point().end_point
			elif terminal_orientation == -1:
				intermediate_pt = terminal_street.final_point().end_point
			else:
				raise ValueError
			if terminal_street != street:
				new_street_orientation = set_orientation_from_starting_point(
						street, intermediate_pt, slot
					)
				self.path.insert(0, [street, new_street_orientation])
			else:
				new_street_orientation = -set_orientation_from_starting_point(
						street, intermediate_pt, slot
					)
				self.path.insert(0, [street, new_street_orientation])

		else:
			raise Exception(
				'It appears that the dash {} cannot be extended along '
				'street {} from the {} endpoint. Note that the growth'
				'restriction is {}'
				.format(
					self.label, street.label, end_pt, self.growth_restriction
				)
			)

	def initial_street(self):
		return self.path[0][0]

	def final_street(self):
		return self.path[-1][0]

	def starting_point(self):
		if self.path[0][1] == 1:
			return DashEndpoint(
				dash=self,
				street_end_pt=self.initial_street().initial_point(),
				orientation='out'
			)
		elif self.path[0][1] == -1:
			return DashEndpoint(
				dash=self,
				street_end_pt=self.initial_street().final_point(),
				orientation='out'
			)
		else:
			raise ValueError

	def ending_point(self):
		if self.path[-1][1] == 1:
			return DashEndpoint(
				dash=self,
				street_end_pt=self.final_street().final_point(),
				orientation='in'
			)
		elif self.path[-1][1] == -1:
			return DashEndpoint(
				dash=self,
				street_end_pt=self.final_street().initial_point(),
				orientation='in'
			)
		else:
			raise ValueError

	def print_endpoints(self):
		print (
			'Path {} starts from {} at slot {}, going out on street {}, '
			'and ends on {} at slot {} arriving on street {}.'
			.format(
				self.label, 
				self.starting_point().end_point.label, 
				self.starting_point().slot,
				self.starting_point().street.label,
				self.ending_point().end_point.label,
				self.ending_point().slot,
				self.ending_point().street.label,
			)
		)


class DashEndpoint:
	def __init__(self, dash=None, street_end_pt=None, orientation=None):
		self.dash = dash
		self.end_point = street_end_pt.end_point
		self.street = street_end_pt.street
		self.slot = street_end_pt.slot
		self.street_end_point = street_end_pt
		self.orientation = orientation


class SolitonPath:
	"""
	The attribute .streets_set just tracks the homology of a soliton.
	But the attribute .path contains full information
	about its path as it winds on the Seiberg-Witten curve.

	A soliton path must be constructed by propagation through the network. 
	One starts with a pair of dashes, above a certain street of the network.
	Then dashes can be extended by propagation, following network rules,
	either across joints or branch points.
	When all dashes end up uniting into a single one, 
	the soliton is completed.
	A completed soliton is characterized by the attribure .is_complete = True.

	Note: the dashes of a soliton will include the lift of the 
	initial street on which it is sourced.
	Therefore, when gluing two solitons to form a closed one,  
	one must remove one copy of the lift of the initial street on which 
	both solitons are supported.

	When growing a soliton, the number of dashes can temporarily increase.
	This will happen for example at a 3-way joint.
	There will be two dashes initially (sheets i and j) and one more dash
	on sheet k eventually.
	It's important to grow dashes in pairs, not singularly: this is how
	the network rules work.
	For example, suppose we have two dashes d1, d2 above one street: 
	their orientations will be opposite. 
	Since they must both grow in the same direction, one of them will grow
	"forward" and the other "in reverse".
	When they split at a Y-joint, a third dash d3 is created.
	The new dash will grow forward in pair with d1, but also
	and backwards in pair with d2.
	Note that now d1 and d2 are no longer paired.
	"""
	def __init__(self, label=None):
		self.label = label
		self.growing_pairs = []
		self.starting_point = None
		self.ending_point = None
		self.is_complete = False

	def create(self, street=None, source_pt=None, slot=None):
		"""
		The source point and slot characterize the direction of growth
		of each dash.
		The variable source_pt must be a branch point or a joint, while
		slot must be an integer.
		
		If a soliton is stretching "from a joint/branch point" e
		of a street p, then we create an initial dash d_i that 
		is given by the street itself with orientation into e, 
		moreover we impose the constraint that such initial dash 
		can only grow forward.
		Likewise we also create a final dash d_f that 
		is given by the street itself with orientation outgoing from e, 
		such final dash can only grow backward.
		"""
		# First of all, check that the slot actually corresponds to
		# one here the street ands on the starting_point
		if (
			street.initial_point().end_point == source_pt and
			street.initial_point().slot == slot
		): 
			other_endpoint = street.final_point().end_point
			other_slot = street.final_point().slot
		elif (
			street.final_point().end_point == source_pt and
			street.final_point().slot == slot
		):
			other_endpoint = street.initial_point().end_point
			other_slot = street.initial_point().slot
		else:
			raise Exception(
				'Street {} doesnt end on {} at slot {}'
				.format(street.label, starting_pt.label, starting_slot)
			)

		d_i = Dash(label='initial_dash', growth_restriction='forward_only')
		d_f = Dash(label='final_dash', growth_restriction='backward_only')
		
		d_i.extend_dash_along_street(
			street=street, end_pt=other_endpoint, slot=other_slot
		)
		d_f.extend_dash_along_street(
			street=street, end_pt=source_pt, slot=slot
		)
		self.growing_pairs.append([d_i.ending_point(), d_f.starting_point()])
		self.check_growing_pairs()
		self.starting_point = d_i.starting_point()
		self.ending_point = d_f.ending_point()

	def dashes(self):
		dashes = []
		for p in self.growing_pairs:
			if not p[0].dash in dashes:
				dashes.append(p[0].dash)
			if not p[1].dash in dashes:
				dashes.append(p[1].dash)
		return dashes

	def street_set(self):
		raise NotImplemented

	def path(self):
		raise NotImplemented

	def print_growing_pairs(self):
		print 'The growing pairs are:'
		for p in self.growing_pairs:
			print (
				'{} at slot {} with orientation {}, '
				'and {} at slot {} with orientation {}'.
				format(
					p[0].end_point.label, p[0].slot, p[0].orientation, 
					p[1].end_point.label, p[1].slot, p[1].orientation, 
				)
			)

	def check_growing_pairs(self):
		for p in self.growing_pairs:
			if (
				p[0].orientation=='in' and p[1].orientation=='out' or
				p[1].orientation=='in' and p[0].orientation=='out'
			):
				pass
			else:
				raise Exception('A growing pair is not of the (in,out) type.')


# TO DO: develop this class.
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


def set_orientation_from_starting_point(
		street, starting_pt, starting_slot, opposite=False
	):
	"""
	The variable starting_pt must be a branch point or a joint,
	while starting_slot must be an integer.
	"""
	# First of all, check that the slot actually corresponds to
	# one here the street ands on the starting_point
	if (
		street.initial_point().end_point == starting_pt and
		street.initial_point().slot == starting_slot
	) or (
		street.final_point().end_point == starting_pt and
		street.final_point().slot == starting_slot
	):
		pass
	else:
		raise Exception(
			'Street {} doesnt end on {} at slot {}'
			.format(street.label, starting_pt.label, starting_slot)
		)

	if opposite is False:
		coeff = +1
	elif opposite is True:
		coeff = -1

	if (
		starting_pt == street.initial_point().end_point and
		starting_slot == street.initial_point().slot
	):
		return +1 * coeff
	elif (
		starting_pt == street.final_point().end_point and
		starting_slot == street.final_point().slot	
	):
		return -1 * coeff
	else:
		raise Exception(
			'The point {} is not an endpoint of street {} with slot {}'
			.format(starting_pt.label, street.label, starting_slot)
		)


w = MCSN()
# w.check_network()
w.attach_streets()
w.check_network()
s1 = w.streets['p_1']
e11 = s1.initial_point().end_point
e12 = s1.final_point().end_point
s2 = w.streets['p_2']
e21 = s2.initial_point().end_point
e22 = s2.final_point().end_point
s3 = w.streets['p_3']
e31 = s3.initial_point().end_point
e32 = s3.final_point().end_point
j1 = w.joints.values()[0]
b1 = w.branch_points.values()[0]
b1.print_type()
j1.print_type()

dash1 = Dash(label='dash_1')
dash1.extend_dash_along_street(street=s1, end_pt=e11, slot=0)
dash1.print_endpoints()
dash1.extend_dash_along_street(
	street=s2, end_pt='last', slot=j1.street_position(s2)[0]
)
dash1.print_endpoints()
dash1.extend_dash_along_street(
	street=s1, end_pt='first', slot=b1.street_position(s1)[0]
)
dash1.print_endpoints()

a1 = SolitonPath(label='a_1')
a1.create(street=s1, source_pt=b1, slot=0)
a1.print_growing_pairs()
print [d.label for d in a1.dashes()]

