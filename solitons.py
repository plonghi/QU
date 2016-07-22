from copy import deepcopy
from mcsn import MCSN, Street, Joint, BranchPoint

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
	backward_only, both (dash is a completed soliton), 
	or None (if can grow both ways).

	The starting point and ending point of a dash are important objects.
	They are attributes of solitons.
	So when a dash is grown, it is important to NOT change the dash
	endpoint, but rather to update its attributes to reflect the 
	new path of the dash.
	"""
	def __init__(self, label='no_label', growth_restriction=None, path=None):
		self.label = label
		if path is None:
			self.path = []
		else:
			self.path=path
		if growth_restriction is not None:
			self.growth_restriction = growth_restriction
		else:
			self.growth_restriction = None
		self.starting_point = self.determine_starting_point()
		self.ending_point = self.determine_ending_point()
	
	def extend_dash_along_street(self, street=None, end_pt=None, slot=None):
		"""
		Extends a dash along a street, starting from a joint/branch point
		corresponding to end_pt, at the given slot (need to specify, 
		as a street may end on more than one slot on the same end_pt)
		"""
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

		elif self.growth_restriction == 'both':
			# The dash cannot be grown anymore, 
			# it should correspond to a completed soliton.
			pass
		
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
				# new_street_orientation = set_orientation_from_starting_point(
				# 		street, intermediate_pt, slot
				# 	)
				# self.path.append([street, new_street_orientation])

				# since 'end_pt' is 'last', in this case the dash goes 
				# INTO the NEW street,
				# therefore the intermediate point is a "starting point"
				# for the new street, according to the orientation
				new_street_orientation = (
					set_orientation_from_starting_point(
						street, intermediate_pt, slot
					)
				)
				self.path.append([street, new_street_orientation])
			else:
				raise Exception('Review code here before proceeding.')
				new_street_orientation = -set_orientation_from_starting_point(
						street, intermediate_pt, slot
					)
				self.path.append([street, new_street_orientation])

		elif end_pt == 'first' and (
				self.growth_restriction == 'backward_only' or 
				self.growth_restriction is None
			):
			# recover info on the terminal street of the path
			# (either the first street or the last one)
			# and whether the street orientation agrees or disagrees
			# with the direction of the dash
			terminal_street, terminal_orientation = self.path[0]
			if terminal_orientation == +1:
				intermediate_pt = terminal_street.initial_point().end_point
			elif terminal_orientation == -1:
				intermediate_pt = terminal_street.final_point().end_point
			else:
				raise ValueError
			if terminal_street != street:
				# if :
				# 	# In this case the dash goes INTO the NEW street,
				# 	# therefore the intermediate point is a "starting point"
				# 	# for the new street, according to the orientation
				# 	new_street_orientation = (
				# 		set_orientation_from_starting_point(
				# 			street, intermediate_pt, slot
				# 		)
				# 	)
				# elif terminal_orientation == -1:
				# 	# In this case the dash goes OUT from the NEW street,
				# 	# therefore the intermediate point is an "ending point"
				# 	# for the new street, according to the orientation
				# 	new_street_orientation = (
				# 		set_orientation_from_starting_point(
				# 			street, intermediate_pt, slot, opposite=True
				# 		)
				# 	)

				# since 'end_pt' is 'first', in this case the dash goes 
				# OUT from the NEW street,
				# therefore the intermediate point is an "ending point"
				# for the new street, according to the orientation
				new_street_orientation = (
					set_orientation_from_starting_point(
						street, intermediate_pt, slot, opposite=True
					)
				)
				self.path.insert(0, [street, new_street_orientation])
			else:
				raise Exception('Review code here before proceeding.')
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
		# update endpoints of the dash
		self.update_endpoints()

	def initial_street(self):
		return self.path[0][0]

	def final_street(self):
		return self.path[-1][0]

	def determine_starting_point(self):
		if len(self.path)==0:
			return None
		else:
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

	def determine_ending_point(self):
		if len(self.path)==0:
			return None
		else:
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

	def update_endpoints(self):
		old_start = self.starting_point
		old_end = self.ending_point
		new_start = self.determine_starting_point()
		new_end = self.determine_ending_point()

		if old_start==None:
			self.starting_point = new_start
		elif (
			old_start.end_point!=new_start.end_point or 
			old_start.slot!=new_start.slot
		):
			# note that we don't replace the DashEndpoint instance,
			# but just update its attribute
			self.starting_point.dash = new_start.dash
			self.starting_point.end_point = new_start.end_point
			self.starting_point.street = new_start.street
			self.starting_point.slot = new_start.slot
			self.starting_point.street_end_point = new_start.street_end_point
			self.starting_point.orientation = new_start.orientation
		else:
			pass

		if old_end==None:
			self.ending_point = new_end
		elif (
			old_end.end_point!=new_end.end_point or 
			old_end.slot!=new_end.slot
		):
			# note that we don't replace the DashEndpoint instance,
			# but just update its attribute
			self.ending_point.dash = new_end.dash
			self.ending_point.end_point = new_end.end_point
			self.ending_point.street = new_end.street
			self.ending_point.slot = new_end.slot
			self.ending_point.street_end_point = new_end.street_end_point
			self.ending_point.orientation = new_end.orientation
		else:
			pass

	def print_endpoints(self, text=False):
		if len(self.path)==0:
			print 'The dash is empty.'
		elif text is True:
			print (
				'Path {} starts from {} at slot {}, going out on street {}, '
				'and ends on {} at slot {} arriving on street {}.'
				.format(
					self.label, 
					self.starting_point.end_point.label, 
					self.starting_point.slot,
					self.starting_point.street.label,
					self.ending_point.end_point.label,
					self.ending_point.slot,
					self.ending_point.street.label,
				)
			)
		else:
			print (
				'Path {} :\n\t({}, slot {}, street {}) ---> ({}, slot {}, street {}), '
				.format(
					self.label, 
					self.starting_point.end_point.label, 
					self.starting_point.slot,
					self.starting_point.street.label,
					self.ending_point.end_point.label,
					self.ending_point.slot,
					self.ending_point.street.label,
				)
			)

	def print_path_info(self):
		print [[street.label, sign] for street, sign in self.path]


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
	# TODO: elimiate growing_pairs form the input, just keep 
	# dashes, and use a function to determine the growing pairs.
	def __init__(
		self, label=None, starting_point=None, 
		ending_point=None, growing_pairs=None,
		dashes=None, complete_dash=None,
	):
		self.label = label
		if growing_pairs is None:
			self.growing_pairs = []
		else:
			self.growing_pairs = growing_pairs
		if dashes is None:
			self.dashes = []
		else:
			self.dashes = dashes
		# self.starting_point = starting_point
		# self.ending_point = ending_point
		self.is_complete = False
		self.complete_dash = complete_dash

	def create(self, street=None, source_pt=None, slot=None):
		"""
		The source point and slot characterize the direction of growth
		of each dash.
		The variable source_pt must be a branch point or a joint, while
		slot must be an integer.
		
		If a soliton is stretching "from a joint/branch point" e
		of a street p, then we create an INITIAL dash d_i that 
		is given by the street itself with orientation into e, 
		moreover we impose the constraint that such initial dash 
		can only grow FORWARD.
		Likewise we also create a FINAL dash d_f that 
		is given by the street itself with orientation OUTGOING from e, 
		such final dash can only grow BACKWARD.

		The convention is that the soliton endpoints are fixed at the 
		beginning, as they determine its relative homology class.
		The growth constraints then indeed preserve the initial point
		and the final point.
		As the soliton is created, it is not a closed soliton, even if 
		it is sourced at a branch point. The growth process
		will take care of closing it, or propagating it further through the
		branch point, if this is of type II or III, following the traffic
		rules of GMN5.

		The dashes of a soliton will be stored in a list, their order 
		must reflect the ptoper time parametrization of the soliton.
		Likewise, the "growing pairs" must reflect the (matching) endpoints
		of consective dashes.
		"""
		# First of all, check that the slot actually corresponds to
		# one where the street ands on the starting_point
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
		self.dashes = [d_i, d_f]
		# self.growing_pairs.append([d_i.ending_point, d_f.starting_point])
		self.growing_pairs = growing_pairs_from_dashes(self.dashes)
		self.check_growing_pairs()

	# def dashes_from_growing_pairs(self):
	# 	dashes = []
	# 	for p in self.growing_pairs:
	# 		if not p[0].dash in dashes:
	# 			dashes.append(p[0].dash)
	# 		if not p[1].dash in dashes:
	# 			dashes.append(p[1].dash)
	# 	if self.is_complete is True:
	# 		dashes.append(self.complete_dash)
	# 	return dashes

	def street_set(self):
		raise NotImplementedError

	def path(self):
		raise NotImplementedError

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

	def print_info(self, full_path=False):
		print 'Dashes of soliton path {} / {}:'.format(self.label, self)
		for d in self.dashes:
			d.print_endpoints()
		if self.is_complete is True:
			print 'This soliton is complete'
		else:
			print 'This soliton is not complete'

		if full_path is True:
			if self.is_complete is True:
				print 'The full path of the soliton is'
				self.complete_dash.print_path_info()

	def check_growing_pairs(self):
		for p in self.growing_pairs:
			if (
				p[0].orientation=='in' and p[1].orientation=='out' or
				p[1].orientation=='in' and p[0].orientation=='out'
			) and (
				p[0].end_point==p[1].end_point
			):
				pass
			else:
				raise Exception('A growing pair is not of the (in,out) type.')

	def check_complete(self):
		dashes = self.dashes()
		for d in dashes:
			if (
				d.starting_point() == self.starting_point() and
				d.ending_point() == self.ending_point()
			):
				if len(dashes) == 1:
					return True
				else: 
					return Exception(
						'The soliton has a complete dash, but also some'
						'disconnected pieces.'
					)
			else:
				return False

	def starting_point(self):
		check_dashes_ordering(self.dashes)
		return self.dashes[0].starting_point

	def ending_point(self):
		check_dashes_ordering(self.dashes)
		return self.dashes[-1].ending_point


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


def copy_of_soliton(soliton, label=None):
	"""
	This function makes a copy of a soliton.
	This includes copying the dashes and preparing a new soliton.
	Note that one cannot just use deeocopy, otherwise the branch points
	and the joints of the network will be copied too.
	That would prevent us from ever joining solitons.
	"""
	dashes = soliton.dashes
	new_dashes = (
		[Dash(
			growth_restriction=d.growth_restriction, path=d.path
		) for d in dashes]
	)

	new_growing_pairs = growing_pairs_from_dashes(new_dashes)

	# now determine the new starting and ending points
	j_0 = dashes.index(soliton.starting_point().dash)
	j_1 = dashes.index(soliton.ending_point().dash)
	new_starting_point = new_dashes[j_0].starting_point
	new_ending_point = new_dashes[j_1].ending_point

	return SolitonPath(
		label=label,
		starting_point=new_starting_point, 
		ending_point=new_ending_point, 
		growing_pairs=new_growing_pairs,
		dashes=new_dashes,
	)
	# new_soliton = deepcopy(soliton)
	# new_soliton.label = label
	# return new_soliton


def join_dashes(growing_pair):
	"""
	Returns a new dash.
	"""
	dash_endpoint_1 = growing_pair[0]
	dash_endpoint_2 = growing_pair[1]

	if dash_endpoint_1.end_point != dash_endpoint_2.end_point:
		raise Exception('Growing pair cannot be joined.')
	# dashes must have a common endpoint to be joined, this is checked above.
	# then, the first dash will be the one whose orientation is 'in'
	# and the second dash will be the one whose orientation is 'out'
	# of the common point
	if (
		dash_endpoint_1.orientation == 'in' and
		dash_endpoint_2.orientation == 'out'
	):
		first_dash = dash_endpoint_1.dash
		second_dash = dash_endpoint_2.dash
	elif (
		dash_endpoint_1.orientation == 'out' and
		dash_endpoint_2.orientation == 'in'
	):
		first_dash = dash_endpoint_2.dash
		second_dash = dash_endpoint_1.dash
	else:
		raise Exception('Growing pair cannot be joined')

	# Now determine growth restrictions
	if (
		first_dash.growth_restriction is 'forward_only' and
		second_dash.growth_restriction is 'backward_only'
	):
		new_growth_restriction = 'both'
	elif (
		first_dash.growth_restriction is 'forward_only' and
		second_dash.growth_restriction is None
	):
		new_growth_restriction = 'forward_only'
	elif (
		first_dash.growth_restriction is None and
		second_dash.growth_restriction is 'backward_only'
	):
		new_growth_restriction = 'backward_only'
	elif (
		first_dash.growth_restriction is None and
		second_dash.growth_restriction is None
	):
		new_growth_restriction = None
	else:
		raise Exception('Cannot determine growth restriction')

	new_dash = Dash(
		# label=first_dash.label+'_'+second_dash.label, 
		growth_restriction=new_growth_restriction, 
		path=(first_dash.path + second_dash.path)
	)

	return new_dash


def check_dashes_ordering(dash_sequence):
	"""
	Checks an ORDERED sequence of dashes (typically from a soliton)
	such that the first one is "forward only", and its ending point 
	is on the same nodal point (branch pt/joint), same slot as the 
	beginning point of the second one, with opposite orientation.
	And so on for the ending point of the 2nd dash with the beginning 
	point of the 3rd dash, and so on.
	The last dash must be "backward only" type.
	"""
	if (
		dash_sequence[0].growth_restriction is not 'forward_only' and 
		dash_sequence[0].growth_restriction is not 'both'
	):
		raise Exception
	if (
		dash_sequence[-1].growth_restriction is not 'backward_only' and 
		dash_sequence[-1].growth_restriction is not 'both'
	):
		raise Exception
	for i in range(len(dash_sequence)-1):
		d_1 = dash_sequence[i]
		d_2 = dash_sequence[i+1]
		if (
			d_1.ending_point.end_point != d_2.starting_point.end_point or
			d_1.ending_point.slot != d_2.starting_point.slot
		):
			print '\nDashes are not ordered properly:'
			for d in dash_sequence:
				d.print_endpoints()
			raise Exception
	pass


def growing_pairs_from_dashes(dash_sequence):
	"""
	This function returns the corresponding set of growing pairs.
	"""
	check_dashes_ordering(dash_sequence)

	growing_pairs = []
	for i in range(len(dash_sequence)-1):
		d_1 = dash_sequence[i]
		d_2 = dash_sequence[i+1]
		growing_pairs.append([d_1.ending_point, d_2.starting_point])

	return growing_pairs


def growing_clusters(soliton):
	"""
	This function returns a list of lists of growing pairs.
	Each sub-list is a set of growing pairs attached to the same 
	joint or branch point.
	It may happen that two or more growing pairs occupy the same slot,
	within a single sub-list.
	"""
	growing_pairs = soliton.growing_pairs
	nodal_points = []	# a list of branch points or joints
	for p in growing_pairs:
		if p[0].end_point not in nodal_points:
			nodal_points.append(p[0].end_point)
		# The following is superfluous, if growing pairs are correctly formed
		# if p[1].end_point not in nodal_points:
		# 	nodal_points.append(p[1].end_point)
	# print 'Nodal points are {}'.format(nodal_points)

	clusters = []
	for n in nodal_points:
		# select all the pairs that end on the nodal point
		nodal_cluster = [p for p in growing_pairs if p[0].end_point == n]
		# print 'For nodal point {}, the pairs are {}'.format(
		# 	n.label, nodal_cluster
		# )
		# note that there may be two or more pairs in a nodal cluster 
		# that end on the same slot.
		clusters.append(nodal_cluster)

	return clusters


def find_corresponding_cluster(soliton, ref_cluster):
	"""
	For a given soliton, identifies which of its growth clusters
	corresponds to a certain 'reference cluster'.
	"""
	# characterize the reference cluster by 
	# - the nodal point (joint/branch point)
	# - the slot on the nodal point
	ref_gr_pair = ref_cluster[0]
	ref_end_point = ref_gr_pair[0].end_point
	ref_slot = ref_gr_pair[0].slot

	clusters = growing_clusters(soliton)
	candidates = []
	# collect all candidate clusters, based on their nodal point and slot
	for cl in clusters:
		sample_gr_pair = cl[0]
		if (
			sample_gr_pair[0].end_point == ref_end_point and
			sample_gr_pair[0].slot == ref_slot
		):
			candidates.append(cl)

	if len(candidates) > 1:
		raise Exception(
			'Cannot find a unique cluster corresponding to the reference one'
		)
	elif len(candidates) == 0:
		raise Exception(
			'Cannot find any cluster corresponding to the reference one'
		)
	else:
		return candidates[0]





