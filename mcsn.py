ALLOWED_JOINTS = (
	# 3-way
	[
	['Street', 'NoneType', 'Street', 'NoneType', 'Street', 'NoneType'],
	['NoneType', 'Street', 'NoneType', 'Street', 'NoneType', 'Street'],
	] +
	# 4-way - peace sign
	[
	['NoneType', 'Street', 'NoneType', 'Street', 'Street', 'Street'],
	['Street', 'NoneType', 'Street', 'Street', 'Street', 'NoneType'],
	['NoneType', 'Street', 'Street', 'Street', 'NoneType', 'Street'],
	['Street', 'Street', 'Street', 'NoneType', 'Street', 'NoneType'],
	['Street', 'Street', 'NoneType', 'Street', 'NoneType', 'Street'],
	['Street', 'NoneType', 'Street', 'NoneType', 'Street', 'Street'],
	] + 
	# 4-way - X sign
	[
	['Street', 'Street', 'NoneType', 'Street', 'Street', 'NoneType'],
	['Street', 'NoneType', 'Street', 'Street', 'NoneType', 'Street'],
	['NoneType', 'Street', 'Street', 'NoneType', 'Street', 'Street'],
	] +
	# 5-way
	[
	['NoneType', 'Street', 'Street', 'Street', 'Street', 'Street'],
	['Street', 'NoneType', 'Street', 'Street', 'Street', 'Street'],
	['Street', 'Street', 'NoneType', 'Street', 'Street', 'Street'],
	['Street', 'Street', 'Street', 'NoneType', 'Street', 'Street'],
	['Street', 'Street', 'Street', 'Street', 'NoneType', 'Street'],
	['Street', 'Street', 'Street', 'Street', 'Street', 'NoneType'],
	] + 
	# 6-way
	[
	['Street', 'Street', 'Street', 'Street', 'Street', 'Street'],
	]
)


class BranchPoint:
	"""
	All streets ending on a branch point must be
	given when the object is instanced.
	The ordering is important, it is understood to be 
	counter-clockwise, as in 1204.4824 figure 47 on page 74.
	"""
	def __init__(self, label=None, streets=None):
		self.label=label
		if streets is None:
			self.streets = [None for i in range(3)]
		elif len(streets) < 3:
			self.streets = streets + [None for i in range(3-len(streets))]
		elif len(streets) == 3:
			self.streets = streets
		else:
			raise Exception(
				'Cannot handle branch points with more than three streets.'
			)
		pass

	def street_position(self, street):
		"""
		Note, there may be cases in which a street ends twice on a 
		branch point or joint.
		"""
		positions = []
		for i, s in enumerate(self.streets):
			if street == s:
				positions.append(i)
		return positions
		
		# try:
		# 	ind = self.streets.index(street)
		# except ValueError:
		# 	ind = None
		# return ind


class Street:
	def __init__(self, label=None):
		self.label=label
		self.endpoints = []
		pass

	def add_endpoint(self, end_point):
		if len(self.endpoints) == 2:
			raise Exception(
				'Street {} cannot have more than two endpoints.'
				.format(self.label)
			)
		else:
			self.endpoints.append(end_point)

	def initial_point(self):
		return self.endpoints[0]
		
	def final_point(self):
		return self.endpoints[1]


class StreetEndpoint:
	def __init__(self, street=None, end_point=None, slot=None):
		self.street = street
		self.end_point = end_point
		self.slot = slot


class Joint:
	"""
	All streets ending on a joint must be
	given when the object is instanced.
	The ordering is important, it is understood to be 
	counter-clockwise, as in 1204.4824 figure 46 on page 73.
	Very important: All 6 streets must be specified! If a street 
	is one-way, the corresponding entry will be the keyword "None".
	"""
	def __init__(self, label=None, streets=None):
		self.label=label
		if streets is None:
			self.streets = [None for i in range(6)]
		elif len(streets) == 6:
			self.streets = streets
		else:
			raise Exception(
				'Must specify 6 streets for each joint.'
			)

	def street_position(self, street):
		"""
		Note, there may be cases in which a street ends twice on a 
		branch point or joint.
		"""
		positions = []
		for i, s in enumerate(self.streets):
			if street == s:
				positions.append(i)

		return positions

		# try:
		# 	ind = self.streets.index(street)
		# except ValueError:
		# 	ind = None
		# return ind


### Important: must develop this class, 
### to extract the charges to assign to each
### lift of solitons.
class StringUnit:
	def __init__(self, label=None):
		self.label=label
		pass


class MCSN:
	"""
	The network is constructed by specifying 
	the streets attached to each joint and branch point.
	Then, the method 'attach_streets()' is called.
	"""
	def __init__(self, label=None):
		self.label=label
		self.streets = {
			'p_1' : Street(label='p_1'),
			'p_2' : Street(label='p_2'),
			'p_3' : Street(label='p_3'),
		}
		self.branch_points = {
			'b_1' : BranchPoint(
				label='b_1', streets=[self.streets['p_1'], None, None]
			), 
			'b_2' : BranchPoint(
				label='b_2', streets=[self.streets['p_2']]
			),
			'b_3' : BranchPoint(
				label='b_3', streets=[self.streets['p_3']]
			),
		}
		self.joints = {'j_1': Joint(
			label='j_1', streets=[
				self.streets['p_1'], 
				None,
				self.streets['p_2'], 
				None,
				self.streets['p_3'], 
				None
			]
		)}

	def attach_streets(self):
		for b_pt in self.branch_points.values():
			for i, street in enumerate(b_pt.streets):
				if street is not None:
					street.add_endpoint(
						StreetEndpoint(
							street=street, end_point=b_pt, slot=i,
						)
					)
		
		for j_pt in self.joints.values():
			for i, street in enumerate(j_pt.streets):
				if street is not None:
					street.add_endpoint(
						StreetEndpoint(
							street=street, end_point=j_pt, slot=i,
						)
					)

	def check_network(self):
		# Check that every street has exactly two endpoints
		for street in self.streets.values():
			if len(street.endpoints) == 2:
				pass
			else:
				raise Exception(
					'Street {} does not have two endpoints, but {}'
					.format(street.label, street.endpoints)
				)
		print 'All streets have two well-defined endpoints.'

		# Check that every joint is of an allowed type
		for j_pt in self.joints.values():
			j_type = [st.__class__.__name__ for  st in j_pt.streets]
			is_allowed = False
			for a_joint in ALLOWED_JOINTS:
				if j_type == a_joint:
					is_allowed = True
					break
				else:
					pass

			if is_allowed is False:
				raise Exception(
					'Joint {} is of type {}'.format(j_pt.label, j_type)
				)
		print 'All joints are of a well-defined type.'
		
	def add_street(self):
		pass

	def add_joint(self):
		pass

	def add_branch_point(self):
		pass


# w = MCSN()
# # w.check_network()
# w.attach_streets()
# w.check_network()
# s1 = w.joints['j_1'].streets[0]
# print w.joints['j_1'].street_position(w.streets['p_2'])
# print w.branch_points['b_1'].street_position(w.streets['p_2'])
# print s1.initial_point().end_point
# print w.joints['j_1'].street_position(w.streets['p_1'])