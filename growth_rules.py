from solitons import copy_of_soliton, join_dashes

"""
Growth of a soliton works by bifurcation of possibilities.
Each soliton has several growing pairs.
Two types of things can happen:
1- two or more growing pairs can merge into a single one (or even more??)
2- a growing pair keeps growing

Option 1 arises when two or more growing pairs arrive on the same 
branch point / joint.
In this case, the rules of a joint /branch point will determine
if the two (or more) growing pairs can be merged in several ways.
In this case the function must:
1a- connect one or more pairs of dashes (pairwise)
1b- form a new set of pairs of endpoints, and return those

Option 2 is simpler to handle.
Focus on a single growing pair of endpoints for a soliton.
A pair of endpoints is given at a certain joint or 
branch point. Then, by applying a certain rule,
which will depdend on the branch point / joint in 
question and on which of its slots the pair ends, 
the function must 
2a- grow the two Dashes of that pair of endpoints
2b- create new Dashes if necessary according to the rules
2c- form a new set of pairs of endpoints, and return those

To determine the various possibilities, one first looks at the
whole set of growing pairs. 
One makes a "map" of which slots are taken, for each branch/point and joint.
Then, a function for each of these branch points or joints will
determine how many possibilities there are, at each site.
Note that more than one pair of endpoint may end on the same slot, 
each must be considered separately in all possible combinations.
All possibilities are considered, generating a number of new solitons.
"""

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
	print 'Nodal points are {}'.format(nodal_points)

	clusters = []
	for n in nodal_points:
		# select all the pairs that end on the nodal point
		nodal_cluster = [p for p in growing_pairs if p[0].end_point == n]
		print 'For nodal point {}, the pairs are {}'.format(
			n.label, nodal_cluster
		)
		# note that there may be two or more pairs in a nodal cluster 
		# that end on the same slot.
		clusters.append(nodal_cluster)

	return clusters


def grow_soliton(soliton):
	"""
	This function must return several new solitons.
	New solitons are created as copies of the original one,
	and their evolution is outsourced to external functions.
	Since growth must be performed on copies, the directions 
	given to external functions must be provided in terms of
	relative data, not actual data.
	Concretely: we cannot ask a certain growing pair to be 
	grown at a certain joint. Because that growing pair
	would belong to the original soliton.
	Instead, we must first copy the soliton, then identify 
	the corresponding growing pair that matches to the one
	in the original soliton, and then ask for that to be grown.
	
	Moreover, instead of growing pairs, we must consider 
	growing clusters: there may be more than a pair of dashes 
	that ends on a joint or branch point, for a certain time
	during soliton growth.
	Two or more pairs of dashes can furthermore join together,
	into a smaller number of pairs. 
	QUESTION: [IS THIS ACTUALLY CORRECT? THINK ABOUT \TAU AND \NU'S]
	The way this can happen is handled by the rules at each 
	joint or branch point.
	"""
	growing_pairs = soliton.growing_pairs
	clusters = growing_clusters(soliton)
	new_solitons = []

	for c_i, c in enumerate(clusters):
		# slot taken by each growing pair
		c_slots = [pair[0].slot for pair in c]
		node = c[0][0].end_point
		node_type = node.type
		av_slots = node.available_slots
		for s in av_slots:
			# check if a slot is taken by more than one growing pair
			if c_slots.count(s) > 1:
				# probably this kind of evolution will be handled directly 
				# by each branch point or joint. To be decided.
				raise NotImplementedError
			else:
				pass

		if node_type == 'type_1_branch_point':
			# here c_i tells which cluster must be grown
			# we must not pass the actual cluster,
			# because growth will not be performed on this soliton, 
			# but rather on a copy of it.
			new_solitons.join(bp_type_1_growth(soliton, c))
		elif node_type == 'type_3_joint':
			new_solitons.join(j_type_3_growth(soliton, c_i))
		else:
			raise NotImplementedError


def bp_type_1_growth(old_soliton, old_cluster):
	"""
	Branch points of type 1 are those with a single two-way street 
	ending on them, so all that can happen is that each pair of
	dashes "caps off" there.
	There can be no mixing between different pairs, this follows
	from the algebra of tau and nu generating functions of 2d solitons.
	"""
	new_solitons = []
	for p in old_cluster:
		new_soliton = copy_of_soliton(old_soliton)
		# Now for the growing pair p of DashEndpoints we 
		# identify the new corresponding growing pair in 
		# the new soliton
		old_pair_index = old_soliton.growing_pairs.index(p)
		new_p = new_soliton.growing_pairs[old_pair_index]
		# Then, we create a new dash, by joining the dashes of this 
		# pair at the branch point. This is where the rule for a branch
		# point of type 1 enters. 
		new_dash = join_dashes(new_p)
		# Remove the growing pair that has been joined at the branch point.
		new_growing_pairs = (
			[new_p_i for new_p_i in new_soliton.growing_pairs 
			if new_p_i != new_p]
		)
		new_soliton.growing_pairs = new_growing_pairs

		print '\nThese are the soliton starting point and ending point'
		print new_soliton.starting_point
		print new_soliton.ending_point
		# print 'But notice that these are the new dashes points'
		# print new_p[0]
		# print new_p[1]

		if len(new_growing_pairs) == 0: 
			# In this case the soliton has been completed
			new_soliton.is_complete = True
			new_soliton.complete_dash = new_dash
			new_soliton.starting_point = new_dash.starting_point
			new_soliton.ending_point = new_dash.ending_point
		else:
			# If the dash that is being replaced is initial or 
			# final dash, must also replace initial and 
			# final points of the soliton
			### TODO: this needs some testing
			dash_endpoint_1 = new_p[0]
			dash_endpoint_2 = new_p[0]
			if (
				dash_endpoint_1 == new_soliton.starting_point or
				dash_endpoint_2 == new_soliton.starting_point
			):
				new_soliton.starting_point = new_dash.starting_point
			if (
				dash_endpoint_1 == new_soliton.ending_point or
				dash_endpoint_2 == new_soliton.ending_point
			):
				new_soliton.ending_point = new_dash.ending_point

		new_solitons.append(new_soliton)
		### QUESTION: in the network rules, there seems to be a precise order
		### of dashes. [REALLY??? No: each dash is a string of \tau and \nu
		### and these strings are just summed together like
		### \tau = \nu\nu\nu+\tau\nu+\nu\tau etc] 
		### That moreover affects how they can compose.
		### Should keep the ordering into account as well, and when creating 
		### a new dash make sure that the new ordering is correct.
		### It's important to pin down carefully the relation between dashes 
		### and the 2d soliton generating functions, tau and nu.
		
	return new_solitons
		


def bp_type_2_growth():
	raise NotImplementedError


def bp_type_3_growth():
	raise NotImplementedError


def j_type_3_growth():
	raise NotImplementedError


def j_type_4_A_growth():
	raise NotImplementedError

from solitons import a1, b1, j1

print '\nthe soliton dashes are \n{}'.format([d.label for d in a1.dashes()])
print 'growth restrictions are {}\n'.format([d.growth_restriction for d in a1.dashes()])

a1.dashes()[0].print_endpoints()
print 'starting point {}'.format(a1.dashes()[0].starting_point)
print 'ending point {}\n'.format(a1.dashes()[0].ending_point)

a1.dashes()[1].print_endpoints()
print 'starting point {}'.format(a1.dashes()[1].starting_point)
print 'ending point {}'.format(a1.dashes()[1].ending_point)

print '\nsoliton starting point {}'.format(a1.starting_point)
print 'soliton ending point {}'.format(a1.ending_point)

# print '\navailable slots at branch point b1'
# print b1.available_slots
# print '\navailable slots at joint j1'
# print j1.available_slots


cls = growing_clusters(a1)
print '\nThe growing clusters:\n{}'.format(cls)
print '\nNow I grow soliton a1 at {}'.format(cls[0][0][0].end_point.label)
new_sols = bp_type_1_growth(a1, cls[0])
# print len(new_sols)
a2=new_sols[0]

print '\nthe soliton dashes are \n{}'.format([d.label for d in a2.dashes()])
print 'growth restrictions are {}\n'.format([d.growth_restriction for d in a2.dashes()])

a2.dashes()[0].print_endpoints()
print 'starting point {}'.format(a2.dashes()[0].starting_point)
print 'ending point {}'.format(a2.dashes()[0].ending_point)


print '\nThe soliton is now complete: {}'.format(a2.is_complete)
print (
	'The soliton starting point equals the dash starting point: {}'
	.format(a2.starting_point==a2.complete_dash.starting_point)
)
print (
	'The soliton ending point equals the dash ending point: {}'
	.format(a2.ending_point==a2.complete_dash.ending_point)
)




