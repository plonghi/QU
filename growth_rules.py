from solitons import (
	copy_of_soliton, join_dashes, SolitonPath, Dash, 
	growing_pairs_from_dashes, growing_clusters, find_corresponding_cluster
)

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

def grow_soliton_once(soliton):
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

	Growth will be organized in 'clusters', i.e. all the 
	growing pairs (endpoints of dashes above the same nodal point
	of the network) of a soliton will be organized according
	to which nodal point they lie on.
	All growing pairs above the same nodal point will belong to
	the same cluster, and will be grown at once.
	
	# NOTE: the following is DEPRECATED.
	# it seems unnecessary at this point. But keep for the future.
	# 
	# Two or more pairs of dashes can furthermore join together,
	# into a smaller number of pairs. 
	# QUESTION: [IS THIS ACTUALLY CORRECT? THINK ABOUT \TAU AND \NU'S]
	# The way this can happen is handled by the rules at each 
	# joint or branch point.
	
	"""

	if soliton.is_complete:
		return [soliton]

	else:
		# growing_pairs = soliton.growing_pairs
		original_clusters = growing_clusters(soliton)
		old_solitons = [soliton]
		new_solitons = []

		# DEPRECATED: all nontrivial interactions at a nodal point
		# are taken into account by the traffic rules, there is no 
		# need to check possible combinations of dashes.
		#
		# for c_i, c in enumerate(clusters):
		# 	# slot taken by each growing pair
		# 	c_slots = [pair[0].slot for pair in c]
		# 	node = c[0][0].end_point
		# 	node_type = node.type
		# 	av_slots = node.available_slots
		# 	for s in av_slots:
		# 		# check if a slot is taken by more than one growing pair
		# 		if c_slots.count(s) > 1:
		# 			# probably this kind of evolution will be handled directly 
		# 			# by each branch point or joint. To be decided.
		# 			raise NotImplementedError
		# 		else:
		# 			pass

		# 	if node_type == 'type_1_branch_point':
		# 		# here c_i tells which cluster must be grown
		# 		# we must not pass the actual cluster,
		# 		# because growth will not be performed on this soliton, 
		# 		# but rather on a copy of it.
		# 		new_solitons.join(bp_type_1_growth(soliton, c))
		# 	elif node_type == 'type_3_joint':
		# 		new_solitons.join(j_type_3_growth(soliton, c_i))
		# 	else:
		# 		raise NotImplementedError

		for cl in original_clusters:
			# update the solitons in the list old_solitons, for each growing
			# cluster, returning a new list called 'new_solitons' which
			# will supersed old_solitons, thus being updated at the
			# next growth step (next cluster)

			node = cl[0][0].end_point
			node_type = node.type

			# if node_type == 'type_1_branch_point':
			# 	new_solitons = []
			# 	for sol in old_solitons:
			# 		# must now identify which growth cluster of each soliton 
			# 		# is the one corresponding to 'original_cluster'
			# 		sol_cl = find_corresponding_cluster(sol, cl)
			# 		new_solitons += bp_type_1_growth(sol, sol_cl)
			# 	old_solitons = new_solitons

			# elif node_type == 'type_3_joint':
			# 	new_solitons = []
			# 	for sol in old_solitons:
			# 		# must now identify which growth cluster of each soliton 
			# 		# is the one corresponding to 'original_cluster'
			# 		sol_cl = find_corresponding_cluster(sol, cl)
			# 		new_solitons += j_type_3_growth(sol, sol_cl)
			# 	old_solitons = new_solitons

			# elif node_type == 'type_2_branch_point':
			# 	new_solitons = []
			# 	for sol in old_solitons:
			# 		# must now identify which growth cluster of each soliton 
			# 		# is the one corresponding to 'original_cluster'
			# 		sol_cl = find_corresponding_cluster(sol, cl)
			# 		new_solitons += bp_type_2_growth(sol, sol_cl)
			# 	old_solitons = new_solitons

			# elif node_type == 'type_3_branch_point':
			# 	new_solitons = []
			# 	for sol in old_solitons:
			# 		# must now identify which growth cluster of each soliton 
			# 		# is the one corresponding to 'original_cluster'
			# 		sol_cl = find_corresponding_cluster(sol, cl)
			# 		new_solitons += bp_type_3_growth(sol, sol_cl)
			# 	old_solitons = new_solitons

			# else:
			# 	raise NotImplementedError

			if node_type == 'type_1_branch_point':
				grow_at_node = bp_type_1_growth
			elif node_type == 'type_3_joint':
				grow_at_node = j_type_3_growth
			elif node_type == 'type_2_branch_point':
				grow_at_node = bp_type_2_growth
			elif node_type == 'type_3_branch_point':
				grow_at_node = bp_type_3_growth
			elif node_type == 'type_4_A_joint':
				grow_at_node = j_type_4_A_growth
			else:
				raise NotImplementedError

			new_solitons = []
			for sol in old_solitons:
				# must now identify which growth cluster of each soliton 
				# is the one corresponding to 'original_cluster'
				sol_cl = find_corresponding_cluster(sol, cl)
				new_solitons += grow_at_node(sol, sol_cl)
			old_solitons = new_solitons

		return new_solitons


def grow_soliton(soliton, n_steps=1):
	"""
	Repeated application of grow_soliton_once(soliton)
	"""
	old_solitons = [soliton]
	for i in range(n_steps):
		new_solitons = []
		for sol in old_solitons:
			new_solitons += grow_soliton_once(sol)
		old_solitons = new_solitons

	return new_solitons


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
		#the two dashes that will be merged are the following
		new_d_1 = new_p[0].dash
		new_d_2 = new_p[1].dash

		# Then, we create a new dash, by joining the dashes of this 
		# pair at the branch point. This is where the rule for a branch
		# point of type 1 enters. 
		new_dash = join_dashes(new_p)

		# collect the new dashes, by replacing the two we just merged
		# with their union
		# the dashes are ordered, so the two that have been joined 
		# are consecutive ones, their indices are
		ind_1 = new_soliton.dashes.index(new_d_1)
		ind_2 = new_soliton.dashes.index(new_d_2)
		if (ind_1 + 1 != ind_2) and (ind_1 - 1 != ind_2):
			raise Exception('Dashes to be concatenated are not consecutive.')
		
		new_dashes = (
			new_soliton.dashes[0:min(ind_1, ind_2)] 
			+ [new_dash] 
			+ new_soliton.dashes[(max(ind_1, ind_2)+1):]
		)
		new_growing_pairs = growing_pairs_from_dashes(new_dashes)

		new_soliton.dashes = new_dashes
		new_soliton.growing_pairs = new_growing_pairs

		if len(new_growing_pairs) == 0: 
			# In this case the soliton has been completed
			new_soliton.is_complete = True
			new_soliton.complete_dash = new_dash
		# else:
		# 	# If the dash that is being replaced is initial or 
		# 	# final dash, must also replace initial and 
		# 	# final points of the soliton
		# 	### TODO: this needs some testing
		# 	dash_endpoint_1 = new_p[0]
		# 	dash_endpoint_2 = new_p[1]
		# 	if (
		# 		dash_endpoint_1 == new_soliton.starting_point or
		# 		dash_endpoint_2 == new_soliton.starting_point
		# 	):
		# 		new_soliton.starting_point = new_dash.starting_point
		# 	if (
		# 		dash_endpoint_1 == new_soliton.ending_point or
		# 		dash_endpoint_2 == new_soliton.ending_point
		# 	):
		# 		new_soliton.ending_point = new_dash.ending_point

		new_solitons.append(new_soliton)
		
	return new_solitons


def bp_type_2_growth(old_soliton, old_cluster):
	"""
	Branch points of type 2 are those with two two-way streets 
	ending on them.
	Two things can happen: 
	1- a soliton 'caps off' there, 
	2- it gets propagated on another street ending on the branch point,
	   as explained in equation (A.7) of 1204.4824
	But, note that the second option is only available to solitons
	propagating on one street, not the other.
	"""
	new_solitons = []
	branch_pt = old_cluster[0][0].end_point
	sol_slot = old_cluster[0][0].slot
	av_slots = branch_pt.available_slots

	# First, consider the capping off of the soliton
	# ths is just the same as for type_1 branch points
	for p in old_cluster:
		new_soliton = copy_of_soliton(old_soliton)
		# Now for the growing pair p of DashEndpoints we 
		# identify the new corresponding growing pair in 
		# the new soliton
		old_pair_index = old_soliton.growing_pairs.index(p)
		new_p = new_soliton.growing_pairs[old_pair_index]
		#the two dashes that will be merged are the following
		new_d_1 = new_p[0].dash
		new_d_2 = new_p[1].dash

		# Then, we create a new dash, by joining the dashes of this 
		# pair at the branch point. This is where the rule for a branch
		# point of type 1 enters. 
		new_dash = join_dashes(new_p)

		# collect the new dashes, by replacing the two we just merged
		# with their union
		# the dashes are ordered, so the two that have been joined 
		# are consecutive ones, their indices are
		ind_1 = new_soliton.dashes.index(new_d_1)
		ind_2 = new_soliton.dashes.index(new_d_2)
		if (ind_1 + 1 != ind_2) and (ind_1 - 1 != ind_2):
			raise Exception('Dashes to be concatenated are not consecutive.')
		
		new_dashes = (
			new_soliton.dashes[0:min(ind_1, ind_2)] 
			+ [new_dash] 
			+ new_soliton.dashes[(max(ind_1, ind_2)+1):]
		)
		new_growing_pairs = growing_pairs_from_dashes(new_dashes)

		new_soliton.dashes = new_dashes
		new_soliton.growing_pairs = new_growing_pairs

		if len(new_growing_pairs) == 0: 
			# In this case the soliton has been completed
			new_soliton.is_complete = True
			new_soliton.complete_dash = new_dash

		new_solitons.append(new_soliton)

	# Then, if pertinent, also handle the growth of the soliton through
	# the branch point.
	# It depends on the branch point's slot after the soliton's occupied 
	# slot, counter-clockwise. If this slot is also attached to a street 
	# (i.e. if it appears in available_slots), then the soliton 
	# can propagate as explained in point 2. Otherwise it cannot.
	next_slot = (sol_slot + 1) % 3
	if next_slot in av_slots:
		next_street = branch_pt.streets[next_slot]
		for p in old_cluster:
			new_soliton = copy_of_soliton(old_soliton)
			# Now for the growing pair p of DashEndpoints we 
			# identify the new corresponding growing pair in 
			# the new soliton
			old_pair_index = old_soliton.growing_pairs.index(p)
			new_p = new_soliton.growing_pairs[old_pair_index]

			# now identify the dashes d_in and d_out, and grow them.
			# it will be important to check the orientation
			# of each growth point in the pair
			if new_p[0].orientation=='in' and new_p[1].orientation=='out':
				d_in = new_p[0].dash
				d_out = new_p[1].dash
			elif new_p[0].orientation=='out' and new_p[1].orientation=='in':
				d_in = new_p[1].dash
				d_out = new_p[0].dash
			else:
				raise ValueError

			# Now grow dashes d_in and d_out
			# note that d_in must be grown forward, so 
			# we specify the 'last' point for growth
			d_in.extend_dash_along_street(
				street=next_street, end_pt='last', slot=next_slot
			)
			# instead d_out must be grown backward, so 
			# we specify the 'first' point for growth
			d_out.extend_dash_along_street(
				street=next_street, end_pt='first', slot=next_slot
			)

			# Then update the soliton's growing pairs
			new_soliton.growing_pairs = (
				growing_pairs_from_dashes(new_soliton.dashes)
			)

			new_solitons.append(new_soliton)
		
	return new_solitons


def bp_type_3_growth(old_soliton, old_cluster):
	"""
	Branch points of type 3 are those with three two-way streets 
	ending on them.
	Two things can happen: 
	1- a soliton 'caps off' there, 
	2- it gets propagated on another street ending on the branch point,
	   as explained in equation (A.7) of 1204.4824
	"""
	new_solitons = []
	branch_pt = old_cluster[0][0].end_point
	sol_slot = old_cluster[0][0].slot
	av_slots = branch_pt.available_slots

	# First, consider the capping off of the soliton
	# ths is just the same as for type_1 branch points
	for p in old_cluster:
		new_soliton = copy_of_soliton(old_soliton)
		# Now for the growing pair p of DashEndpoints we 
		# identify the new corresponding growing pair in 
		# the new soliton
		old_pair_index = old_soliton.growing_pairs.index(p)
		new_p = new_soliton.growing_pairs[old_pair_index]
		#the two dashes that will be merged are the following
		new_d_1 = new_p[0].dash
		new_d_2 = new_p[1].dash

		# Then, we create a new dash, by joining the dashes of this 
		# pair at the branch point. This is where the rule for a branch
		# point of type 1 enters. 
		new_dash = join_dashes(new_p)

		# collect the new dashes, by replacing the two we just merged
		# with their union
		# the dashes are ordered, so the two that have been joined 
		# are consecutive ones, their indices are
		ind_1 = new_soliton.dashes.index(new_d_1)
		ind_2 = new_soliton.dashes.index(new_d_2)
		if (ind_1 + 1 != ind_2) and (ind_1 - 1 != ind_2):
			raise Exception('Dashes to be concatenated are not consecutive.')
		
		new_dashes = (
			new_soliton.dashes[0:min(ind_1, ind_2)] 
			+ [new_dash] 
			+ new_soliton.dashes[(max(ind_1, ind_2)+1):]
		)
		new_growing_pairs = growing_pairs_from_dashes(new_dashes)

		new_soliton.dashes = new_dashes
		new_soliton.growing_pairs = new_growing_pairs

		if len(new_growing_pairs) == 0: 
			# In this case the soliton has been completed
			new_soliton.is_complete = True
			new_soliton.complete_dash = new_dash

		new_solitons.append(new_soliton)

	# Then, if pertinent, also handle the growth of the soliton through
	# the branch point.
	# It depends on the branch point's slot after the soliton's occupied 
	# slot, counter-clockwise. If this slot is also attached to a street 
	# (i.e. if it appears in available_slots), then the soliton 
	# can propagate as explained in point 2. Otherwise it cannot.
	next_slot = (sol_slot + 1) % 3
	next_street = branch_pt.streets[next_slot]
	for p in old_cluster:
		new_soliton = copy_of_soliton(old_soliton)
		# Now for the growing pair p of DashEndpoints we 
		# identify the new corresponding growing pair in 
		# the new soliton
		old_pair_index = old_soliton.growing_pairs.index(p)
		new_p = new_soliton.growing_pairs[old_pair_index]

		# now identify the dashes d_in and d_out, and grow them.
		# it will be important to check the orientation
		# of each growth point in the pair
		if new_p[0].orientation=='in' and new_p[1].orientation=='out':
			d_in = new_p[0].dash
			d_out = new_p[1].dash
		elif new_p[0].orientation=='out' and new_p[1].orientation=='in':
			d_in = new_p[1].dash
			d_out = new_p[0].dash
		else:
			raise ValueError

		# Now grow dashes d_in and d_out
		# note that d_in must be grown forward, so 
		# we specify the 'last' point for growth
		d_in.extend_dash_along_street(
			street=next_street, end_pt='last', slot=next_slot
		)
		# instead d_out must be grown backward, so 
		# we specify the 'first' point for growth
		d_out.extend_dash_along_street(
			street=next_street, end_pt='first', slot=next_slot
		)

		# Then update the soliton's growing pairs
		new_soliton.growing_pairs = (
			growing_pairs_from_dashes(new_soliton.dashes)
		)

		new_solitons.append(new_soliton)
		
	return new_solitons


def j_type_3_growth(old_soliton, old_cluster):
	"""
	Joints of type 3 are those with a three two-way streets p1, p2, p3
	ending on them, so when a soliton ends on them from one of 
	the streets (say p1), it will keep growing on both the other two, p2, p3.
	For definiteness, take fig 46 of GMN5 (spectral networks), then
	p1 is the bottom steet, p2 will be the one on the top-left, 
	and p3 on the top-right. 
	The relavant relation here is (eq. A.3)
		" tau1 = nu6 nu2 "
	In terms of dashes, this means that the two incoming dashes 
	(supported on street p1) will be grown along different streets: 
	one called d_21, flowing from p2 into p1 (on sheet j in the picture), 
	and another called d_13 flowing from p1 into p3 
	(on sheet k in the picture)	(also, note the ordering!).
	Then a third dash d_32 will be added to the soliton, 
	flowing from p3 into p2 (on sheet i of the picture). 
	The new growing pairs will then be: 
	[..., [end(d_32), start(d_21)], [end(d_13), start(d_32)] ...]
	while the previous growing pair corresponding to the joint will
	of course be removed.
	"""
	new_solitons = []
	for p in old_cluster:
		new_soliton = copy_of_soliton(old_soliton)
		# Now for the growing pair p of DashEndpoints we 
		# identify the new corresponding growing pair in 
		# the new soliton
		old_pair_index = old_soliton.growing_pairs.index(p)
		new_p = new_soliton.growing_pairs[old_pair_index]

		# The joint in question. 
		# NOTE: must tart using new_p from here on!
		joint = new_p[0].end_point
		# The slot of street p_1 at the joint in question
		s_1 = new_p[0].slot
		# the slots of streets 2 and 3 
		# (recall that ordering of slots is CCW, as in GMN)
		# so hwre we first have s_1, then an empty slot, then s_3, 
		# then empty, then s_2.
		if s_1 + 4 < 6:
			s_2 = s_1 + 4
		else:
			s_2 = s_1 - 2
		if s_1 + 2 < 6:
			s_3 = s_1 + 2
		else:
			s_3 = s_1 - 4

		# the three streets
		p1 = joint.streets[s_1]
		p2 = joint.streets[s_2]
		p3 = joint.streets[s_3]

		# now identify d_21 and d_13, and grow them.
		# it will be important to check the orientation
		# of each growth point in the pair
		if new_p[0].orientation=='in' and new_p[1].orientation=='out':
			d_13 = new_p[0].dash
			d_21 = new_p[1].dash
		elif new_p[0].orientation=='out' and new_p[1].orientation=='in':
			d_13 = new_p[1].dash
			d_21 = new_p[0].dash
		else:
			raise ValueError

		# check that all three slots are actually 
		# the ones available at the joint
		if (
			s_1 in joint.available_slots and 
			s_2 in joint.available_slots and
			s_3 in joint.available_slots
		):
			pass
		else:
			raise ValueError
		
		# Now grow dashes d_13 and d_21
		# note that d_13 must be grown forward, so 
		# we specify the 'last' point for growth
		d_13.extend_dash_along_street(street=p3, end_pt='last', slot=s_3)
		# instead d_21 must be grown backward, so 
		# we specify the 'first' point for growth
		d_21.extend_dash_along_street(street=p2, end_pt='first', slot=s_2)
		

		# Then, we create a new dash, 
		d_32 = Dash(
			label='joint_new_dash_'+joint.label, 
			growth_restriction=None
		)
		# and extend it first along p2
		# NOTE: this will fix the overall orientation of the dash!
		d_32.extend_dash_along_street(
			street=p2, end_pt=joint, slot=s_2
		)
		# then also extend along p3, by growing backwards 
		# (hence specify the 'first' point for growth)
		d_32.extend_dash_along_street(
			street=p3, end_pt='first', slot=s_3
		)

		# At this point we add the new dash to the soliton's dashes
		# it must be inserted just after d_13 and before d_21
		indx = new_soliton.dashes.index(d_21)
		new_soliton.dashes.insert(indx, d_32)

		# And the update the soliton's growing pairs
		new_soliton.growing_pairs = (
			growing_pairs_from_dashes(new_soliton.dashes)
		)

		new_solitons.append(new_soliton)
		
	return new_solitons


def j_type_4_A_growth(old_soliton, old_cluster):
	"""
	Joints of type 4A are those with four two-way streets p1, p2, p3, p4
	ending on them. 

		p4		p3 		 p2
		  `      |      '
		   `     |     '
		    `    |    '
		     `   |   '
		      `  |  '
		       ` | '
		         |
		         |
		         |
		         |
		         |
		         |
		         p1

	The relavant relations here are 
	(from eq. A.3 of 1204.4824, but with our labels)
		tau1 = nu3 + nu2 nu4
		tau2 = nu4 nu1
		tau3 = nu1
		tau4 = nu1 nu2
	
	In terms of dashes, this means that there will be a complex behavior,
	involving both the one seen in type_3 joints, 
	and simple straight propagation.
	"""

	# first of all, determine on which street 
	# the soliton is approaching the joint

	new_solitons = []
	joint = old_cluster[0][0].end_point
	sol_slot = old_cluster[0][0].slot
	av_slots = joint.available_slots

	previous_slot = (sol_slot - 1) % 6
	next_slot = (sol_slot + 1) % 6

	incoming_street = None

	if (
		previous_slot not in av_slots and
		next_slot not in av_slots
	):
		incoming_street = 'p1'
	elif (
		previous_slot not in av_slots and
		next_slot in av_slots
	):
		incoming_street = 'p2'
	elif (
		previous_slot in av_slots and
		next_slot in av_slots
	):
		incoming_street = 'p3'
	elif (
		previous_slot in av_slots and
		next_slot not in av_slots
	):
		incoming_street = 'p4'
	else:
		raise Exception('Cannot determine position of soliton on joint.')

	# Now handle various cases separately.

	if incoming_street=='p1':
		# In this case there is both straight propagation, and Y-propagation

		# Start with straight propagation
		new_slot = (sol_slot + 3) % 6
		new_street = joint.streets[new_slot]
		for p in old_cluster:
			new_soliton_straight = copy_of_soliton(old_soliton)
			# Now for the growing pair p of DashEndpoints we 
			# identify the new corresponding growing pair in 
			# the new soliton
			old_pair_index = old_soliton.growing_pairs.index(p)
			new_p = new_soliton_straight.growing_pairs[old_pair_index]

			# now identify the dashes d_in and d_out, and grow them.
			# it will be important to check the orientation
			# of each growth point in the pair
			if new_p[0].orientation=='in' and new_p[1].orientation=='out':
				d_in = new_p[0].dash
				d_out = new_p[1].dash
			elif new_p[0].orientation=='out' and new_p[1].orientation=='in':
				d_in = new_p[1].dash
				d_out = new_p[0].dash
			else:
				raise ValueError

			# Now grow dashes d_in and d_out
			# note that d_in must be grown forward, so 
			# we specify the 'last' point for growth
			d_in.extend_dash_along_street(
				street=new_street, end_pt='last', slot=new_slot
			)

			# instead d_out must be grown backward, so 
			# we specify the 'first' point for growth
			d_out.extend_dash_along_street(
				street=new_street, end_pt='first', slot=new_slot
			)

			# Then update the soliton's growing pairs
			new_soliton_straight.growing_pairs = (
				growing_pairs_from_dashes(new_soliton_straight.dashes)
			)

			new_solitons.append(new_soliton_straight)

		# Then proceed with Y-propagation (like in type_3 joint)
		# in this section we use the labeling conventions explained
		# in the function 'j_type_3_growth'
		for p in old_cluster:
			new_soliton_Y = copy_of_soliton(old_soliton)
			# Now for the growing pair p of DashEndpoints we 
			# identify the new corresponding growing pair in 
			# the new soliton
			old_pair_index = old_soliton.growing_pairs.index(p)
			new_p = new_soliton_Y.growing_pairs[old_pair_index]

			# The joint in question. 
			# NOTE: must tart using new_p from here on!
			# The slot of street p_1 at the joint in question
			s_1 = new_p[0].slot
			# the slots of streets 2 and 3 
			# (recall that ordering of slots is CCW, as in GMN)
			# so hwre we first have s_1, then an empty slot, then s_3, 
			# then empty, then s_2.
			if s_1 + 4 < 6:
				s_2 = s_1 + 4
			else:
				s_2 = s_1 - 2
			if s_1 + 2 < 6:
				s_3 = s_1 + 2
			else:
				s_3 = s_1 - 4

			# the three streets
			p1 = joint.streets[s_1]
			p2 = joint.streets[s_2]
			p3 = joint.streets[s_3]

			# now identify d_21 and d_13, and grow them.
			# it will be important to check the orientation
			# of each growth point in the pair
			if new_p[0].orientation=='in' and new_p[1].orientation=='out':
				d_13 = new_p[0].dash
				d_21 = new_p[1].dash
			elif new_p[0].orientation=='out' and new_p[1].orientation=='in':
				d_13 = new_p[1].dash
				d_21 = new_p[0].dash
			else:
				raise ValueError

			# check that all three slots are actually 
			# the ones available at the joint
			if (
				s_1 in joint.available_slots and 
				s_2 in joint.available_slots and
				s_3 in joint.available_slots
			):
				pass
			else:
				raise ValueError
			
			# Now grow dashes d_13 and d_21
			# note that d_13 must be grown forward, so 
			# we specify the 'last' point for growth
			d_13.print_path_info()
			d_13.extend_dash_along_street(street=p3, end_pt='last', slot=s_3)
			# instead d_21 must be grown backward, so 
			# we specify the 'first' point for growth
			d_21.extend_dash_along_street(street=p2, end_pt='first', slot=s_2)
			
			# Then, we create a new dash, 
			d_32 = Dash(
				label='joint_new_dash_'+joint.label, 
				growth_restriction=None
			)
			# and extend it first along p2
			# NOTE: this will fix the overall orientation of the dash!
			d_32.extend_dash_along_street(
				street=p2, end_pt=joint, slot=s_2
			)
			# then also extend along p3, by growing backwards 
			# (hence specify the 'first' point for growth)
			d_32.extend_dash_along_street(
				street=p3, end_pt='first', slot=s_3
			)

			# At this point we add the new dash to the soliton's dashes
			# it must be inserted just after d_13 and before d_21
			indx = new_soliton_Y.dashes.index(d_21)
			new_soliton_Y.dashes.insert(indx, d_32)

			# And the update the soliton's growing pairs
			new_soliton_Y.growing_pairs = (
				growing_pairs_from_dashes(new_soliton_Y.dashes)
			)

			new_solitons.append(new_soliton_Y)

		return new_solitons
	
	else:
		return NotImplementedError






