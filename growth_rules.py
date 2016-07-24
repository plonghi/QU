from solitons import (
	copy_of_soliton, join_dashes_at_branch_point, SolitonPath, Dash, 
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

		for cl in original_clusters:
			# update the solitons in the list old_solitons, for each growing
			# cluster, returning a new list called 'new_solitons' which
			# will supersed old_solitons, thus being updated at the
			# next growth step (next cluster)

			node = cl[0][0].end_point
			node_type = node.type

			# if node_type == 'type_1_branch_point':
			# 	grow_at_node = bp_type_1_growth
			# elif node_type == 'type_3_joint':
			# 	grow_at_node = j_type_3_growth
			# elif node_type == 'type_2_branch_point':
			# 	grow_at_node = bp_type_2_growth
			# elif node_type == 'type_3_branch_point':
			# 	grow_at_node = bp_type_3_growth
			# elif node_type == 'type_4_A_joint':
			# 	grow_at_node = j_type_4_A_growth
			# elif node_type == 'type_4_B_joint':
			# 	grow_at_node = j_type_4_B_growth
			# else:
			# 	raise NotImplementedError

			if node_type == 'type_1_branch_point':
				grow_at_node = bp_type_1_growth
			elif node_type == 'type_2_branch_point':
				grow_at_node = bp_type_2_growth
			elif node_type == 'type_3_branch_point':
				grow_at_node = bp_type_3_growth
			else:
				grow_at_node = j_type_six_way

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


### TODO: merge the handling of the three branch-point types.

def bp_type_1_growth(old_soliton, old_cluster):
	"""
	Branch points of type 1 are those with a single two-way street 
	ending on them, so all that can happen is that each pair of
	dashes "caps off" there.
	There can be no mixing between different pairs, this follows
	from the algebra of tau and nu generating functions of 2d solitons.
	"""
	new_solitons = soliton_capping_off(old_soliton, old_cluster)
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
	new_solitons += soliton_capping_off(old_soliton, old_cluster)
	
	# Then, if pertinent, also handle the growth of the soliton through
	# the branch point.
	# It depends on the branch point's slot after the soliton's occupied 
	# slot, counter-clockwise. If this slot is also attached to a street 
	# (i.e. if it appears in available_slots), then the soliton 
	# can propagate as explained in point 2. Otherwise it cannot.
	next_slot = (sol_slot + 1) % 3
	if next_slot in av_slots:
		next_street = branch_pt.streets[next_slot]
		# We use 'straight propagation'
		new_solitons += soliton_propagation_1(
			old_soliton, old_cluster, next_street, next_slot
		)		
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
	new_solitons += soliton_capping_off(old_soliton, old_cluster)

	# Then also handle the growth of the soliton through
	# the branch point.
	next_slot = (sol_slot + 1) % 3
	next_street = branch_pt.streets[next_slot]

	# We use 'straight propagation'
	new_solitons += soliton_propagation_1(
		old_soliton, old_cluster, next_street, next_slot
	)		
		
	return new_solitons


### The single joint types have been deprecated.
### now they are handled by the 6-way joint.

# def j_type_3_growth(old_soliton, old_cluster):
# 	"""
# 	Joints of type 3 are those with a three two-way streets p1, p2, p3
# 	ending on them, so when a soliton ends on them from one of 
# 	the streets (say p1), it will keep growing on both the other two, p2, p3.
# 	For definiteness, take fig 46 of GMN5 (spectral networks), then
# 	p1 is the bottom steet, p2 will be the one on the top-left, 
# 	and p3 on the top-right. 
# 	The relavant relation here is (eq. A.3)
# 		" tau1 = nu6 nu2 "
# 	In terms of dashes, this means that the two incoming dashes 
# 	(supported on street p1) will be grown along different streets: 
# 	one called d_21, flowing from p2 into p1 (on sheet j in the picture), 
# 	and another called d_13 flowing from p1 into p3 
# 	(on sheet k in the picture)	(also, note the ordering!).
# 	Then a third dash d_32 will be added to the soliton, 
# 	flowing from p3 into p2 (on sheet i of the picture). 
# 	The new growing pairs will then be: 
# 	[..., [end(d_32), start(d_21)], [end(d_13), start(d_32)] ...]
# 	while the previous growing pair corresponding to the joint will
# 	of course be removed.

# 	Now all of this is simply handled by the function 
# 	'soliton_propagation_2'
# 	"""
# 	new_solitons = soliton_propagation_2(old_soliton, old_cluster)
# 	return new_solitons


# def j_type_4_A_growth(old_soliton, old_cluster):
# 	"""
# 	Joints of type 4A are those with four two-way streets p1, p2, p3, p4
# 	ending on them, formin a 'peace-sign' 

# 		p4		p3 		 p2
# 		  `      |      '
# 		   `     |     '
# 		    `    |    '
# 		     `   |   '
# 		      `  |  '
# 		       ` | '
# 		         |
# 		         |
# 		         |
# 		         |
# 		         |
# 		         |
# 		         p1

# 	The relavant relations here are 
# 	(from eq. A.3 of 1204.4824, but with our labels)
# 		tau1 = nu3 + nu2 nu4
# 		tau2 = nu4 nu1
# 		tau3 = nu1
# 		tau4 = nu1 nu2
	
# 	In terms of dashes, this means that there will be a complex behavior,
# 	involving both the one seen in type_3 joints, 
# 	and simple straight propagation.
# 	"""

# 	# first of all, determine on which street 
# 	# the soliton is approaching the joint

# 	new_solitons = []
# 	joint = old_cluster[0][0].end_point
# 	sol_slot = old_cluster[0][0].slot
# 	av_slots = joint.available_slots

# 	previous_slot = (sol_slot - 1) % 6
# 	next_slot = (sol_slot + 1) % 6

# 	incoming_street = None

# 	if (
# 		previous_slot not in av_slots and
# 		next_slot not in av_slots
# 	):
# 		incoming_street = 'p1'
# 	elif (
# 		previous_slot not in av_slots and
# 		next_slot in av_slots
# 	):
# 		incoming_street = 'p2'
# 	elif (
# 		previous_slot in av_slots and
# 		next_slot in av_slots
# 	):
# 		incoming_street = 'p3'
# 	elif (
# 		previous_slot in av_slots and
# 		next_slot not in av_slots
# 	):
# 		incoming_street = 'p4'
# 	else:
# 		raise Exception('Cannot determine position of soliton on joint.')

# 	# Now handle various cases separately.

# 	if incoming_street=='p1':
# 		# In this case there is both straight propagation, and Y-propagation

# 		# Start with straight propagation
# 		new_slot = (sol_slot + 3) % 6
# 		new_street = joint.streets[new_slot]
# 		new_solitons = soliton_propagation_1(
# 			old_soliton, old_cluster, new_street, new_slot
# 		)

# 		# Then proceed with Y-propagation (like in type_3 joint)
# 		new_solitons += soliton_propagation_2(old_soliton, old_cluster)
# 		return new_solitons

# 	elif incoming_street=='p2' or incoming_street=='p4':
# 		# In these cases there is only Y-propagation
# 		new_solitons = soliton_propagation_2(old_soliton, old_cluster)
# 		return new_solitons

# 	elif incoming_street=='p3':
# 		# In these cases there is only straight propagation
# 		new_slot = (sol_slot + 3) % 6
# 		new_street = joint.streets[new_slot]
# 		new_solitons = soliton_propagation_1(
# 			old_soliton, old_cluster, new_street, new_slot
# 		)
# 		return new_solitons
	
# 	else:
# 		return NotImplementedError


# def j_type_4_B_growth(old_soliton, old_cluster):
# 	"""
# 	Joints of type 4B are those with four two-way streets p1, p2, p3, p4
# 	ending on them, forming a cross

# 		p4		 		 p3
# 		  `             '
# 		   `           '
# 		    `         '
# 		     `       '
# 		      `     '
# 		       `   '
# 		        ` '
# 		        ' `
# 		       '   `
# 		      '     `
# 		     '       `
# 		    '         `
# 		   '           `
# 		 p1             p2

# 	The relavant relations here are 
# 	(from eq. A.3 of 1204.4824, but with our labels)
# 		tau1 = nu3
# 		tau2 = nu4 + nu3 nu1 nu4
# 		tau3 = nu1
# 		tau4 = nu2 + nu1 nu1 nu2
	
# 	In terms of dashes, this means that there will be a simple behavior
# 	for p1 and p3, involving just simple straight propagation.
# 	However for p2 and p4 there will be an additional contribution,
# 	coming from detours.
# 	"""

# 	# first of all, determine on which street 
# 	# the soliton is approaching the joint
# 	# due to the symmetry of the picture, there 
# 	# are really just two cases: p1 or p2.

# 	new_solitons = []
# 	joint = old_cluster[0][0].end_point
# 	sol_slot = old_cluster[0][0].slot
# 	av_slots = joint.available_slots

# 	previous_slot = (sol_slot - 1) % 6
# 	next_slot = (sol_slot + 1) % 6

# 	incoming_street = None

# 	if (
# 		previous_slot not in av_slots and
# 		next_slot in av_slots
# 	):
# 		incoming_street = 'p1'
# 	elif (
# 		previous_slot in av_slots and
# 		next_slot not in av_slots
# 	):
# 		incoming_street = 'p2'
# 	else:
# 		raise Exception('Cannot determine position of soliton on joint.')

# 	# Now handle various cases separately.

# 	if incoming_street=='p1':
# 		new_slot = (sol_slot + 3) % 6
# 		new_street = joint.streets[new_slot]

# 		new_solitons += soliton_propagation_1(
# 			old_soliton, old_cluster, new_street, new_slot
# 		)

# 	elif incoming_street=='p2':
# 		# In this case there is both straight propagation, and Y-propagation

# 		# Start with straight propagation
# 		new_slot = (sol_slot + 3) % 6
# 		new_street = joint.streets[new_slot]
# 		new_solitons = soliton_propagation_1(
# 			old_soliton, old_cluster, new_street, new_slot
# 		)

# 		# Then proceed with Y-propagation (like in type_3 joint)
# 		new_solitons += soliton_propagation_3(
# 			old_soliton, old_cluster
# 		)

# 	return new_solitons


def j_type_six_way(old_soliton, old_cluster):
	"""
	This function handles soliton propagation on a generic joint,
	which could be of type 3, 4_A, 4_B, 5 or 6.

		p5		p4		 p3
		  `      |      '
		   `     |     '
		    `    |    '
		     `   |   '
		      `  |  '
		       ` | '
		        `|'
		        '|`
		       ' | `
		      '  |  `
		     '   |   `
		    '    |    `
		   '     |     `
		  '      |      `
		p6      p1       p2

	The relavant relations here are 
	(from eq. A.3 of 1204.4824, but with our labels)
		tau1 = nu4 
		     + nu3 nu5 
		     + nu3 nu6 nu4 
		     + nu3 nu5 nu1 nu4 
		     + nu3 nu5 nu2 nu6 nu4 
		     + nu3 nu5 nu1 nu3 nu6 nu4
		     + ...

	NOTE: the six-way joint is very tricky, it doesn't have a 
	Z_6 symmetry, but a semi-direct product Z_2 x Z_3 symmetry.
	To see how it works, let tau1 be of type kj, tau2 of type ij, etc
	as in GMN's figure 46. Then the above formula reflects correctly 
	the path composition. (note that our nu's differ from GMN's)
	tau1  >   nu3   nu6   nu4
 	(kj)      (ki)  (ik)  (kj)

	However, if tau1 were of type ki, and tau2 of type kj, tau3 of type ij,
	etc, i.e. identifying tau1 with GMN's tau6 and so on, the formula 
	will be very different.
	In this case the order of the factors of each summand must be REVERSED 
	because now we would have:
	tau1  >   nu3   nu6   nu4
 	(ki)      (ji)  (ij)  (ki)

 	This subtle symmetry was noted in GMN's paper, in the comment just 
 	below equation (A.2)

 	NOTE: When we build a network, we should probably be careful in specifying 
 	how joints exactly look, i.e. to which slot (ODD/EVEN) we assign 
 	each street.
 	TODO: STUDY THIS ISSUE FURTHER, GLOBAL CONSTRAINTS?!!!

 	For propagation at a joint, we will thus have the ORIGINAL equation 
 	working for streets at slots 0, 2, 4; while the REVERSED equations must
 	be employed for streets at slots 1, 3, 5.
	
	For EVEN-slot streets, the propagation will thus look like:
		> For the first piece there will be simple straight propagation.
		> For the second piece there will be Y-propagation
			- the incoming dash of p1 will grow along p3, call it d_13
			- a new dash from p3 to p5, call it d_35
			- the outgoing dash of p1 will grow (backwards) along p5, call it d_51
			then the contribution to the soliton will be to replace the 
			original dashes of the growing pair 
			[..., d1_in, d1_out ,...]
			with
			[..., d_13, d_35, d_51,...]
		> The third piece will involve four dashes now, 
			in order they are:
			- the incoming dash of p1 will grow along p3, call it d_13
			- a new dash from p3 to p6, call it d_36
			- a new dash from p6 to p4, call it d_64
			- the outgoing dash of p1 will grow (backwards) along p4, 
			call it d_41 then the contribution to the soliton will 
			be to replace the original dashes of the growing pair 
			[..., d1_in, d1_out ,...]
			with
			[..., d_13, d_36, d_64, d_41,...]

	While for ODD-slot streets, the propagation will thus look like:
		> For the first piece there will be simple straight propagation.
		> For the second piece there will be Y-propagation
			- the incoming dash of p1 will grow along p5, call it d_15
			- a new dash from p5 to p3, call it d_53
			- the outgoing dash of p1 will grow (backwards) along p3, 
			call it d_31 then the contribution to the soliton will 
			be to replace the original dashes of the growing pair 
			[..., d1_in, d1_out ,...]
			with
			[..., d_15, d_53, d_31,...]
		> The third piece will involve four dashes now, 
			in order they are:
			- the incoming dash of p1 will grow along p4, call it d_14
			- a new dash from p4 to p6, call it d_46
			- a new dash from p6 to p3, call it d_63
			- the outgoing dash of p1 will grow (backwards) along p3, call it d_31
			then the contribution to the soliton will be to replace the 
			original dashes of the growing pair 
			[..., d1_in, d1_out ,...]
			with
			[..., d_14, d_46, d_63, d_31,...]
	"""	
	new_solitons = []
	joint = old_cluster[0][0].end_point
	sol_slot = old_cluster[0][0].slot
	if sol_slot % 2 == 0:
		parity_reversal = False
	else:
		parity_reversal = True
	av_slots = joint.available_slots
	# let p1 be the street from which 
	# the soliton is approaching the joint.

	# the slot of each street
	[s_1, s_2, s_3, s_4, s_5, s_6] = [(sol_slot + i) % 6 for i in range(6)]
	s_k = [s_1, s_2, s_3, s_4, s_5, s_6] 
	
	# the streets
	[p1, p2, p3, p4, p5, p6] = [joint.streets[s_k[i]] for i in range(6)]
	
	# then map which streets are 'available'
	pk_is_av = [False for i in range(6)]
	for k in range(6):
		if s_k[k] in av_slots:
			pk_is_av[k] = True 

	[p1_is_av, p2_is_av, p3_is_av, p4_is_av, p5_is_av, p6_is_av] = pk_is_av

	# # tau1 > nu4
	# if (
	# 	p4_is_av is True
	# ):
	# 	# Just straight propagation
	# 	new_solitons += soliton_propagation_1(
	# 		old_soliton, old_cluster, p4, s_4,
	# 	)

	# # tau1 > nu3 nu5
	# if (
	# 	p3_is_av is True and
	# 	p5_is_av is True
	# ):
	# 	# Y-propagation
	# 	new_solitons += soliton_propagation_2(
	# 		old_soliton, old_cluster, parity,
	# 	)

	# # tau1 > nu3 nu6 nu4 
	# if (
	# 	p3_is_av is True and
	# 	p6_is_av is True and
	# 	p4_is_av is True
	# ):
	# 	# X-propagation
	# 	new_solitons += soliton_propagation_3(
	# 		old_soliton, old_cluster, parity,
	# 	)

	# # tau1 > nu3 nu5 nu1 nu4 
	# if (
	# 	p3_is_av is True and
	# 	p5_is_av is True and
	# 	p1_is_av is True and
	# 	p4_is_av is True
	# ):
	# 	# type-4 propagation
	# 	new_solitons += soliton_propagation_4(
	# 		old_soliton, old_cluster, parity,
	# 	)

	# # tau1 > nu3 nu5 nu2 nu6 nu4 
	# if (
	# 	p3_is_av is True and
	# 	p5_is_av is True and
	# 	p2_is_av is True and
	# 	p6_is_av is True and
	# 	p4_is_av is True
	# ):
	# 	# type-5 propagation
	# 	new_solitons += soliton_propagation_5(
	# 		old_soliton, old_cluster, parity,
	# 	)

	# # tau1 > nu3 nu5 nu1 nu3 nu6 nu4
	# if (
	# 	p3_is_av is True and
	# 	p5_is_av is True and
	# 	p1_is_av is True and
	# 	p6_is_av is True and
	# 	p4_is_av is True
	# ):
	# 	# type-6 propagation
	# 	new_solitons += soliton_propagation_6(
	# 		old_soliton, old_cluster, parity,
	# 	)
	
	# ### TODO: if all slots are available, add more iterations!

	# tau1 > nu4
	if (
		p4_is_av is True
	):
		# Just straight propagation, along street at the opposite slot (+3)
		new_solitons += soliton_propagation_street_sequence(
			old_soliton, old_cluster, [3], 
			parity_reversed=parity_reversal
		)

	# tau1 > nu3 nu5
	if (
		p3_is_av is True and
		p5_is_av is True
	):
		# Y-propagation on streets 3 (+2) then 5 (+4)
		new_solitons += soliton_propagation_street_sequence(
			old_soliton, old_cluster, [2, 4], 
			parity_reversed=parity_reversal
		)

	# tau1 > nu3 nu6 nu4 
	if (
		p3_is_av is True and
		p6_is_av is True and
		p4_is_av is True
	):
		# X-propagation on streets 3 (+2), 6 (+5), 4 (+3)
		new_solitons += soliton_propagation_street_sequence(
			old_soliton, old_cluster, [2, 5, 3], 
			parity_reversed=parity_reversal
		)

	# tau1 > nu3 nu5 nu1 nu4 
	if (
		p3_is_av is True and
		p5_is_av is True and
		p1_is_av is True and
		p4_is_av is True
	):
		# propagation on streets 3 (+2), 5 (+4), 1 (+0), 4 (+3)
		new_solitons += soliton_propagation_street_sequence(
			old_soliton, old_cluster, [2, 4, 0, 3], 
			parity_reversed=parity_reversal
		)

	# tau1 > nu3 nu5 nu2 nu6 nu4 
	if (
		p3_is_av is True and
		p5_is_av is True and
		p2_is_av is True and
		p6_is_av is True and
		p4_is_av is True
	):
		# propagation on streets 3, 5, 2, 6, 4 i.e. relatively 
		# [+2, +4, +1, +5, +3]
		new_solitons += soliton_propagation_street_sequence(
			old_soliton, old_cluster, [2, 4, 1, 5, 3], 
			parity_reversed=parity_reversal
		)

	# tau1 > nu3 nu5 nu1 nu3 nu6 nu4
	if (
		p3_is_av is True and
		p5_is_av is True and
		p1_is_av is True and
		p6_is_av is True and
		p4_is_av is True
	):
		# propagation on streets 3, 5, 1, 3, 6, 4 i.e. relatively 
		# [+2, +4, +0, +2, +5, +3]
		new_solitons += soliton_propagation_street_sequence(
			old_soliton, old_cluster, [2, 4, 0, 2, 5, 3], 
			parity_reversed=parity_reversal
		)
	
	### TODO: if all slots are available, add more iterations!

	return new_solitons


def soliton_capping_off(old_soliton, old_cluster):
	"""
	Handles capping-off of a soliton at a branch point.
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
		new_dash = join_dashes_at_branch_point(new_p)

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


def create_dash_through_joint(
	joint, in_street, out_street, in_slot, out_slot
):
	"""
	Creates a dash stretching from 'in_street', through 'joint' and out on
	'out_street'.
	"""
	### TODO:
	### - do some checks on the consistency of streets being attached to joint
	# NOTE: it's necessary that slots are specified externally, 
	# they cannot be computed here. Because a street may end on a 
	# branch point or joint more than once (twice) leading to ambiguities.

	new_dash = Dash(
		label='joint_new_dash_'+joint.label, 
		growth_restriction=None
	)
	# and extend it first along the out_street
	# NOTE: this will fix the overall orientation of the dash!
	new_dash.extend_dash_along_street(
		street=out_street, end_pt=joint, slot=out_slot
	)
	# then also extend along in_street, by growing backwards 
	# (hence specify the 'first' point for growth)
	new_dash.extend_dash_along_street(
		street=in_street, end_pt='first', slot=in_slot
	)

	return new_dash


# def soliton_propagation_1(old_soliton, old_cluster, new_street, new_slot):
# 	"""
# 	Handles propagation of a soliton along a single new street,
# 	throught the specified slot (sometimes a street ends on a nodal
# 	point on more than one slot).
# 	"""
# 	new_solitons = []
# 	for p in old_cluster:
# 		new_soliton_straight = copy_of_soliton(old_soliton)
# 		# Now for the growing pair p of DashEndpoints we 
# 		# identify the new corresponding growing pair in 
# 		# the new soliton
# 		old_pair_index = old_soliton.growing_pairs.index(p)
# 		new_p = new_soliton_straight.growing_pairs[old_pair_index]

# 		# now identify the dashes d_in and d_out, and grow them.
# 		# it will be important to check the orientation
# 		# of each growth point in the pair
# 		if new_p[0].orientation=='in' and new_p[1].orientation=='out':
# 			d_in = new_p[0].dash
# 			d_out = new_p[1].dash
# 		elif new_p[0].orientation=='out' and new_p[1].orientation=='in':
# 			d_in = new_p[1].dash
# 			d_out = new_p[0].dash
# 		else:
# 			raise ValueError

# 		# Now grow dashes d_in and d_out
# 		# note that d_in must be grown forward, so 
# 		# we specify the 'last' point for growth
# 		d_in.extend_dash_along_street(
# 			street=new_street, end_pt='last', slot=new_slot
# 		)

# 		# instead d_out must be grown backward, so 
# 		# we specify the 'first' point for growth
# 		d_out.extend_dash_along_street(
# 			street=new_street, end_pt='first', slot=new_slot
# 		)

# 		# Then update the soliton's growing pairs
# 		new_soliton_straight.growing_pairs = (
# 			growing_pairs_from_dashes(new_soliton_straight.dashes)
# 		)

# 		new_solitons.append(new_soliton_straight)

# 	return new_solitons


# def soliton_propagation_2(old_soliton, old_cluster):
# 	"""
# 	Handles propagation of a soliton along a Y-type bifurcation, like for 
# 	joints of type 3.

# 	Joints of type 3 are those with a three two-way streets p1, p2, p3
# 	ending on them, so when a soliton ends on them from one of 
# 	the streets (say p1), it will keep growing on both the other two, p2, p3.
# 	For definiteness, take fig 46 of GMN5 (spectral networks), then
# 	p1 is the bottom steet, p2 will be the one on the top-left, 
# 	and p3 on the top-right. 
# 	The relavant relation here is (eq. A.3)
# 		" tau1 = nu6 nu2 "
# 	In terms of dashes, this means that the two incoming dashes 
# 	(supported on street p1) will be grown along different streets: 
# 	one called d_21, flowing from p2 into p1 (on sheet j in the picture), 
# 	and another called d_13 flowing from p1 into p3 
# 	(on sheet k in the picture)	(also, note the ordering!).
# 	Then a third dash d_32 will be added to the soliton, 
# 	flowing from p3 into p2 (on sheet i of the picture). 
# 	The new growing pairs will then be: 
# 	[..., [end(d_32), start(d_21)], [end(d_13), start(d_32)] ...]
# 	while the previous growing pair corresponding to the joint will
# 	of course be removed.
# 	"""
# 	new_solitons = []
# 	joint = old_cluster[0][0].end_point

# 	for p in old_cluster:
# 		new_soliton_Y = copy_of_soliton(old_soliton)
# 		# Now for the growing pair p of DashEndpoints we 
# 		# identify the new corresponding growing pair in 
# 		# the new soliton
# 		old_pair_index = old_soliton.growing_pairs.index(p)
# 		new_p = new_soliton_Y.growing_pairs[old_pair_index]

# 		# The joint in question. 
# 		# NOTE: must tart using new_p from here on!
# 		# The slot of street p_1 at the joint in question
# 		s_1 = new_p[0].slot
# 		# the slots of streets 2 and 3 
# 		# (recall that ordering of slots is CCW, as in GMN)
# 		# so hwre we first have s_1, then an empty slot, then s_3, 
# 		# then empty, then s_2.
# 		if s_1 + 4 < 6:
# 			s_2 = s_1 + 4
# 		else:
# 			s_2 = s_1 - 2
# 		if s_1 + 2 < 6:
# 			s_3 = s_1 + 2
# 		else:
# 			s_3 = s_1 - 4

# 		# the three streets
# 		p1 = joint.streets[s_1]
# 		p2 = joint.streets[s_2]
# 		p3 = joint.streets[s_3]

# 		# now identify d_21 and d_13, and grow them.
# 		# it will be important to check the orientation
# 		# of each growth point in the pair
# 		if new_p[0].orientation=='in' and new_p[1].orientation=='out':
# 			d_13 = new_p[0].dash
# 			d_21 = new_p[1].dash
# 		elif new_p[0].orientation=='out' and new_p[1].orientation=='in':
# 			d_13 = new_p[1].dash
# 			d_21 = new_p[0].dash
# 		else:
# 			raise ValueError

# 		# check that all three slots are actually 
# 		# the ones available at the joint
# 		if (
# 			s_1 in joint.available_slots and 
# 			s_2 in joint.available_slots and
# 			s_3 in joint.available_slots
# 		):
# 			pass
# 		else:
# 			raise ValueError
		
# 		# Now grow dashes d_13 and d_21
# 		# note that d_13 must be grown forward, so 
# 		# we specify the 'last' point for growth
# 		d_13.extend_dash_along_street(street=p3, end_pt='last', slot=s_3)
# 		# instead d_21 must be grown backward, so 
# 		# we specify the 'first' point for growth
# 		d_21.extend_dash_along_street(street=p2, end_pt='first', slot=s_2)
		
# 		# Then, we create a new dash, 
# 		d_32 = create_dash_through_joint(joint, p3, p2, s_3, s_2)

# 		# At this point we add the new dash to the soliton's dashes
# 		# it must be inserted just after d_13 and before d_21
# 		indx = new_soliton_Y.dashes.index(d_21)
# 		new_soliton_Y.dashes.insert(indx, d_32)

# 		# And the update the soliton's growing pairs
# 		new_soliton_Y.growing_pairs = (
# 			growing_pairs_from_dashes(new_soliton_Y.dashes)
# 		)

# 		new_solitons.append(new_soliton_Y)

# 	return new_solitons


# def soliton_propagation_3(old_soliton, old_cluster):
# 	"""
# 	For a joint of type 4B with four two-way streets p1, p2, p3, p4

# 		p4		 		 p3
# 		  `             '
# 		   `           '
# 		    `         '
# 		     `       '
# 		      `     '
# 		       `   '
# 		        ` '
# 		        ' `
# 		       '   `
# 		      '     `
# 		     '       `
# 		    '         `
# 		   '           `
# 		 p1             p2

# 	The relavant relations here are 
# 	(from eq. A.3 of 1204.4824, but with our labels)
# 		tau1 = nu3
# 		tau2 = nu4 + nu3 nu1 nu4
# 		tau3 = nu1
# 		tau4 = nu2 + nu1 nu1 nu2
	
# 	This propagation describes the joining of dashes for the piece
# 	tau2 > nu3 nu1 nu4

# 	There will be four dashes now involved:, in order they are:
# 	- the incoming dash of p2 will grow along p3, call it d_23
# 	- a new dash from p3 to p1, call it d_31
# 	- a new dash from p1 to p4, call it d_14
# 	- the outgoing dash of p2 will grow (backwards) along p4, call it d_42

# 	"""	
# 	new_solitons = []
# 	joint = old_cluster[0][0].end_point

# 	for p in old_cluster:
# 		new_soliton = copy_of_soliton(old_soliton)
# 		# Now for the growing pair p of DashEndpoints we 
# 		# identify the new corresponding growing pair in 
# 		# the new soliton
# 		old_pair_index = old_soliton.growing_pairs.index(p)
# 		new_p = new_soliton.growing_pairs[old_pair_index]

# 		# The joint in question. 
# 		# NOTE: must tart using new_p from here on!
# 		# The slot of street p_2 at the joint in question
# 		s_2 = new_p[0].slot
# 		# the slots of streets 1, 3 and 4 
# 		s_1 = (s_2 - 1) % 6
# 		s_3 = (s_2 + 2) % 6
# 		s_4 = (s_2 + 3) % 6

# 		# the four streets
# 		p1 = joint.streets[s_1]
# 		p2 = joint.streets[s_2]
# 		p3 = joint.streets[s_3]
# 		p4 = joint.streets[s_4]

# 		# now identify d_23 and d_42, and grow them.
# 		# It will be important to check the orientation
# 		# of each growth point in the pair
# 		if new_p[0].orientation=='in' and new_p[1].orientation=='out':
# 			d_23 = new_p[0].dash
# 			d_42 = new_p[1].dash
# 		elif new_p[0].orientation=='out' and new_p[1].orientation=='in':
# 			d_23 = new_p[1].dash
# 			d_42 = new_p[0].dash
# 		else:
# 			raise ValueError

# 		# check that all four slots are actually 
# 		# the ones available at the joint
# 		if (
# 			s_1 in joint.available_slots and 
# 			s_2 in joint.available_slots and
# 			s_3 in joint.available_slots and
# 			s_4 in joint.available_slots
# 		):
# 			pass
# 		else:
# 			raise ValueError
		
# 		# Now grow dashes d_23 and d_42
# 		# note that d_23 must be grown forward, so 
# 		# we specify the 'last' point for growth
# 		d_23.extend_dash_along_street(street=p3, end_pt='last', slot=s_3)
# 		# instead d_42 must be grown backward, so 
# 		# we specify the 'first' point for growth
# 		d_42.extend_dash_along_street(street=p4, end_pt='first', slot=s_4)
		
# 		# Then, we create a new dash, 
# 		d_31 = create_dash_through_joint(joint, p3, p1, s_3, s_1)

# 		# Then, we create a new dash, 
# 		d_14 = create_dash_through_joint(joint, p1, p4, s_1, s_4)

# 		# At this point we add the new dashes to the soliton's dashes
# 		# they must be inserted just after d_23 and before d_42
# 		# moreover d_31 comes first, while d_14 comes later
# 		indx = new_soliton.dashes.index(d_42)
# 		# at the beginnig it's [... d_23, d_42 ...]
# 		new_soliton.dashes.insert(indx, d_31)
# 		# now it's [... d_23, d_31, d_42 ...]
# 		new_soliton.dashes.insert(indx + 1, d_14)
# 		# now it's [... d_23, d_31, d_14, d_42 ...]

# 		# Finally update the soliton's growing pairs
# 		new_soliton.growing_pairs = (
# 			growing_pairs_from_dashes(new_soliton.dashes)
# 		)

# 		new_solitons.append(new_soliton)

# 	return new_solitons



# def soliton_propagation_4(old_soliton, old_cluster):
# 	"""
# 	The relavant relation here is 
# 	tau1 > nu3 nu5 nu1 nu4

# 	There will be five dashes now involved:, in order they are:
# 	- the incoming dash of p1 will grow along p3, call it d_13
# 	- a new dash from p3 to p5, call it d_35
# 	- a new dash from p5 to p1, call it d_51
# 	- a new dash from p1 to p4, call it d_14
# 	- the outgoing dash of p1 will grow (backwards) along p4, call it d_41

# 	"""	
# 	new_solitons = []
# 	joint = old_cluster[0][0].end_point

# 	for p in old_cluster:
# 		new_soliton = copy_of_soliton(old_soliton)
# 		# Now for the growing pair p of DashEndpoints we 
# 		# identify the new corresponding growing pair in 
# 		# the new soliton
# 		old_pair_index = old_soliton.growing_pairs.index(p)
# 		new_p = new_soliton.growing_pairs[old_pair_index]

# 		# The joint in question. 
# 		# NOTE: must tart using new_p from here on!
# 		# The slot of street p_1 at the joint in question

# 		# the slot of each street
# 		[s_1, s_2, s_3, s_4, s_5, s_6] = (
# 			[(new_p[0].slot + i) % 6 for i in range(6)]
# 		)
# 		s_k = [s_1, s_2, s_3, s_4, s_5, s_6] 
		
# 		# the streets
# 		[p1, p2, p3, p4, p5, p6] = [joint.streets[s_k[i]] for i in range(6)]

# 		# now identify d_13 and d_41, and grow them.
# 		# It will be important to check the orientation
# 		# of each growth point in the pair
# 		if new_p[0].orientation=='in' and new_p[1].orientation=='out':
# 			d_13 = new_p[0].dash
# 			d_41 = new_p[1].dash
# 		elif new_p[0].orientation=='out' and new_p[1].orientation=='in':
# 			d_13 = new_p[1].dash
# 			d_41 = new_p[0].dash
# 		else:
# 			raise ValueError
		
# 		# Now grow dashes d_13 and d_41
# 		# note that the former must be grown forward, so 
# 		# we specify the 'last' point for growth
# 		d_13.extend_dash_along_street(street=p3, end_pt='last', slot=s_3)
# 		# instead the latter must be grown backward, so 
# 		# we specify the 'first' point for growth
# 		d_41.extend_dash_along_street(street=p4, end_pt='first', slot=s_4)
		
# 		# Then, we create a new dash, 
# 		d_35 = create_dash_through_joint(joint, p3, p5, s_3, s_5)

# 		# Then, we create a new dash, 
# 		d_51 = create_dash_through_joint(joint, p5, p1, s_5, s_1)

# 		# Then, we create a new dash, 
# 		d_14 = create_dash_through_joint(joint, p1, p4, s_1, s_4)

# 		# At this point we add the new dashes to the soliton's dashes
# 		# they must be inserted just after d_13 and before d_41
# 		# in the obvious order: d_ij, d_jk, d_kl, etc.
# 		indx = new_soliton.dashes.index(d_13) + 1
# 		# start adding
# 		new_soliton.dashes.insert(indx, d_35)
# 		new_soliton.dashes.insert(indx + 1, d_51)
# 		new_soliton.dashes.insert(indx + 2, d_14)

# 		# Finally update the soliton's growing pairs
# 		new_soliton.growing_pairs = (
# 			growing_pairs_from_dashes(new_soliton.dashes)
# 		)

# 		new_solitons.append(new_soliton)

# 	return new_solitons


# def soliton_propagation_5(old_soliton, old_cluster):
# 	"""
# 	The relavant relation here is 
# 	tau1 > nu3 nu5 nu2 nu6 nu4 

# 	There will be six dashes now involved:, in order they are:
# 	- the incoming dash of p1 will grow along p3, call it d_13
# 	- a new dash from p3 to p5, call it d_35
# 	- a new dash from p5 to p2, call it d_52
# 	- a new dash from p2 to p6, call it d_26
# 	- a new dash from p6 to p4, call it d_64
# 	- the outgoing dash of p1 will grow (backwards) along p4, call it d_41

# 	"""	
# 	new_solitons = []
# 	joint = old_cluster[0][0].end_point

# 	for p in old_cluster:
# 		new_soliton = copy_of_soliton(old_soliton)
# 		# Now for the growing pair p of DashEndpoints we 
# 		# identify the new corresponding growing pair in 
# 		# the new soliton
# 		old_pair_index = old_soliton.growing_pairs.index(p)
# 		new_p = new_soliton.growing_pairs[old_pair_index]

# 		# The joint in question. 
# 		# NOTE: must tart using new_p from here on!
# 		# The slot of street p_1 at the joint in question

# 		# the slot of each street
# 		[s_1, s_2, s_3, s_4, s_5, s_6] = (
# 			[(new_p[0].slot + i) % 6 for i in range(6)]
# 		)
# 		s_k = [s_1, s_2, s_3, s_4, s_5, s_6] 
		
# 		# the streets
# 		[p1, p2, p3, p4, p5, p6] = [joint.streets[s_k[i]] for i in range(6)]

# 		# now identify d_13 and d_41, and grow them.
# 		# It will be important to check the orientation
# 		# of each growth point in the pair
# 		if new_p[0].orientation=='in' and new_p[1].orientation=='out':
# 			d_13 = new_p[0].dash
# 			d_41 = new_p[1].dash
# 		elif new_p[0].orientation=='out' and new_p[1].orientation=='in':
# 			d_13 = new_p[1].dash
# 			d_41 = new_p[0].dash
# 		else:
# 			raise ValueError
		
# 		# Now grow dashes d_13 and d_41
# 		# note that the former must be grown forward, so 
# 		# we specify the 'last' point for growth
# 		d_13.extend_dash_along_street(street=p3, end_pt='last', slot=s_3)
# 		# instead the latter must be grown backward, so 
# 		# we specify the 'first' point for growth
# 		d_41.extend_dash_along_street(street=p4, end_pt='first', slot=s_4)
		
# 		# Then, we create a new dash, 
# 		d_35 = create_dash_through_joint(joint, p3, p5, s_3, s_5)

# 		# Then, we create a new dash, 
# 		d_52 = create_dash_through_joint(joint, p5, p2, s_5, s_2)

# 		# Then, we create a new dash, 
# 		d_26 = create_dash_through_joint(joint, p2, p6, s_2, s_6)

# 		# Then, we create a new dash, 
# 		d_64 = create_dash_through_joint(joint, p6, p4, s_6, s_4)

# 		# At this point we add the new dashes to the soliton's dashes
# 		# they must be inserted just after d_13 and before d_41
# 		# in the obvious order: d_ij, d_jk, d_kl, etc.
# 		indx = new_soliton.dashes.index(d_13) + 1
# 		# start adding
# 		new_soliton.dashes.insert(indx, d_35)
# 		new_soliton.dashes.insert(indx + 1, d_52)
# 		new_soliton.dashes.insert(indx + 2, d_26)
# 		new_soliton.dashes.insert(indx + 3, d_64)

# 		# Finally update the soliton's growing pairs
# 		new_soliton.growing_pairs = (
# 			growing_pairs_from_dashes(new_soliton.dashes)
# 		)

# 		new_solitons.append(new_soliton)

# 	return new_solitons


# def soliton_propagation_6(old_soliton, old_cluster):
# 	"""
# 	The relavant relation here is 
# 	tau1 > nu3 nu5 nu1 nu3 nu6 nu4

# 	There will be seven dashes now involved:, in order they are:
# 	- the incoming dash of p1 will grow along p3, call it d_13
# 	- a new dash from p3 to p5, call it d_35
# 	- a new dash from p5 to p1, call it d_51
# 	- a new dash from p1 to p3, call it d_13_a
# 	- a new dash from p3 to p6, call it d_36
# 	- a new dash from p6 to p4, call it d_64
# 	- the outgoing dash of p1 will grow (backwards) along p4, call it d_41

# 	"""	
# 	new_solitons = []
# 	joint = old_cluster[0][0].end_point

# 	for p in old_cluster:
# 		new_soliton = copy_of_soliton(old_soliton)
# 		# Now for the growing pair p of DashEndpoints we 
# 		# identify the new corresponding growing pair in 
# 		# the new soliton
# 		old_pair_index = old_soliton.growing_pairs.index(p)
# 		new_p = new_soliton.growing_pairs[old_pair_index]

# 		# The joint in question. 
# 		# NOTE: must tart using new_p from here on!
# 		# The slot of street p_1 at the joint in question

# 		# the slot of each street
# 		[s_1, s_2, s_3, s_4, s_5, s_6] = (
# 			[(new_p[0].slot + i) % 6 for i in range(6)]
# 		)
# 		s_k = [s_1, s_2, s_3, s_4, s_5, s_6] 
		
# 		# the streets
# 		[p1, p2, p3, p4, p5, p6] = [joint.streets[s_k[i]] for i in range(6)]

# 		# now identify d_13 and d_41, and grow them.
# 		# It will be important to check the orientation
# 		# of each growth point in the pair
# 		if new_p[0].orientation=='in' and new_p[1].orientation=='out':
# 			d_13 = new_p[0].dash
# 			d_41 = new_p[1].dash
# 		elif new_p[0].orientation=='out' and new_p[1].orientation=='in':
# 			d_13 = new_p[1].dash
# 			d_41 = new_p[0].dash
# 		else:
# 			raise ValueError
		
# 		# Now grow dashes d_13 and d_41
# 		# note that the former must be grown forward, so 
# 		# we specify the 'last' point for growth
# 		d_13.extend_dash_along_street(street=p3, end_pt='last', slot=s_3)
# 		# instead the latter must be grown backward, so 
# 		# we specify the 'first' point for growth
# 		d_41.extend_dash_along_street(street=p4, end_pt='first', slot=s_4)

# 		# Then, we create a new dash, 
# 		d_35 = create_dash_through_joint(joint, p3, p5, s_3, s_5)

# 		# Then, we create a new dash, 
# 		d_51 = create_dash_through_joint(joint, p5, p1, s_5, s_1)

# 		# Then, we create a new dash, 
# 		d_13_a = create_dash_through_joint(joint, p1, p3, s_1, s_3)
		
# 		# Then, we create a new dash, 
# 		d_36 = create_dash_through_joint(joint, p3, p6, s_3, s_6)

# 		# Then, we create a new dash, 
# 		d_64 = create_dash_through_joint(joint, p6, p4, s_6, s_4)

# 		# At this point we add the new dashes to the soliton's dashes
# 		# they must be inserted just after d_13 and before d_41
# 		# in the obvious order: d_ij, d_jk, d_kl, etc.
# 		indx = new_soliton.dashes.index(d_13) + 1
# 		# start adding
# 		new_soliton.dashes.insert(indx, d_35)
# 		new_soliton.dashes.insert(indx + 1, d_51)
# 		new_soliton.dashes.insert(indx + 2, d_13_a)
# 		new_soliton.dashes.insert(indx + 3, d_36)
# 		new_soliton.dashes.insert(indx + 4, d_64)

# 		# Finally update the soliton's growing pairs
# 		new_soliton.growing_pairs = (
# 			growing_pairs_from_dashes(new_soliton.dashes)
# 		)

# 		new_solitons.append(new_soliton)

# 	return new_solitons


def soliton_propagation_street_sequence(
	old_soliton, old_cluster, rel_street_sequence, parity_reversed=False
):
	"""
	Handles the propagaiton of a soliton at a joint.
	The relative street sequence should be a list like
	[2, 4, 0, 2, 5, 3]
	standing for a soliton arriving on street1 and propagating 
	on street 1+k for k an element of the list, i.e. for
	tau1 > nu3 nu5 nu1 nu3 nu6 nu4

	There will be a number of dashes involved.
	In this particular example they are:
	- the incoming dash of p1 will grow along p3, call it d_13
	- a new dash from p3 to p5, call it d_35
	- a new dash from p5 to p1, call it d_51
	- a new dash from p1 to p3, call it d_13_a
	- a new dash from p3 to p6, call it d_36
	- a new dash from p6 to p4, call it d_64
	- the outgoing dash of p1 will grow (backwards) along p4, call it d_41

	if the 'reversed' option is True, then propagation will be 
	handled backwards, i.e. starting from p1 then on p4, then p6, etc

	"""	
	if parity_reversed is True:
		rev_sequence = list(reversed(rel_street_sequence))
		rel_street_sequence = rev_sequence

	new_solitons = []
	joint = old_cluster[0][0].end_point

	for p in old_cluster:
		new_soliton = copy_of_soliton(old_soliton)
		# Now for the growing pair p of DashEndpoints we 
		# identify the new corresponding growing pair in 
		# the new soliton
		old_pair_index = old_soliton.growing_pairs.index(p)
		new_p = new_soliton.growing_pairs[old_pair_index]

		# The joint in question. 
		# NOTE: must tart using new_p from here on!
		# The slot of street p_1 at the joint in question is s_1
		# but it's not necessarily 0, of course.

		# the slot of each street
		[s_1, s_2, s_3, s_4, s_5, s_6] = (
			[(new_p[0].slot + i) % 6 for i in range(6)]
		)
		s_k = [s_1, s_2, s_3, s_4, s_5, s_6] 
		
		# the streets
		[p1, p2, p3, p4, p5, p6] = [joint.streets[s_k[i]] for i in range(6)]
		p_k = [p1, p2, p3, p4, p5, p6]

		# the street sequence corresponding to the path
		street_sequence = [p_k[k] for k in rel_street_sequence]

		# the slot sequence corresponding to the path
		slot_sequence = [s_k[k] for k in rel_street_sequence]

		# now identify d_in (the analog of d_13 in the example) 
		# and d_out (the analog of d_41 in the example), and grow them.
		# It will be important to check the orientation
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
		# note that the former must be grown forward on the first 
		# street of the sequence, so we specify the 'last' point for growth	
		d_in.extend_dash_along_street(
			street=street_sequence[0], 
			end_pt='last', 
			slot=slot_sequence[0],
		)
		# instead the latter must be grown backward on the last
		# street, so we specify the 'first' point for growth
		d_out.extend_dash_along_street(
			street=street_sequence[-1], 
			end_pt='first', 
			slot=slot_sequence[-1],
		)

		# Then, for each pair of streets on which propagation takes place,
		# we create a new dash
		intermediate_dashes = []
		for i in range(len(street_sequence) - 1):
			p_a = street_sequence[i]
			p_b = street_sequence[i+1]
			s_a = slot_sequence[i]
			s_b = slot_sequence[i+1]

			d_ab = create_dash_through_joint(joint, p_a, p_b, s_a, s_b)
			intermediate_dashes.append(d_ab)

		# At this point we add the new dashes to the soliton's dashes
		# they must be inserted just after d_in and before d_out
		# in the obvious order: d_ij, d_jk, d_kl, etc.
		indx = new_soliton.dashes.index(d_in)
		# the dashes up to d_in, included
		dashes_before = new_soliton.dashes[0:indx+1]
		# the dashes after to d_out, included
		dashes_after = new_soliton.dashes[indx+1:]
		# update the set of dashes by adding the new ones
		new_soliton.dashes = (
				dashes_before + intermediate_dashes + dashes_after
			)

		# Finally update the soliton's growing pairs
		new_soliton.growing_pairs = (
			growing_pairs_from_dashes(new_soliton.dashes)
		)

		new_solitons.append(new_soliton)

	return new_solitons





