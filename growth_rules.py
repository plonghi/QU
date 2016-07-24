from solitons import (
	copy_of_soliton, join_dashes_at_branch_point, SolitonPath, Dash, 
	growing_pairs_from_dashes, growing_clusters, find_corresponding_cluster,
	find_corresponding_pair,
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

	# A soliton must be grown at each one of its growth clusters.
	# Some of these evolutions may return more than one outcome, 
	# for example, this would happen at actual 6-way joints.
	# So we go through the growth clusters, and grow each of them,
	# while allowing for the possibility that, as we evolve, new 
	# solitons will pop up.
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
	if n_steps == 0:
		return old_solitons

	elif n_steps > 0:
		for i in range(n_steps):
			new_solitons = []
			for sol in old_solitons:
				new_solitons += grow_soliton_once(sol)
			old_solitons = new_solitons

		return new_solitons

	else:
		ValueError


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
	# A soliton must be grown at each one of its growth clusters.
	# Some of these evolutions may return more than one outcome, 
	# for example, this would happen at actual 6-way joints.
	# So we go through the growth clusters, and grow each of them,
	# while allowing for the possibility that, as we evolve, new 
	# solitons will pop up.
	old_solitons = [old_soliton]
	new_solitons = []
	joint = old_cluster[0][0].end_point

	for old_growing_pair in old_cluster:
		# for each growing pair of the cluster we do the evolution procedure
		# at the joint. Each time, this may return more than one new soliton.
		# Then we must keep evolving all of the new ones for the next growing 
		# pair.
		new_solitons = []
		for sol in old_solitons:
			# must now identify which growing pair of each soliton 
			# is the one corresponding to 'old_growing_pair'
			sol_growing_pair = find_corresponding_pair(
				sol, old_growing_pair, multi=False
			)[0]
			sol_slot = sol_growing_pair[0].slot
			if sol_slot % 2 == 0:
				parity_reversal = False
			else:
				parity_reversal = True
			av_slots = joint.available_slots
			# let p1 be the street from which 
			# the soliton is approaching the joint (for this growing pair).

			# the slot of each street
			[s_1, s_2, s_3, s_4, s_5, s_6] = (
				[(sol_slot + i) % 6 for i in range(6)]
			)
			s_k = [s_1, s_2, s_3, s_4, s_5, s_6] 
			
			# the streets
			[p1, p2, p3, p4, p5, p6] = (
				[joint.streets[s_k[i]] for i in range(6)]
			)
			
			# then map which streets are 'available'
			pk_is_av = [False for i in range(6)]
			for k in range(6):
				if s_k[k] in av_slots:
					pk_is_av[k] = True 

			[p1_is_av, p2_is_av, p3_is_av, p4_is_av, p5_is_av, p6_is_av] = (
				pk_is_av
			)
			# tau1 > nu4
			if (
				p4_is_av is True
			):
				# Just straight propagation, along street 
				# at the opposite slot (+3)
				new_solitons.append(soliton_propagation_street_sequence(
					sol, sol_growing_pair, [3], 
					parity_reversed=parity_reversal
				))

			# tau1 > nu3 nu5
			if (
				p3_is_av is True and
				p5_is_av is True
			):
				# Y-propagation on streets 3 (+2) then 5 (+4)
				new_solitons.append(soliton_propagation_street_sequence(
					sol, sol_growing_pair, [2, 4], 
					parity_reversed=parity_reversal
				))

			# tau1 > nu3 nu6 nu4 
			if (
				p3_is_av is True and
				p6_is_av is True and
				p4_is_av is True
			):
				# X-propagation on streets 3 (+2), 6 (+5), 4 (+3)
				new_solitons.append(soliton_propagation_street_sequence(
					sol, sol_growing_pair, [2, 5, 3], 
					parity_reversed=parity_reversal
				))

			# tau1 > nu3 nu5 nu1 nu4 
			if (
				p3_is_av is True and
				p5_is_av is True and
				p1_is_av is True and
				p4_is_av is True
			):
				# propagation on streets 3 (+2), 5 (+4), 1 (+0), 4 (+3)
				new_solitons.append(soliton_propagation_street_sequence(
					sol, sol_growing_pair, [2, 4, 0, 3], 
					parity_reversed=parity_reversal
				))

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
				new_solitons.append(soliton_propagation_street_sequence(
					sol, sol_growing_pair, [2, 4, 1, 5, 3], 
					parity_reversed=parity_reversal
				))

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
				new_solitons.append(soliton_propagation_street_sequence(
					sol, sol_growing_pair, [2, 4, 0, 2, 5, 3], 
					parity_reversed=parity_reversal
				))
			
			### TODO: if all slots are available, add more iterations!
			old_solitons = new_solitons

	return new_solitons


def soliton_capping_off(old_soliton, old_cluster):
	"""
	Handles capping-off of a soliton at a branch point, 
	for one or more growing pairs.
	Can only return one soliton, not more, but it will return a list
	for convenience of integration of this function with the rest.
	"""

	new_soliton = copy_of_soliton(old_soliton)
	n_growing_pairs = len(old_cluster)
	reference_pair = old_cluster[0]
	# for the new soliton, we now look for all growing pairs
	# that end on this branch point, and cap them off one by one.
	for i in range(n_growing_pairs):
		# Now for the growing pair p of DashEndpoints we 
		# identify the new corresponding growing pair in 
		# the new soliton. 
		# At first there may be more than one, but as we cap 
		# them off, the number will start decreasing.
		# We always cap off the first one in the list at each step.
		new_p = find_corresponding_pair(
			new_soliton, reference_pair, multi=True
		)[0] 
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

	return [new_soliton]


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


def soliton_propagation_street_sequence(
	old_soliton, old_growing_pair, rel_street_sequence, parity_reversed=False
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

	joint = old_growing_pair[0].end_point

	new_soliton = copy_of_soliton(old_soliton)
	# Now for the growing pair p of DashEndpoints we 
	# identify the new corresponding growing pair in 
	# the new soliton
	old_pair_index = old_soliton.growing_pairs.index(old_growing_pair)
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

	return new_soliton





