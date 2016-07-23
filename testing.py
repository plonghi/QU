from mcsn import MCSN, Street, Joint, BranchPoint
from solitons import (
	copy_of_soliton, join_dashes, SolitonPath, Dash, 
	growing_pairs_from_dashes, growing_clusters
)
from growth_rules import (
	bp_type_1_growth, j_type_3_growth, grow_soliton_once, grow_soliton,
)
from soliton_data import SolitonData

# # ----- Create a spectral network ------
# #
# w = MCSN()
# w.streets = {
# 	'p_1' : Street(label='p_1'),
# 	'p_2' : Street(label='p_2'),
# 	'p_3' : Street(label='p_3'),
# }
# w.branch_points = {
# 	'b_1' : BranchPoint(
# 		label='b_1', streets=[w.streets['p_1'], None, None]
# 	), 
# 	'b_2' : BranchPoint(
# 		label='b_2', streets=[w.streets['p_2']]
# 	),
# 	'b_3' : BranchPoint(
# 		label='b_3', streets=[w.streets['p_3']]
# 	),
# }
# w.joints = {'j_1': Joint(
# 	label='j_1', streets=[
# 		w.streets['p_1'], 
# 		None,
# 		w.streets['p_2'], 
# 		None,
# 		w.streets['p_3'], 
# 		None
# 	]
# )}
# w.attach_streets()
# w.check_network()
# #------ Finished creating network -------

# s1 = w.streets['p_1']
# e11 = s1.initial_point().end_point
# e12 = s1.final_point().end_point
# s2 = w.streets['p_2']
# e21 = s2.initial_point().end_point
# e22 = s2.final_point().end_point
# s3 = w.streets['p_3']
# e31 = s3.initial_point().end_point
# e32 = s3.final_point().end_point
# j1 = w.joints.values()[0]
# b1 = w.branch_points.values()[0]
# b1.print_type()
# j1.print_type()

# print '\nCreating a dash, and growing it by hand'
# dash1 = Dash(label='dash_1')
# dash1.extend_dash_along_street(street=s1, end_pt=e11, slot=0)
# dash1.print_endpoints()
# dash1.extend_dash_along_street(
# 	street=s2, end_pt='last', slot=j1.street_position(s2)[0]
# )
# dash1.print_endpoints()
# dash1.extend_dash_along_street(
# 	street=s1, end_pt='first', slot=b1.street_position(s1)[0]
# )
# dash1.print_endpoints()


# print '\nCreating a soliton, and growing it automatically'
# a1 = SolitonPath(label='a_1')
# a1.print_growing_pairs()

# a1.create(street=s1, source_pt=b1, slot=0)
# a1.print_growing_pairs()



# # ### Check how a soliton is copied:
# # a2 = copy_of_soliton(a1)
# # print a1.growing_pairs[0][0]
# # print a2.growing_pairs[0][0]
# # # for a given growing point, the streets should be the same
# # print a1.growing_pairs[0][0].street==a2.growing_pairs[0][0].street
# # # also the streets andpoint (joint/branch point) should be the same
# # print a1.growing_pairs[0][0].street_end_point==a2.growing_pairs[0][0].street_end_point
# # # also the slot on the andpoint (joint/branch point) should be the same
# # print a1.growing_pairs[0][0].slot==a2.growing_pairs[0][0].slot
# # # but the dash should be different, as it will be grown differntly for different solitons
# # print a1.growing_pairs[0][0].dash==a2.growing_pairs[0][0].dash


# print '\nthe soliton dashes are \n{}'.format([d.label for d in a1.dashes])
# print 'growth restrictions are {}\n'.format([d.growth_restriction for d in a1.dashes])

# a1.dashes[0].print_endpoints()
# print 'starting point {}'.format(a1.dashes[0].starting_point)
# print 'ending point {}\n'.format(a1.dashes[0].ending_point)

# a1.dashes[1].print_endpoints()
# print 'starting point {}'.format(a1.dashes[1].starting_point)
# print 'ending point {}'.format(a1.dashes[1].ending_point)

# print '\nsoliton starting point {}'.format(a1.starting_point)
# print 'soliton ending point {}'.format(a1.ending_point)


# cls = growing_clusters(a1)
# print '\nThe growing clusters:\n{}'.format(cls)
# ### THE FOLLOWING NOW SUBSTITUTED 
# # print '\nNow I grow soliton a1 at {}'.format(cls[0][0][0].end_point.label)
# # new_sols = bp_type_1_growth(a1, cls[0])
# new_sols = grow_soliton_once(a1)
# # print len(new_sols)
# a2=new_sols[0]

# print '\nthe soliton dashes are \n{}'.format([d.label for d in a2.dashes])
# print 'growth restrictions are {}\n'.format([d.growth_restriction for d in a2.dashes])

# a2.dashes[0].print_endpoints()
# print 'starting point {}'.format(a2.dashes[0].starting_point)
# print 'ending point {}'.format(a2.dashes[0].ending_point)


# print '\nThe soliton is now complete: {}'.format(a2.is_complete)
# print (
# 	'The soliton starting point equals the dash starting point: {}'
# 	.format(a2.starting_point==a2.complete_dash.starting_point)
# )
# print (
# 	'The soliton ending point equals the dash ending point: {}'
# 	.format(a2.ending_point==a2.complete_dash.ending_point)
# )

# # -------------------------------------------------------

# print '\n\n-------------------------------------------------------'
# print '\nCreating a second soliton, and growing it automatically'
# a3 = SolitonPath(label='a_3')
# a3.print_growing_pairs()

# a3.create(street=s1, source_pt=j1, slot=0)
# a3.print_growing_pairs()

# print '\nthe soliton dashes are \n{}'.format([d.label for d in a3.dashes])
# print 'growth restrictions are {}\n'.format([d.growth_restriction for d in a3.dashes])

# a3.dashes[0].print_endpoints()
# print 'starting point {}'.format(a3.dashes[0].starting_point)
# print 'ending point {}\n'.format(a3.dashes[0].ending_point)

# a3.dashes[1].print_endpoints()
# print 'starting point {}'.format(a3.dashes[1].starting_point)
# print 'ending point {}'.format(a3.dashes[1].ending_point)

# print '\nsoliton starting point {}'.format(a3.starting_point)
# print 'soliton ending point {}'.format(a3.ending_point)

# cls_3 = growing_clusters(a3)
# print '\nThe growing clusters:\n{}'.format(cls_3)
# ### THE FOLLOWING NOW SUBSTITUTED 
# new_sols = grow_soliton_once(a3)
# # print '\nNow I grow soliton a3 at {}'.format(cls_3[0][0][0].end_point.label)
# # new_sols = j_type_3_growth(a3, cls_3[0])
# print '\nThere are {} new solitons after growth'.format(len(new_sols))
# a4=new_sols[0]

# print '\nthe soliton dashes are \n{}'.format([d.label for d in a4.dashes])
# print 'with growth restrictions {}\n'.format([d.growth_restriction for d in a4.dashes])

# a4.dashes[0].print_endpoints()
# print 'starting point {}'.format(a4.dashes[0].starting_point)
# print 'ending point {}\n'.format(a4.dashes[0].ending_point)

# a4.dashes[1].print_endpoints()
# print 'starting point {}'.format(a4.dashes[1].starting_point)
# print 'ending point {}\n'.format(a4.dashes[1].ending_point)

# a4.dashes[2].print_endpoints()
# print 'starting point {}'.format(a4.dashes[2].starting_point)
# print 'ending point {}\n'.format(a4.dashes[2].ending_point)

# print '\nThe soliton is now complete: {}'.format(a4.is_complete)


# print '\n\n-------------------------------------------------------'
# print '\nKeep on growing automatically, with one more step'

# cls_4 = growing_clusters(a4)
# print '\nThe growing clusters:\n{}'.format(cls_4)
# ### THE FOLLOWING NOW SUBSTITUTED 
# new_sols_5 = grow_soliton_once(a4)
# # print '\nNow I grow soliton a4 at {}'.format(cls_4[0][0][0].end_point.label)
# # new_sols_5 = bp_type_1_growth(a4, cls_4[0])
# print '\nThere are {} new solitons after growth'.format(len(new_sols))
# a5=new_sols_5[0]

# print '\nthe soliton dashes are \n{}'.format([d.label for d in a5.dashes])
# print 'with growth restrictions {}\n'.format([d.growth_restriction for d in a5.dashes])

# a5.dashes[0].print_endpoints()
# print 'starting point {}'.format(a5.dashes[0].starting_point)
# print 'ending point {}\n'.format(a5.dashes[0].ending_point)

# a5.dashes[1].print_endpoints()
# print 'starting point {}'.format(a5.dashes[1].starting_point)
# print 'ending point {}\n'.format(a5.dashes[1].ending_point)

# print '\nThe soliton is now complete: {}'.format(a5.is_complete)



# print '\n\n-------------------------------------------------------'
# print '\nKeep on growing automatically, with one more step'

# cls_5 = growing_clusters(a5)
# print '\nThe growing clusters:\n{}'.format(cls_5)
# ### THE FOLLOWING NOW SUBSTITUTED 
# new_sols_6 = grow_soliton_once(a5)
# # print '\nNow I grow soliton a5 at {}'.format(cls_5[0][0][0].end_point.label)
# # new_sols_6 = bp_type_1_growth(a5, cls_5[0])
# print '\nThere are {} new solitons after growth'.format(len(new_sols))
# a6=new_sols_6[0]

# print '\nthe soliton dashes are \n{}'.format([d.label for d in a6.dashes])
# print 'with growth restrictions {}\n'.format([d.growth_restriction for d in a6.dashes])

# a6.dashes[0].print_endpoints()
# print 'starting point {}'.format(a6.dashes[0].starting_point)
# print 'ending point {}\n'.format(a6.dashes[0].ending_point)

# print '\nThe soliton is now complete: {}'.format(a6.is_complete)

# print (
# 	'The soliton starting point equals the dash starting point: {}'
# 	.format(a6.starting_point()==a6.complete_dash.starting_point)
# )
# print (
# 	'The soliton ending point equals the dash ending point: {}'
# 	.format(a6.ending_point()==a6.complete_dash.ending_point)
# )

# print '\n\n-------------------------------------------------------'
# print '\nCheck that older solitons havent been accidentally grown'
# print '\nBack to a5 soliton\n'
# a5.dashes[0].print_endpoints()
# print 'starting point {}'.format(a5.dashes[0].starting_point)
# print 'ending point {}\n'.format(a5.dashes[0].ending_point)

# a5.dashes[1].print_endpoints()
# print 'starting point {}'.format(a5.dashes[1].starting_point)
# print 'ending point {}\n'.format(a5.dashes[1].ending_point)



# print '\n\n-------------------------------------------------------'
# print 'Automatic growth method'
# a7 = SolitonPath(label='a_7')

# a7.create(street=s1, source_pt=j1, slot=0)
# new_sols = grow_soliton(a7, n_steps=3)
# print '\nThere are {} new solitons after growth'.format(len(new_sols))
# a8=new_sols[0]

# print '\nthe soliton dashes are \n{}'.format([d.label for d in a8.dashes])
# print 'with growth restrictions {}\n'.format([d.growth_restriction for d in a8.dashes])

# a8.dashes[0].print_endpoints()
# print 'starting point {}'.format(a8.dashes[0].starting_point)
# print 'ending point {}'.format(a8.dashes[0].ending_point)
# print '\nThe soliton is now complete: {}'.format(a8.is_complete)



#### TESTING LEVEL 2 ####

# # ----- Create a spectral network ------
# # 			Y-shaped network	
# #
# #	b3 ---(p3)--- j1 ---(p2)--- b2
# #				   |
# #				   |
# #				  (p1)
# #				   |
# #				   |
# #				   b1

# w = MCSN()
# w.streets = {
# 	'p_1' : Street(label='p_1'),
# 	'p_2' : Street(label='p_2'),
# 	'p_3' : Street(label='p_3'),
# }
# w.branch_points = {
# 	'b_1' : BranchPoint(
# 		label='b_1', streets=[w.streets['p_1'], None, None]
# 	), 
# 	'b_2' : BranchPoint(
# 		label='b_2', streets=[w.streets['p_2']]
# 	),
# 	'b_3' : BranchPoint(
# 		label='b_3', streets=[w.streets['p_3']]
# 	),
# }
# w.joints = {'j_1': Joint(
# 	label='j_1', streets=[
# 		w.streets['p_1'], 
# 		None,
# 		w.streets['p_2'], 
# 		None,
# 		w.streets['p_3'], 
# 		None
# 	]
# )}
# w.attach_streets()
# w.check_network()
# #------ Finished creating network -------



# # ----- Create a spectral network ------
# # 			H-shaped network	
# #
# #	b3 ---(p3)--- j2 ---(p4)--- b4
# #				   |
# #				   |
# #				  (p5)
# #				   |
# #				   |
# #	b1 ---(p1)--- j1 ---(p2)--- b2

# w = MCSN()
# w.streets = {
# 	'p_1' : Street(label='p_1'),
# 	'p_2' : Street(label='p_2'),
# 	'p_3' : Street(label='p_3'),
# 	'p_4' : Street(label='p_4'),
# 	'p_5' : Street(label='p_5'),
# }
# w.branch_points = {
# 	'b_1' : BranchPoint(
# 		label='b_1', streets=[w.streets['p_1'], None, None]
# 	), 
# 	'b_2' : BranchPoint(
# 		label='b_2', streets=[w.streets['p_2']]
# 	),
# 	'b_3' : BranchPoint(
# 		label='b_3', streets=[w.streets['p_3']]
# 	),
# 	'b_4' : BranchPoint(
# 		label='b_4', streets=[w.streets['p_4']]
# 	),
# }
# w.joints = {
# 	'j_1': Joint(
# 		label='j_1', streets=[
# 			w.streets['p_1'], 
# 			None,
# 			w.streets['p_2'], 
# 			None,
# 			w.streets['p_5'], 
# 			None
# 		]
# 	),
# 	'j_2': Joint(
# 		label='j_2', streets=[
# 			w.streets['p_5'], 
# 			None,
# 			w.streets['p_4'], 
# 			None,
# 			w.streets['p_3'], 
# 			None
# 		]
# 	),
# }

# w.attach_streets()
# w.check_network()
# #------ Finished creating network -------

# s1 = w.streets['p_1']
# s2 = w.streets['p_2']
# s3 = w.streets['p_3']
# s4 = w.streets['p_4']
# s5 = w.streets['p_5']

# print '\n\n-------------------------------------------------------'
# print 'Soliton Data'
# Q1 = SolitonData(label='Q_1', network=w , street=s1)

# Q1.initialize()
# # Q1.print_info(full_path=True)

# Q1.grow(n_steps=3)
# Q1.print_info(full_path=True)


# Q5 = SolitonData(label='Q_5', network=w , street=s5)

# Q5.initialize()
# # Q1.print_info(full_path=True)

# Q5.grow(n_steps=3)
# Q5.print_info(full_path=True)


# ### TEST TYPE 2 BRANCH POINTS

# # ----- Create a spectral network ------
# # 			H-shaped network	
# #
# #	b1 ---(p1)--- b2 ---(p2)--- b3

# w = MCSN()
# w.streets = {
# 	'p_1' : Street(label='p_1'),
# 	'p_2' : Street(label='p_2'),
# }
# w.branch_points = {
# 	'b_1' : BranchPoint(
# 		label='b_1', streets=[w.streets['p_1'], None, None]
# 	), 
# 	'b_2' : BranchPoint(
# 		label='b_2', streets=[w.streets['p_1'], w.streets['p_2'], None]
# 	),
# 	'b_3' : BranchPoint(
# 		label='b_3', streets=[w.streets['p_2']]
# 	),
# }

# w.joints = {}

# w.attach_streets()
# w.check_network()
# #------ Finished creating network -------

# s1 = w.streets['p_1']
# s2 = w.streets['p_2']

# print '\n\n-------------------------------------------------------'
# print 'Soliton Data'
# Q1 = SolitonData(label='Q_1', network=w , street=s1)

# Q1.initialize()
# # Q1.print_info(full_path=True)

# Q1.grow(n_steps=3)
# Q1.print_info(full_path=True)


# Q2 = SolitonData(label='Q_2', network=w , street=s2)

# Q2.initialize()
# # Q1.print_info(full_path=True)

# Q2.grow(n_steps=3)
# Q2.print_info(full_path=True)




# ### TEST TYPE 3 BRANCH POINTS

# # # ----- Create a spectral network ------
# # # 			H-shaped network	
# # #
# # #	b1 ---(p1)--- b4 ---(p3)--- b3
# # #				   |
# # #				   |
# # #				  (p2)
# # #				   |
# # #				   |
# # #	 			   b2

# # w = MCSN()
# # w.streets = {
# # 	'p_1' : Street(label='p_1'),
# # 	'p_2' : Street(label='p_2'),
# # 	'p_3' : Street(label='p_3'),
# # }
# # w.branch_points = {
# # 	'b_1' : BranchPoint(
# # 		label='b_1', streets=[w.streets['p_1'], None, None]
# # 	), 
# # 	'b_2' : BranchPoint(
# # 		label='b_2', streets=[w.streets['p_2']]
# # 	),
# # 	'b_3' : BranchPoint(
# # 		label='b_3', streets=[w.streets['p_3']]
# # 	),
# # 	'b_4' : BranchPoint(
# # 		label='b_3', 
# # 		streets=[w.streets['p_1'], w.streets['p_2'], w.streets['p_3']]
# # 	),
# # }

# # w.joints = {}

# # w.attach_streets()
# # w.check_network()
# # #------ Finished creating network -------

# # ----- Create a spectral network ------
# # 			  T2 network	
# #
# #			  .--- b1 --.
# #			  |    |	|
# #			  |    |	|
# #			(p3)  (p1)	(p2)
# #			  |    |	|
# #			  |    |	|
# #	 		  '--- b2 --'

# w = MCSN()
# w.streets = {
# 	'p_1' : Street(label='p_1'),
# 	'p_2' : Street(label='p_2'),
# 	'p_3' : Street(label='p_3'),
# }
# w.branch_points = {
# 	'b_1' : BranchPoint(
# 		label='b_1', 
# 		streets=[w.streets['p_1'], w.streets['p_2'], w.streets['p_3']]
# 	), 
# 	'b_2' : BranchPoint(
# 		label='b_2', 
# 		streets=[w.streets['p_1'], w.streets['p_3'], w.streets['p_2']]
# 	),
# }

# w.joints = {}

# w.attach_streets()
# w.check_network()
# #------ Finished creating network -------

# s1 = w.streets['p_1']
# s2 = w.streets['p_2']

# print '\n\n-------------------------------------------------------'
# print 'Soliton Data'
# Q1 = SolitonData(label='Q_1', network=w , street=s1)

# Q1.initialize()
# # Q1.print_info(full_path=True)

# Q1.grow(n_steps=3)
# Q1.print_info(full_path=True)






# ### TEST TYPE 4_A JOINTS

# # ----- Create a spectral network ------
# # 	  	 Peace-sign network
# #
# # 		b4		b3 		 b2
# # 		  `      |      '
# # 		   `     |     '
# # 		 (p4)  (p3)   (p2)
# # 		     `   |   '
# # 		      `  |  '
# # 		       ` | '
# # 		        j1
# # 		         |
# # 		         |
# # 		        (p1)
# # 		         |
# # 		         |
# # 		         b1

# w = MCSN()
# w.streets = {
# 	'p_1' : Street(label='p_1'),
# 	'p_2' : Street(label='p_2'),
# 	'p_3' : Street(label='p_3'),
# 	'p_4' : Street(label='p_4'),
# }
# w.branch_points = {
# 	'b_1' : BranchPoint(
# 		label='b_1', 
# 		streets=[w.streets['p_1']]
# 	), 
# 	'b_2' : BranchPoint(
# 		label='b_2', 
# 		streets=[w.streets['p_2']]
# 	),
# 	'b_3' : BranchPoint(
# 		label='b_3', 
# 		streets=[w.streets['p_3']]
# 	),
# 	'b_4' : BranchPoint(
# 		label='b_4', 
# 		streets=[w.streets['p_4']]
# 	),
# }

# w.joints = {
# 	'j_1': Joint(
# 		label='j_1', streets=[
# 			w.streets['p_1'], 
# 			None,
# 			w.streets['p_2'], 	
# 			w.streets['p_3'], 
# 			w.streets['p_4'], 
# 			None
# 		]
# 	),
# }

# w.attach_streets()
# w.check_network()
# #------ Finished creating network -------

# s1 = w.streets['p_1']
# s2 = w.streets['p_2']
# s3 = w.streets['p_3']
# s4 = w.streets['p_4']


# print '\n\n-------------------------------------------------------'
# print 'Soliton Data'
# Q1 = SolitonData(label='Q_1', network=w , street=s1)
# Q1.initialize()
# # Q1.print_info(full_path=True)
# Q1.grow(n_steps=3)
# Q1.print_info(full_path=True)


# Q4 = SolitonData(label='Q_4', network=w , street=s4)
# Q4.initialize()
# # Q1.print_info(full_path=True)
# Q4.grow(n_steps=3)
# Q4.print_info(full_path=True)

# Q3 = SolitonData(label='Q_3', network=w , street=s3)
# Q3.initialize()
# # Q1.print_info(full_path=True)
# Q3.grow(n_steps=3)
# Q3.print_info(full_path=True)




# ### TEST TYPE 4_B JOINTS

# # ----- Create a spectral network ------
# # 	  	 Cross-sign network
# #
# # 		b4		 		 b3
# # 		  `             '
# # 		   `           '
# # 		   (p4)     (p3)
# # 		     `       '
# # 		      `     '
# # 		       `   '
# # 		        ` '
# #				j1
# # 		        ' `
# # 		       '   `
# # 		      '     `
# # 		   (p1)     (p2)
# # 		    '         `
# # 		   '           `
# # 		 b1             b2

# w = MCSN()
# w.streets = {
# 	'p_1' : Street(label='p_1'),
# 	'p_2' : Street(label='p_2'),
# 	'p_3' : Street(label='p_3'),
# 	'p_4' : Street(label='p_4'),
# }
# w.branch_points = {
# 	'b_1' : BranchPoint(
# 		label='b_1', 
# 		streets=[w.streets['p_1']]
# 	), 
# 	'b_2' : BranchPoint(
# 		label='b_2', 
# 		streets=[w.streets['p_2']]
# 	),
# 	'b_3' : BranchPoint(
# 		label='b_3', 
# 		streets=[w.streets['p_3']]
# 	),
# 	'b_4' : BranchPoint(
# 		label='b_4', 
# 		streets=[w.streets['p_4']]
# 	),
# }

# w.joints = {
# 	'j_1': Joint(
# 		label='j_1', streets=[
# 			w.streets['p_1'], 
# 			w.streets['p_2'], 
# 			None,	
# 			w.streets['p_3'], 
# 			w.streets['p_4'], 
# 			None
# 		]
# 	),
# }

# w.attach_streets()
# w.check_network()
# #------ Finished creating network -------

# s1 = w.streets['p_1']
# s2 = w.streets['p_2']
# s3 = w.streets['p_3']
# s4 = w.streets['p_4']


# print '\n\n-------------------------------------------------------'
# print 'Soliton Data'
# Q1 = SolitonData(label='Q_1', network=w , street=s1)
# Q1.initialize()
# # Q1.print_info(full_path=True)
# Q1.grow(n_steps=3)
# Q1.print_info(full_path=True)

# Q2 = SolitonData(label='Q_2', network=w , street=s2)
# Q2.initialize()
# # Q1.print_info(full_path=True)
# Q2.grow(n_steps=3)
# Q2.print_info(full_path=True)

