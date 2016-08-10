from mcsn import MCSN, Street, Joint, BranchPoint
from solitons import (
    copy_of_soliton, join_dashes_at_branch_point, SolitonPath, Dash, 
    growing_pairs_from_dashes, growing_clusters
)
from growth_rules import (
    # bp_type_1_growth, j_type_3_growth, 
    grow_soliton_once, grow_soliton,
)
from soliton_data import SolitonData

from intersections import get_dash_nodes, compute_self_intersections

# # ----- Create a spectral network ------
# #
# w = MCSN()
# w.streets = {
#   'p_1' : Street(label='p_1'),
#   'p_2' : Street(label='p_2'),
#   'p_3' : Street(label='p_3'),
# }
# w.branch_points = {
#   'b_1' : BranchPoint(
#       label='b_1', streets=[w.streets['p_1'], None, None]
#   ), 
#   'b_2' : BranchPoint(
#       label='b_2', streets=[w.streets['p_2']]
#   ),
#   'b_3' : BranchPoint(
#       label='b_3', streets=[w.streets['p_3']]
#   ),
# }
# w.joints = {'j_1': Joint(
#   label='j_1', streets=[
#       w.streets['p_1'], 
#       None,
#       w.streets['p_2'], 
#       None,
#       w.streets['p_3'], 
#       None
#   ]
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
#   street=s2, end_pt='last', slot=j1.street_position(s2)[0]
# )
# dash1.print_endpoints()
# dash1.extend_dash_along_street(
#   street=s1, end_pt='first', slot=b1.street_position(s1)[0]
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
#   'The soliton starting point equals the dash starting point: {}'
#   .format(a2.starting_point==a2.complete_dash.starting_point)
# )
# print (
#   'The soliton ending point equals the dash ending point: {}'
#   .format(a2.ending_point==a2.complete_dash.ending_point)
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
#   'The soliton starting point equals the dash starting point: {}'
#   .format(a6.starting_point()==a6.complete_dash.starting_point)
# )
# print (
#   'The soliton ending point equals the dash ending point: {}'
#   .format(a6.ending_point()==a6.complete_dash.ending_point)
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
# #             Y-shaped network    
# #
# # b3 ---(p3)--- j1 ---(p2)--- b2
# #                |
# #                |
# #               (p1)
# #                |
# #                |
# #                b1

# w = MCSN()
# w.streets = {
#   'p_1' : Street(label='p_1'),
#   'p_2' : Street(label='p_2'),
#   'p_3' : Street(label='p_3'),
# }
# w.branch_points = {
#   'b_1' : BranchPoint(
#       label='b_1', streets=[w.streets['p_1'], None, None]
#   ), 
#   'b_2' : BranchPoint(
#       label='b_2', streets=[w.streets['p_2']]
#   ),
#   'b_3' : BranchPoint(
#       label='b_3', streets=[w.streets['p_3']]
#   ),
# }
# w.joints = {'j_1': Joint(
#   label='j_1', streets=[
#       w.streets['p_1'], 
#       None,
#       w.streets['p_2'], 
#       None,
#       w.streets['p_3'], 
#       None
#   ]
# )}
# w.attach_streets()
# w.check_network()
# #------ Finished creating network -------



# # ----- Create a spectral network ------
# #             H-shaped network    
# #
# # b3 ---(p3)--- j2 ---(p4)--- b4
# #                |
# #                |
# #               (p5)
# #                |
# #                |
# # b1 ---(p1)--- j1 ---(p2)--- b2

# w = MCSN()
# w.streets = {
#   'p_1' : Street(label='p_1'),
#   'p_2' : Street(label='p_2'),
#   'p_3' : Street(label='p_3'),
#   'p_4' : Street(label='p_4'),
#   'p_5' : Street(label='p_5'),
# }
# w.branch_points = {
#   'b_1' : BranchPoint(
#       label='b_1', streets=[w.streets['p_1'], None, None]
#   ), 
#   'b_2' : BranchPoint(
#       label='b_2', streets=[w.streets['p_2']]
#   ),
#   'b_3' : BranchPoint(
#       label='b_3', streets=[w.streets['p_3']]
#   ),
#   'b_4' : BranchPoint(
#       label='b_4', streets=[w.streets['p_4']]
#   ),
# }
# w.joints = {
#   'j_1': Joint(
#       label='j_1', streets=[
#           w.streets['p_1'], 
#           None,
#           w.streets['p_2'], 
#           None,
#           w.streets['p_5'], 
#           None
#       ]
#   ),
#   'j_2': Joint(
#       label='j_2', streets=[
#           w.streets['p_5'], 
#           None,
#           w.streets['p_4'], 
#           None,
#           w.streets['p_3'], 
#           None
#       ]
#   ),
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
# #             H-shaped network    
# #
# # b1 ---(p1)--- b2 ---(p2)--- b3

# w = MCSN()
# w.streets = {
#   'p_1' : Street(label='p_1'),
#   'p_2' : Street(label='p_2'),
# }
# w.branch_points = {
#   'b_1' : BranchPoint(
#       label='b_1', streets=[w.streets['p_1'], None, None]
#   ), 
#   'b_2' : BranchPoint(
#       label='b_2', streets=[w.streets['p_1'], w.streets['p_2'], None]
#   ),
#   'b_3' : BranchPoint(
#       label='b_3', streets=[w.streets['p_2']]
#   ),
# }

# w.joints = {}

# w.homology_classes = {
#   'gamma_1' : ['p_1'],
#   'gamma_2' : ['p_2'],
# }

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

# print '\nHomology classes of closed solitons in Q1'
# Q1.compute_closed_solitons()
# print [sol.homology_class for sol in Q1.closed_solitons]

# print '\nHomology classes of closed solitons in Q2'
# Q2.compute_closed_solitons()
# print [sol.homology_class for sol in Q2.closed_solitons]



# ## TEST TYPE 3 BRANCH POINTS

# # # ----- Create a spectral network ------
# # #           H-shaped network    
# # #
# # #   b1 ---(p1)--- b4 ---(p3)--- b3
# # #                  |
# # #                  |
# # #                 (p2)
# # #                  |
# # #                  |
# # #                  b2

# # w = MCSN()
# # w.streets = {
# #     'p_1' : Street(label='p_1'),
# #     'p_2' : Street(label='p_2'),
# #     'p_3' : Street(label='p_3'),
# # }
# # w.branch_points = {
# #     'b_1' : BranchPoint(
# #         label='b_1', streets=[w.streets['p_1'], None, None]
# #     ), 
# #     'b_2' : BranchPoint(
# #         label='b_2', streets=[w.streets['p_2']]
# #     ),
# #     'b_3' : BranchPoint(
# #         label='b_3', streets=[w.streets['p_3']]
# #     ),
# #     'b_4' : BranchPoint(
# #         label='b_3', 
# #         streets=[w.streets['p_1'], w.streets['p_2'], w.streets['p_3']]
# #     ),
# # }

# # w.joints = {}

# # w.homology_classes = {
# #     'gamma_1' : ['p_1'],
# #     'gamma_2' : ['p_2'],
# #     'gamma_3' : ['p_3'],
# # }

# # w.attach_streets()
# # w.check_network()

# # #------ Finished creating network -------

# # ----- Create a spectral network ------
# #               T2 network    
# #
# #           .--- b1 --.
# #           |    |    |
# #           |    |    |
# #         (p3)  (p1)  (p2)
# #           |    |    |
# #           |    |    |
# #           '--- b2 --'


# streets = ['p_1', 'p_2', 'p_3']
# branch_points = (
#   {
#       'b_1': ['p_1', 'p_2', 'p_3'],
#       'b_2': ['p_1', 'p_3', 'p_2'],
#   }
# )

# joints = {}

# homology_classes = {
#   'gamma_1' : ['p_1'],
#   'gamma_2' : ['p_2'],
#   'gamma_3' : ['p_3'],
# }

# w = MCSN(
#   branch_points=branch_points, 
#   streets=streets, 
#   joints=joints, 
#   homology_classes=homology_classes
# )

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

# print '\nHomology classes of closed solitons in Q1'
# Q1.compute_closed_solitons()
# print [sol.homology_class for sol in Q1.closed_solitons]

# tr = w.trivial_homology_class
# h1 = w.basis_homology_classes['gamma_1']
# print h1.atomic_labels
# ### TEST TYPE 4_A JOINTS

# # ----- Create a spectral network ------
# #          Peace-sign network
# #
# #         b4      b3       b2
# #           `      |      '
# #            `     |     '
# #          (p4)  (p3)   (p2)
# #              `   |   '
# #               `  |  '
# #                ` | '
# #                 j1
# #                  |
# #                  |
# #                 (p1)
# #                  |
# #                  |
# #                  b1

# w = MCSN()
# w.streets = {
#   'p_1' : Street(label='p_1'),
#   'p_2' : Street(label='p_2'),
#   'p_3' : Street(label='p_3'),
#   'p_4' : Street(label='p_4'),
# }
# w.branch_points = {
#   'b_1' : BranchPoint(
#       label='b_1', 
#       streets=[w.streets['p_1']]
#   ), 
#   'b_2' : BranchPoint(
#       label='b_2', 
#       streets=[w.streets['p_2']]
#   ),
#   'b_3' : BranchPoint(
#       label='b_3', 
#       streets=[w.streets['p_3']]
#   ),
#   'b_4' : BranchPoint(
#       label='b_4', 
#       streets=[w.streets['p_4']]
#   ),
# }

# w.joints = {
#   'j_1': Joint(
#       label='j_1', streets=[
#           w.streets['p_1'], 
#           None,
#           w.streets['p_2'],   
#           w.streets['p_3'], 
#           w.streets['p_4'], 
#           None
#       ]
#   ),
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
# #          Cross-sign network
# #
# #         b4               b3
# #           `             '
# #            `           '
# #            (p4)     (p3)
# #              `       '
# #               `     '
# #                `   '
# #                 ` '
# #             j1
# #                 ' `
# #                '   `
# #               '     `
# #            (p1)     (p2)
# #             '         `
# #            '           `
# #          b1             b2

# w = MCSN()
# w.streets = {
#   'p_1' : Street(label='p_1'),
#   'p_2' : Street(label='p_2'),
#   'p_3' : Street(label='p_3'),
#   'p_4' : Street(label='p_4'),
# }
# w.branch_points = {
#   'b_1' : BranchPoint(
#       label='b_1', 
#       streets=[w.streets['p_1']]
#   ), 
#   'b_2' : BranchPoint(
#       label='b_2', 
#       streets=[w.streets['p_2']]
#   ),
#   'b_3' : BranchPoint(
#       label='b_3', 
#       streets=[w.streets['p_3']]
#   ),
#   'b_4' : BranchPoint(
#       label='b_4', 
#       streets=[w.streets['p_4']]
#   ),
# }

# w.joints = {
#   'j_1': Joint(
#       label='j_1', streets=[
#           w.streets['p_1'], 
#           w.streets['p_2'], 
#           None,   
#           w.streets['p_3'], 
#           w.streets['p_4'], 
#           None
#       ]
#   ),
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




# ### TEST TYPE 5 JOINTS

# # ----- Create a spectral network ------
# #          Cross-sign network
# #
# #         b4               b3
# #           `             '
# #            `           '
# #            (p4)     (p3)
# #              `       '
# #               `     '
# #                `   '
# #                 ` '
# #             j1
# #                 '|`
# #                ' | `
# #               '  |  `
# #            (p1) (p5)  (p2)
# #             '    |    `
# #            '     |     `
# #          b1      b5     b2

# w = MCSN()
# w.streets = {
#   'p_1' : Street(label='p_1'),
#   'p_2' : Street(label='p_2'),
#   'p_3' : Street(label='p_3'),
#   'p_4' : Street(label='p_4'),
#   'p_5' : Street(label='p_5'),
# }
# w.branch_points = {
#   'b_1' : BranchPoint(
#       label='b_1', 
#       streets=[w.streets['p_1']]
#   ), 
#   'b_2' : BranchPoint(
#       label='b_2', 
#       streets=[w.streets['p_2']]
#   ),
#   'b_3' : BranchPoint(
#       label='b_3', 
#       streets=[w.streets['p_3']]
#   ),
#   'b_4' : BranchPoint(
#       label='b_4', 
#       streets=[w.streets['p_4']]
#   ),
#   'b_5' : BranchPoint(
#       label='b_5', 
#       streets=[w.streets['p_5']]
#   ),
# }

# w.joints = {
#   'j_1': Joint(
#       label='j_1', streets=[
#           w.streets['p_1'], 
#           w.streets['p_5'], 
#           w.streets['p_2'], 
#           w.streets['p_3'], 
#           None,
#           w.streets['p_4'], 
#       ]
#   ),
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

# Q2 = SolitonData(label='Q_2', network=w , street=s2)
# Q2.initialize()
# # Q1.print_info(full_path=True)
# Q2.grow(n_steps=3)
# Q2.print_info(full_path=True)

# Q5 = SolitonData(label='Q_5', network=w , street=s5)
# Q5.initialize()
# # Q1.print_info(full_path=True)
# Q5.grow(n_steps=3)
# Q5.print_info(full_path=True)

# Q3 = SolitonData(label='Q_3', network=w , street=s3)
# Q3.initialize()
# # Q1.print_info(full_path=True)
# Q3.grow(n_steps=3)
# Q3.print_info(full_path=True)





# ### TEST the 1-herd

# # ----- Create a spectral network ------
# #          1-herd network



# streets = ['p_1','p_2','q_1','q_2']
# branch_points = {
#   'b_1' : ['p_1'],
#   'b_2' : ['q_1'],
#   'b_3' : ['p_2'],
#   'b_4' : ['q_2']
# }

# joints = {
#   'j_1': ['p_1', None, 'q_1', 'p_2', None, 'q_2' ]
# }

# homology_classes = {
#   'gamma_1' : ['p_1', 'p_2'],
#   'gamma_2' : ['q_1', 'q_2'],
# }

# w = MCSN(
#   branch_points=branch_points, 
#   streets=streets, 
#   joints=joints, 
#   homology_classes=homology_classes
# )

# #------ Finished creating network -------

# s1 = w.streets['p_1']
# s2 = w.streets['p_2']
# r1 = w.streets['q_1']
# r2 = w.streets['q_2']



# print '\n\n-------------------------------------------------------'
# print 'Soliton Data'
# Q1 = SolitonData(label='Q_1', network=w , street=s1)
# Q1.initialize()
# # Q1.print_info(full_path=True)
# Q1.grow(n_steps=10)
# Q1.print_info(full_path=True)

# Q2 = SolitonData(label='Q_2', network=w , street=r1)
# Q2.initialize()
# # Q1.print_info(full_path=True)
# Q2.grow(n_steps=10)
# Q2.print_info(full_path=True)


# Q3 = SolitonData(label='Q_3', network=w , street=s2)
# Q3.initialize()
# # Q1.print_info(full_path=True)
# Q3.grow(n_steps=10)
# Q3.print_info(full_path=True)





# ### TEST the 2-herd

# # ----- Create a spectral network ------
# #          2-herd network



# streets = ['p_1','p_2','p_3','q_1','q_2','q_3']
# branch_points = {
#   'b_1' : ['p_1'],
#   'b_2' : ['q_1'],
#   'b_3' : ['p_3'],
#   'b_4' : ['q_3']
# }

# joints = {
#   'j_1': ['p_1', None, 'q_1', 'p_2', None, 'q_2' ],
#   'j_2': ['p_2', None, 'q_2', 'p_3', None, 'q_3' ]
# }

# homology_classes = {
#   'gamma_1' : ['p_1', 'p_2', 'p_3'],
#   'gamma_2' : ['q_1', 'q_2', 'q_3'],
# }

# w = MCSN(
#   branch_points=branch_points, 
#   streets=streets, 
#   joints=joints, 
#   homology_classes=homology_classes
# )

# #------ Finished creating network -------

# s1 = w.streets['p_1']
# s2 = w.streets['p_2']
# r1 = w.streets['q_1']
# r2 = w.streets['q_2']



# print '\n\n-------------------------------------------------------'
# print 'Soliton Data'
# Q1 = SolitonData(label='Q_1', network=w , street=s1)
# Q1.initialize()
# # Q1.print_info(full_path=True)
# Q1.grow(n_steps=10)
# Q1.print_info(full_path=True)

# Q2 = SolitonData(label='Q_2', network=w , street=r1)
# Q2.initialize()
# # Q1.print_info(full_path=True)
# Q2.grow(n_steps=10)
# Q2.print_info(full_path=True)


# Q3 = SolitonData(label='Q_3', network=w , street=s2)
# Q3.initialize()
# # Q1.print_info(full_path=True)
# Q3.grow(n_steps=10)
# Q3.print_info(full_path=True)





# ### TEST the 3-herd

# # ----- Create a spectral network ------
# #          3-herd network



# streets = ['p_1','p_2','p_3','p_4','q_1','q_2','q_3','q_4']
# branch_points = {
#   'b_1' : ['p_1'],
#   'b_2' : ['q_1'],
#   'b_3' : ['p_4'],
#   'b_4' : ['q_4']
# }

# joints = {
#   'j_1': ['p_1', None, 'q_1', 'p_2', None, 'q_2' ],
#   'j_2': ['p_2', None, 'q_2', 'p_3', None, 'q_3' ],
#   'j_3': ['p_3', None, 'q_3', 'p_4', None, 'q_4' ]
# }

# homology_classes = {
#   'gamma_1' : ['p_1', 'p_2', 'p_3', 'p_4'],
#   'gamma_2' : ['q_1', 'q_2', 'q_3', 'q_4'],
# }

# w = MCSN(
#   branch_points=branch_points, 
#   streets=streets, 
#   joints=joints, 
#   homology_classes=homology_classes
# )

# #------ Finished creating network -------

# s1 = w.streets['p_1']
# s2 = w.streets['p_2']
# r1 = w.streets['q_1']
# r2 = w.streets['q_2']



# print '\n\n-------------------------------------------------------'
# print 'Soliton Data'
# Q1 = SolitonData(label='Q_1', network=w , street=s1)
# Q1.initialize()
# # Q1.print_info(full_path=True)
# Q1.grow(n_steps=10)
# Q1.print_info(full_path=True)

# Q2 = SolitonData(label='Q_2', network=w , street=r1)
# Q2.initialize()
# # Q1.print_info(full_path=True)
# Q2.grow(n_steps=10)
# Q2.print_info(full_path=True)


# Q3 = SolitonData(label='Q_3', network=w , street=s2)
# Q3.initialize()
# # Q1.print_info(full_path=True)
# Q3.grow(n_steps=10)
# Q3.print_info(full_path=True)





# ### T3

# # ----- Create a spectral network ------
# #                T3 network



# streets = ['p_'+str(i) for i in range(1, 13)]
# branch_points = {
#   'b_1' : ['p_1', 'p_7', 'p_8'],
#   'b_2' : ['p_2', 'p_10', 'p_12'],
#   'b_3' : ['p_3', 'p_11', 'p_9'],
#   'b_4' : ['p_4','p_8','p_7'],
#   'b_5' : ['p_5', 'p_9', 'p_11'],
#   'b_6' : ['p_6', 'p_12', 'p_10']
# }

# joints = {
#   'j_1': ['p_1', None, 'p_2', None, 'p_3', None ],
#   'j_2': ['p_4', None, 'p_5', None, 'p_6', None],
# }

# homology_classes = {
#   'gamma_1' : ['p_1', 'p_2', 'p_3'],
#   'gamma_2' : ['p_4', 'p_5', 'p_6'],
#   'gamma_3' : ['p_7'],
#   'gamma_4' : ['p_8'],
#   'gamma_5' : ['p_9'],
#   'gamma_6' : ['p_10'],
#   'gamma_7' : ['p_11'],
#   'gamma_8' : ['p_12'],
# }

# w = MCSN(
#   branch_points=branch_points, 
#   streets=streets, 
#   joints=joints, 
#   homology_classes=homology_classes
# )

# #------ Finished creating network -------

# s1 = w.streets['p_1']
# s2 = w.streets['p_2']
# s7 = w.streets['p_7']


# # print '\n\n-------------------------------------------------------'
# print '\nSoliton Data'
# Q1 = SolitonData(label='Q_1', network=w , street=s1)
# Q1.initialize()
# # # Q1.print_info(full_path=True)
# Q1.grow(n_steps=8)
# Q1.print_info(full_path=True)


# # # print '\n\n-------------------------------------------------------'
# # print '\nSoliton Data'
# # Q7 = SolitonData(label='Q_7', network=w , street=s7)
# # Q7.initialize()
# # # # Q1.print_info(full_path=True)
# # Q7.grow(n_steps=9)
# # Q7.print_info(full_path=True)




# ### SU(2) N_f=4

# # ----- Create a spectral network ------
# #            SU(2) N_f=4 network


# streets = ['p_'+str(i+1) for i in range(6)]
# branch_points = {
#   'b_1': ['p_1', 'p_2', 'p_3'],
#   'b_2': ['p_2', 'p_4', 'p_5'],
#   'b_3': ['p_1', 'p_6', 'p_4'],
#   'b_4': ['p_3', 'p_5', 'p_6'],
# }

# joints = {}

# homology_classes = {
#   'gamma_1' : ['p_1'],
#   'gamma_2' : ['p_2'],
#   'gamma_3' : ['p_3'],
#   'gamma_4' : ['p_4'],
#   'gamma_5' : ['p_5'],
#   'gamma_6' : ['p_6'],
# }

# w = MCSN(
#   branch_points=branch_points, 
#   streets=streets, 
#   joints=joints, 
#   homology_classes=homology_classes
# )

# #------ Finished creating network -------

# s1 = w.streets['p_1']
# s2 = w.streets['p_2']


# # print '\n\n-------------------------------------------------------'
# print '\nSoliton Data'
# Q1 = SolitonData(label='Q_1', network=w , street=s1)
# Q1.initialize()
# # # Q1.print_info(full_path=True)
# Q1.grow(n_steps=8)
# Q1.print_info(full_path=True)


# # # print '\n\n-------------------------------------------------------'
# # print '\nSoliton Data'
# # Q7 = SolitonData(label='Q_7', network=w , street=s7)
# # Q7.initialize()
# # # # Q1.print_info(full_path=True)
# # Q7.grow(n_steps=9)
# # Q7.print_info(full_path=True)




### AD_3 theory

# ----- Create a spectral network ------
#              AD_3 theory


streets = ['p_' + str(i + 1) for i in range(2)]
branch_points = {
    'b_1': ['p_1'],
    'b_2': ['p_1', 'p_2'],
    'b_3': ['p_2'],
}

joints = {}

homology_classes = {
    'gamma_1' : ['p_2'],
    'gamma_2' : ['p_1'],
}

w = MCSN(
    branch_points=branch_points, 
    streets=streets, 
    joints=joints, 
    homology_classes=homology_classes
)

#------ Finished creating network -------

s1 = w.streets['p_1']
s2 = w.streets['p_2']


# print '\n\n-------------------------------------------------------'
print '\nSoliton Data'
Q1 = SolitonData(label='Q_1', network=w , street=s1)
Q1.initialize()
# # Q1.print_info(full_path=True)
Q1.grow(n_steps=8)
Q1.print_info(full_path=True)

Q2 = SolitonData(label='Q_2', network=w, street=s2)
Q2.initialize()
# # Q1.print_info(full_path=True)
Q2.grow(n_steps=8)
Q2.print_info(full_path=True)


w.print_info()

print '\n\n'
print Q1.closed_solitons[1].dash.path
dash_nodes = get_dash_nodes(Q1.closed_solitons[1].dash, is_closed_soliton=True)
print dash_nodes
print [d_n.node.__class__.__name__ for d_n in dash_nodes]
for node in dash_nodes:
    node.print_info()

k = compute_self_intersections(Q1.closed_solitons[1].dash, is_closed_soliton=True)
print k
# # print '\n\n-------------------------------------------------------'
# print '\nSoliton Data'
# Q7 = SolitonData(label='Q_7', network=w , street=s7)
# Q7.initialize()
# # # Q1.print_info(full_path=True)
# Q7.grow(n_steps=9)
# Q7.print_info(full_path=True)
