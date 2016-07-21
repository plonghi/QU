from mcsn import MCSN, Street, Joint, BranchPoint
from solitons import (
	copy_of_soliton, join_dashes, SolitonPath, Dash, 
	growing_pairs_from_dashes, growing_clusters
)
from growth_rules import bp_type_1_growth, j_type_3_growth

# ----- Create a spectral network ------
#
w = MCSN()
w.streets = {
	'p_1' : Street(label='p_1'),
	'p_2' : Street(label='p_2'),
	'p_3' : Street(label='p_3'),
}
w.branch_points = {
	'b_1' : BranchPoint(
		label='b_1', streets=[w.streets['p_1'], None, None]
	), 
	'b_2' : BranchPoint(
		label='b_2', streets=[w.streets['p_2']]
	),
	'b_3' : BranchPoint(
		label='b_3', streets=[w.streets['p_3']]
	),
}
w.joints = {'j_1': Joint(
	label='j_1', streets=[
		w.streets['p_1'], 
		None,
		w.streets['p_2'], 
		None,
		w.streets['p_3'], 
		None
	]
)}
w.attach_streets()
w.check_network()
#------ Finished creating network -------

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


print '\nCreating a soliton, and growing it automatically'
a1 = SolitonPath(label='a_1')
a1.print_growing_pairs()

a1.create(street=s1, source_pt=b1, slot=0)
a1.print_growing_pairs()



# ### Check how a soliton is copied:
# a2 = copy_of_soliton(a1)
# print a1.growing_pairs[0][0]
# print a2.growing_pairs[0][0]
# # for a given growing point, the streets should be the same
# print a1.growing_pairs[0][0].street==a2.growing_pairs[0][0].street
# # also the streets andpoint (joint/branch point) should be the same
# print a1.growing_pairs[0][0].street_end_point==a2.growing_pairs[0][0].street_end_point
# # also the slot on the andpoint (joint/branch point) should be the same
# print a1.growing_pairs[0][0].slot==a2.growing_pairs[0][0].slot
# # but the dash should be different, as it will be grown differntly for different solitons
# print a1.growing_pairs[0][0].dash==a2.growing_pairs[0][0].dash


print '\nthe soliton dashes are \n{}'.format([d.label for d in a1.dashes])
print 'growth restrictions are {}\n'.format([d.growth_restriction for d in a1.dashes])

a1.dashes[0].print_endpoints()
print 'starting point {}'.format(a1.dashes[0].starting_point)
print 'ending point {}\n'.format(a1.dashes[0].ending_point)

a1.dashes[1].print_endpoints()
print 'starting point {}'.format(a1.dashes[1].starting_point)
print 'ending point {}'.format(a1.dashes[1].ending_point)

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

print '\nthe soliton dashes are \n{}'.format([d.label for d in a2.dashes])
print 'growth restrictions are {}\n'.format([d.growth_restriction for d in a2.dashes])

a2.dashes[0].print_endpoints()
print 'starting point {}'.format(a2.dashes[0].starting_point)
print 'ending point {}'.format(a2.dashes[0].ending_point)


print '\nThe soliton is now complete: {}'.format(a2.is_complete)
print (
	'The soliton starting point equals the dash starting point: {}'
	.format(a2.starting_point==a2.complete_dash.starting_point)
)
print (
	'The soliton ending point equals the dash ending point: {}'
	.format(a2.ending_point==a2.complete_dash.ending_point)
)

# -------------------------------------------------------

print '\n\n-------------------------------------------------------'
print '\nCreating a second soliton, and growing it automatically'
a3 = SolitonPath(label='a_3')
a3.print_growing_pairs()

a3.create(street=s1, source_pt=j1, slot=0)
a3.print_growing_pairs()

print '\nthe soliton dashes are \n{}'.format([d.label for d in a3.dashes])
print 'growth restrictions are {}\n'.format([d.growth_restriction for d in a3.dashes])

a3.dashes[0].print_endpoints()
print 'starting point {}'.format(a3.dashes[0].starting_point)
print 'ending point {}\n'.format(a3.dashes[0].ending_point)

a3.dashes[1].print_endpoints()
print 'starting point {}'.format(a3.dashes[1].starting_point)
print 'ending point {}'.format(a3.dashes[1].ending_point)

print '\nsoliton starting point {}'.format(a3.starting_point)
print 'soliton ending point {}'.format(a3.ending_point)

cls_3 = growing_clusters(a3)
print '\nThe growing clusters:\n{}'.format(cls_3)
print '\nNow I grow soliton a3 at {}'.format(cls_3[0][0][0].end_point.label)
new_sols = j_type_3_growth(a3, cls_3[0])
print '\nThere are {} new solitons after growth'.format(len(new_sols))
a4=new_sols[0]

print '\nthe soliton dashes are \n{}'.format([d.label for d in a4.dashes])
print 'with growth restrictions {}\n'.format([d.growth_restriction for d in a4.dashes])

a4.dashes[0].print_endpoints()
print 'starting point {}'.format(a4.dashes[0].starting_point)
print 'ending point {}\n'.format(a4.dashes[0].ending_point)

a4.dashes[1].print_endpoints()
print 'starting point {}'.format(a4.dashes[1].starting_point)
print 'ending point {}\n'.format(a4.dashes[1].ending_point)

a4.dashes[2].print_endpoints()
print 'starting point {}'.format(a4.dashes[2].starting_point)
print 'ending point {}\n'.format(a4.dashes[2].ending_point)

print '\nThe soliton is now complete: {}'.format(a4.is_complete)


print '\n\n-------------------------------------------------------'
print '\nKeep on growing automatically, with one more step'

cls_4 = growing_clusters(a4)
print '\nThe growing clusters:\n{}'.format(cls_4)
print '\nNow I grow soliton a4 at {}'.format(cls_4[0][0][0].end_point.label)
new_sols_5 = bp_type_1_growth(a4, cls_4[0])
print '\nThere are {} new solitons after growth'.format(len(new_sols))
a5=new_sols_5[0]

print '\nthe soliton dashes are \n{}'.format([d.label for d in a4.dashes])
print 'with growth restrictions {}\n'.format([d.growth_restriction for d in a4.dashes])

a4.dashes[0].print_endpoints()
print 'starting point {}'.format(a4.dashes[0].starting_point)
print 'ending point {}\n'.format(a4.dashes[0].ending_point)

a4.dashes[1].print_endpoints()
print 'starting point {}'.format(a4.dashes[1].starting_point)
print 'ending point {}\n'.format(a4.dashes[1].ending_point)

a4.dashes[2].print_endpoints()
print 'starting point {}'.format(a4.dashes[2].starting_point)
print 'ending point {}\n'.format(a4.dashes[2].ending_point)

print '\nThe soliton is now complete: {}'.format(a4.is_complete)




# print (
# 	'The soliton starting point equals the dash starting point: {}'
# 	.format(a4.starting_point==a4.complete_dash.starting_point)
# )
# print (
# 	'The soliton ending point equals the dash ending point: {}'
# 	.format(a4.ending_point==a4.complete_dash.ending_point)
# )
