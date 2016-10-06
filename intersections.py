import sympy
from sympy import Rational

"""
The branch_point_intersection dictionary encodes the intersections
of chords across a branch point.
Let slot 0 be SW, slot 1 be SE, slot 2 be N.
Then we choose the out-going lanes of 2-way streets to be of types
s_0 : ij, s_1 : ji, s_2 : ji
with branch cut running between s_1 and s_2.
So we have chords of type i running from slot 0 to 1, 
from 0 to 2, etc. Likewise we have chords of type j running from 2 to 0
and from 1 to 0. 
There are also mixed-type chords, from 1 to 2, and from 2 to 1, 
but also from 0 to 0 and from 1 to 1 and 2 to 2.
For example, we have intersections 
< [0, 0], [1, 1] > = -1
or also
< [0, 1], [0, 2] > = 1/2
So we encode these into a dictionary:
{'0011': -1, '0102' : 1/2, ...}
and so on.
We only specify < [0, 1], [0, 2] > and not its negative < [0, 2], [0, 1] > 
to minimize errors. That will be taken care of by a function.
We also don't specify zero intersections, those will be understood.
"""
branch_point_intersections_british = {
    '0121': Rational(-1, 2),
    '0102': Rational(1, 2),
    '0212': Rational(1, 2),
    '0100': Rational(1, 2),
    '0200': Rational(-1, 2),
    '0111': Rational(-1, 2),
    '2111': Rational(1, 2),
    '0222': Rational(1, 2),
    '1222': Rational(-1, 2),
    '1020': Rational(1, 2),
    '1012': Rational(-1, 2),
    '2021': Rational(1, 2),
    '2000': Rational(-1, 2),
    '1000': Rational(1, 2),
    '2022': Rational(1, 2),
    '2122': Rational(-1, 2),
    '1211': Rational(1, 2),
    '1011': Rational(-1, 2),
    '0011': -1,
    '1122': -1,
    '2200': -1,
}

"""
The joint_intersection dictionary encodes the intersections
of chords across a joint.
Let slot 0 be S, slot 1 be SE, slot 2 be NE, slot 3 be N,
slot 4 be NW, slot 5 be SW.
Then we choose the out-going lanes of 2-way streets to be of types
s_0 : kj, s_1 : ij, s_2 : ik, s_3: jk, s_4: ji, s_5: ki
as in figure 46 of GMN5 (Spectral Networks).
So we have chords of type i running from slot 1 to 5, 
from 1 to 4, from 2 to 4 and from 2 to 5.
Likewise we have chords of type j running from 3 to 1, etc,
and so on for type k chords.
Of course, type i chords only intersect other type i chords, and 
similarly for j and k types.
For example, we have intersections 
< [1, 4], [1, 5] > = 1/2
or also
< [1, 4], [2, 5] > = 1
So we encode these into a dictionary:
{'1415': 1/2, '1425' : 1, ...}
and so on.
We only specify < [1, 4], [1, 5] > and not its negative < [1, 5], [1, 4] > 
to minimize errors. That will be taken care of by a function.
We also don't specify zero intersections, those will be understood.
"""
joint_intersections_british = {
    '1415': Rational(1, 2),
    '1424': Rational(1, 2),
    '1425': 1,
    '1525': Rational(1, 2),
    '2425': Rational(1, 2),
    '3031': Rational(1, 2),
    '3040': Rational(1, 2),
    '3041': 1,
    '3141': Rational(1, 2),
    '4041': Rational(1, 2),
    '5253': Rational(1, 2),
    '5202': Rational(1, 2),
    '5203': 1,
    '5303': Rational(1, 2),
    '0203': Rational(1, 2),
}

"""
As it turns out the intersection dictionaries for the british
and american resolutions are the same with our conventions. 
See the separate notes for the details.
"""
branch_point_intersections_american = branch_point_intersections_british
joint_intersections_american = joint_intersections_british


class DashNode:
    """
    Class which stores nodal data of a dash.
    Given two pieces of the Dash's path, i.e. two lists of the type
    [Street, +/-1]
    it constructs a class which contains information about the node: 
    which branch point or joint is the node between the two streets,
    and which chord the dash traces through the node.

    The chord type is just a pair of slots, like [0,2] for a 
    chord going from slot #0 to slot #2.
    """
    def __init__(self, path_piece_0, path_piece_1):
        self.node = None
        self.street_in = path_piece_0[0]
        self.street_out = path_piece_1[0]

        if path_piece_0[1] == 1:
            # the orientation of the dash along the street_in
            # coincides with that of the street
            node_in = self.street_in.final_point().end_point
            slot_in = self.street_in.final_point().slot
        elif path_piece_0[1] == -1:
            # the orientation of the dash along the street_in
            # is opposite to that of the street
            node_in = self.street_in.initial_point().end_point
            slot_in = self.street_in.initial_point().slot
        else:
            raise ValueError

        if path_piece_1[1] == 1:
            # the orientation of the dash along the street_out
            # coincides with that of the street
            node_out = self.street_out.initial_point().end_point
            slot_out = self.street_out.initial_point().slot
        elif path_piece_1[1] == -1:
            # the orientation of the dash along the street_out
            # is opposite to that of the street
            node_out = self.street_out.final_point().end_point
            slot_out = self.street_out.final_point().slot
        else:
            raise ValueError

        # check that node_in coincides with node_out
        if node_in == node_out:
            pass
        else:
            raise Exception

        self.node = node_in
        self.chord = [slot_in, slot_out]

    def print_info(self):
        print 'Node {}, chord {}'.format(
            self.node.label,
            self.chord
        )


def get_dash_nodes(dash, is_closed_soliton=False):
    """
    Given a dash (either one corresponding to an open or a closed path)
    this returns the node data. 
    If the Dash's path is 
    [[street_0, +/-1], [street_1, +/-1], [street_2, +/-1], ...]
    it returns a list of DashNode objects
    [DashNode_0, DashNode_1, ...]
    corresponding to the nodes between street_0 and street_1, 
    then between street_1 and street_2, and so on.

    For dashes of closed solitons, must add by hand the node 
    between the last dash and the first one.
    """

    nodes = []
    for i in range(len(dash.path) - 1):
        nodes.append(DashNode(dash.path[i], dash.path[i + 1]))

    if is_closed_soliton is True:
        nodes.append(DashNode(dash.path[-1], dash.path[0]))    

    return nodes


def compute_self_intersections(dash, is_closed_soliton=None, resolution=None):
    """
    Computes the self intersections of a dash. 
    The computation depends on the resolution chosen for the network,
    as this determines how the dashes self-intersect at joints 
    and branch points.
    """
    self_int = 0

    if is_closed_soliton is True or is_closed_soliton is False:
        pass
    else:
        raise ValueError('Must specify whether it is a closed soliton or not')

    dash.nodes = get_dash_nodes(dash, is_closed_soliton=is_closed_soliton)

    # build a list of branch points and joints 
    # through which the dash develops.
    distinct_nodes = []
    for d_n in dash.nodes:
        if d_n.node not in distinct_nodes:
            distinct_nodes.append(d_n.node)

    # compute the intersections for each branch point or joint
    for n in distinct_nodes:
        # n is a branch point or joint, while 
        # dn below is a dash node object
        chords = [dn.chord for dn in dash.nodes if dn.node == n]
        # compute intersections of each chord with subsequent ones
        for i_c, c in enumerate(chords):
            self_int += intersections_with_later_chords(
                i_c, chords, n, resolution=resolution,
            )

    return self_int


def compute_self_intersections_2(dash, is_closed_soliton=None, resolution=None):
    """
    Computes the self intersections of a dash. 
    The computation depends on the resolution chosen for the network,
    as this determines how the dashes self-intersect at joints 
    and branch points.
    """
    self_int = 0

    if is_closed_soliton is True or is_closed_soliton is False:
        pass
    else:
        raise ValueError('Must specify whether it is a closed soliton or not')

    dash.nodes = get_dash_nodes(dash, is_closed_soliton=is_closed_soliton)

    # build a list of branch points and joints 
    # through which the dash develops.
    distinct_nodes = []
    for d_n in dash.nodes:
        if d_n.node not in distinct_nodes:
            distinct_nodes.append(d_n.node)

    print 'distinct nodes:'
    print distinct_nodes

    # compute the intersections for each branch point or joint
    for n in distinct_nodes:
        print '\nat node {}'.format([n.label])
        # n is a branch point or joint, while 
        # dn below is a dash node object
        chords = [dn.chord for dn in dash.nodes if dn.node == n]
        print 'chords are {}'.format(chords)
        # compute intersections of each chord with subsequent ones
        for i_c, c in enumerate(chords):
            print 'chord {} has intersection {} with later chords'.format(
                i_c, intersections_with_later_chords(
                i_c, chords, n, resolution=resolution,)
                )
            self_int += intersections_with_later_chords(
                i_c, chords, n, resolution=resolution,
            )

    return self_int


def intersections_with_later_chords(c_index, chords, node, resolution=None):
    tot_int = 0

    # the current chord
    c_c = chords[c_index]
    later_chords = chords[c_index + 1:]

    # the intersection dictionary
    if resolution == 'british':
        if node.__class__.__name__ == 'BranchPoint':
            int_dict = branch_point_intersections_british
        elif node.__class__.__name__ == 'Joint':
            int_dict = joint_intersections_british
        else:
            raise Exception
    elif resolution == 'american':
        if node.__class__.__name__ == 'BranchPoint':
            int_dict = branch_point_intersections_american
        elif node.__class__.__name__ == 'Joint':
            int_dict = joint_intersections_american
        else:
            raise Exception
    else:
        raise ValueError('Unknown resolution type: {}'.format(resolution))

    for l_c in later_chords:
        intersection_string = (
            str(c_c[0]) + str(c_c[1]) + str(l_c[0]) + str(l_c[1]) 
        )
        opposite_intersection_string = (
            str(l_c[0]) + str(l_c[1]) + str(c_c[0]) + str(c_c[1]) 
        )
        # see if the intersection string appears in the ones 
        # defined in the dictionaries, or whether the opposite string does.
        if intersection_string in int_dict.keys():
            tot_int += int_dict[intersection_string]
        elif opposite_intersection_string in int_dict.keys():
            tot_int += (-1) * int_dict[opposite_intersection_string]
    
    return tot_int
