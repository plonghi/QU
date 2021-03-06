### This module handles operations with the mapping class group on BPS graphs

import sys, itertools, re, numpy, os, datetime, matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mcsn import MCSN, Street, Joint, BranchPoint
from decimal import Decimal

from config import MCSNConfig
from cluster import Seed

# an auxiliary parameter for checking whether, modulo [1, 1]
# two vectors in the plane are the same
EPSILON = 0.0001

SEED_WARNING_MESSAGE = False


class BPSgraph(MCSN):
    """
    A class for BPS graphs.
    The attribute 'mutable_edges' collects (labels of) edges
    which can be flipped.
    The attribute 'mutable_faces' collects (labels of) faces
    on which a cootie move can be performed.

    'e_coordinates' is a dictionary related to edges; 
    it contains coordinates of edges

    'mutable_edges' contains both the I-webs that can be flipped, 
    and the middle-edges of H-webs, which can be 'H-flipped'.

    'faces' contains data that instructs a BPS graph on how to label 
    its own faces (for compatibility thoughout graph moves)
    if its value is None, it lets the graph label its faces automatically.

    In fact, remove altogether the info about positions of BP and JOINTS
    just retain info about edges, and extract from it the infor about 
    BP and Joints. Need to make checks of compatibility, 
    and to allow for integer shifts aroudn the torus.
    """
    def __init__(
        self, 
        label='no_label', 
        branch_points={}, 
        streets=[], 
        joints={}, 
        homology_classes={},
        faces=None,
        bp_coordinates={},
        j_coordinates={},
        e_coordinates={},
        seed=None,
    ):
        MCSN.__init__(
            self,
            label, 
            branch_points, 
            streets, 
            joints, 
            homology_classes,
            quiet=True
        )
        self.bp_coordinates = bp_coordinates
        self.j_coordinates = j_coordinates
        if e_coordinates is not None:
            self.e_coordinates = e_coordinates
        else:
            self.determine_e_coordinates_automatically()
        self.faces = {}
        self.mutable_faces = []
        self.mutable_edges = []
        # keep a check on whether all coordinates are transformed OK after 
        # a flip, for nodes and edges together.
        self.wrong_coordinates = False
        self.check_edge_node_coordinates()
        self.determine_faces(faces)
        for f in self.faces.values():
            f.determine_neighbors(self)
        self.tessellation_data = {
            f.label : f.neighbors for f in self.faces.values()
        }
        self.determine_mutable_elements()
        if seed is None:
            self.seed = Seed(
                self.compute_quiver(), 
                sorted(homology_classes.keys()),
            )
        else:
            self.seed = seed
    

    def determine_faces(self, face_data=None):
        # an array of all initial points
        # for building face paths (these are instances of 
        # StreetEndPoint)
        initial_points = []
        for s in self.streets.values():
            initial_points += s.endpoints

        # record all face paths while they are build to avoid
        # building same face twice from different initial points
        all_face_paths = []

        # build all face paths: if there is no face_data specified,
        # we just label faces arbitrarily.
        # Otherwise we pay attention to match the face labels 
        # with those of the data
        i_f = 0 # a counter
        for ep in initial_points:
            if ep in all_face_paths:
                continue
            else:
                if face_data is None:
                    f_label = 'f_'+str(i_f)
                else:
                    f_label = determine_face_label(ep, self, face_data)
                new_face = Face(ep, label=f_label)
                self.faces[f_label] = new_face
                all_face_paths += new_face.full_sequence
                i_f += 1
        ### DEBUG
        # print '\n\nGRPH FACE DATA'
        # print {f_k : f_v.street_sequence for f_k, f_v in self.faces.iteritems()}
        ###

    def print_face_info(self):
        print '\nThere are {} faces, they are:'.format(len(self.faces))
        for f_k, f_v in self.faces.iteritems():
            print '\nFace {}.'.format(f_k)
            f_v.print_info()

    def print_mutable_info(self):
        print '\nThe are {} mutable edges, they are: {}'.format(
            len(self.mutable_edges), self.mutable_edges
        )
        print 'The are {} mutable faces, they are: {}\n'.format(
            len(self.mutable_faces), self.mutable_faces
        )

    def determine_mutable_elements(self):
        mutable_edges = []
        for s_k, s_v in self.streets.iteritems():
            if (
                can_be_flipped(s_v) is True or
                can_be_H_flipped(s_v) is True
            ):
                mutable_edges.append(s_k)

        mutable_faces = []
        for f_k, f_v in self.faces.iteritems():
            if can_cootie(self, f_v) is True:
                mutable_faces.append(f_k)

        self.mutable_edges = mutable_edges
        self.mutable_faces = mutable_faces

    def check_edge_node_coordinates(self):
        for e_k, e_v in self.streets.iteritems():
            xy_0 = node_coordinates(e_v.endpoints[0].end_point.label, self)
            xy_1 = node_coordinates(e_v.endpoints[1].end_point.label, self)
            edge_xy_0 = self.e_coordinates[e_k][0]
            edge_xy_1 = self.e_coordinates[e_k][1]
            if (
                are_within_range(xy_0, edge_xy_0, EPSILON) is True and
                are_within_range(xy_1, edge_xy_1, EPSILON) is True
                ) or (
                are_within_range(xy_0, edge_xy_1, EPSILON) is True and
                are_within_range(xy_1, edge_xy_0, EPSILON) is True
            ):
                pass
            else:
                ### USEFUL for DEBUGGING: just comment out the exception
                self.wrong_coordinates = True
                print [xy_0, xy_1]
                print [edge_xy_0, edge_xy_1]
                print are_within_range(xy_0, edge_xy_0, EPSILON)
                print are_within_range(xy_1, edge_xy_1, EPSILON) 
                
                raise Exception(
                    'Coordinates of {} seems not to match ' 
                    'those of its endpoints. The edge has {} '
                    'while the endpoints are {}'.format(
                        e_k, [edge_xy_0, edge_xy_1], [xy_0, xy_1])
                )

    def copy(self):
        new_streets = [k for k in self.streets.keys()]
        new_branch_points = {
            b.label : [s.label for s in b.streets] 
            for b in self.branch_points.values()
        }
        new_joints = {
            j.label : [get_label(s) for s in j.streets] 
            for j in self.joints.values()
        }
        new_homology_classes = {
            h.label : [s.label for s in h.streets] 
            for h in self.basis_homology_classes.values()
        }
        new_bp_coordinates = self.bp_coordinates.copy()
        new_j_coordinates = self.j_coordinates.copy()
        new_e_coordinates = self.e_coordinates.copy()
        new_seed = self.seed.copy()
        new_face_data = (
            {f_k : f_v.street_sequence for f_k, f_v in graph.faces.iteritems()}
        )
        return BPSgraph(
            branch_points=new_branch_points, 
            streets=new_streets, 
            joints=new_joints, 
            homology_classes=new_homology_classes,
            faces=new_face_data,
            bp_coordinates=new_bp_coordinates,
            j_coordinates=new_j_coordinates,
            e_coordinates=new_e_coordinates,
            seed=new_seed
        )

    def homology_dictionary(self):
        """
        returns the homology dictionary as was originally used to 
        specify basis homology classes
        """
        return (
            {hc_k : [s.label for s in hc_v.streets] for 
            hc_k, hc_v in self.basis_homology_classes.iteritems()}
        )

    def compute_quiver(self):
        """
        Compute the intersection pairing matrix for the quiver,
        working in the basis of the graph's own homology dictionary.
        We order the dictionary according to lexicographic ordering of keys:
        {'gamma_1' : ..., 'gamma_2' : ...,  etc.}
        Returns a numpy matrix.
        """
        hom_dic = self.homology_dictionary()
        hc = sorted(hom_dic.keys())
        dim = len(hc)
        matrix = numpy.matrix([[0 for i in range(dim)] for j in range(dim)])
        for r in range(dim):
            for c in range(dim):
                matrix[r, c] = intersection_pairing(self, hc[r], hc[c])

        return matrix

    def determine_e_coordinates_automatically(self):
        """
        Only use if you don't care about correct shifts of coordinates
        around cycles of the torus. This brutally gives a set of consistent 
        coordinates.
        """
        self.e_coordinates = {}
        for s_k, s_v in self.streets.iteritems():
            self.e_coordinates[s_k] = (
                [node_coordinates(ep.end_point.label, self) 
                for ep in s_v.endpoints]
            )


class Face():
    """
    A class for faces of BPS graphs
    """
    def __init__(self, endpoint, label=None):
        self.full_sequence = []
        self.node_sequence = []
        self.street_sequence = []
        self.neighbors = []
        self.face_type = []
        self.label = label

        [
            self.full_sequence, 
            self.node_sequence, 
            self.street_sequence
        ] = build_face(endpoint)

        # the face type is a list of letters like ['b', 'b', 'j',...]
        # which records the types of nodes at the feca boundary
        # (b is for branch point, j for joint)
        for el in self.full_sequence:
            if el.__class__.__name__ == 'StreetEndpoint':
                if el.end_point.__class__.__name__ == 'BranchPoint':
                    self.face_type += ['b']
                elif el.end_point.__class__.__name__ == 'Joint':
                    self.face_type += ['j']
                else:
                    pass

    def determine_neighbors(self, graph):
        self.neighbors = face_neighbors(graph, self.label)

    def print_info(self):
        print 'Edges of the face\n{}'.format(self.street_sequence)
        print 'Nodes of the face\n{}'.format(self.node_sequence)
        print 'Face type: {}'.format(self.face_type)


class MoveSequence():
    """
    A container class for the sequences of flips and cootie moves 
    that we find relating a graph to itself.

    - The sequence of moves is given as a list of elements of the graph 
    namely streets for flips and faces for cooties. These are stored in the
    attribute 'moves_sequence'
    - To a given move, there can correspond several permutations of streets.
    These are stored in an attribute 'edge_permutations'
    - Likewise, the permutations of homology classes that result from these 
    edge permutations are stored in the attribute 'homology_permutations'.
    - The quiver mutations are also stored here into the attribute 
    'quiver_mutations'

    """
    def __init__(
        self, moves=None, edge_permutations=None, 
        quiver_mutations=None, homology_permutations=None,
    ):
        self.moves_sequence = moves
        self.edge_permutations = edge_permutations
        self.quiver_mutations = quiver_mutations
        self.homology_permutations = homology_permutations

    def print_info(self):
        print '\nSequence of graph moves\n{}'.format(self.moves_sequence)
        print '\nEdge permutations\n{}'.format(self.edge_permutations)
        print '\nQuiver node mutations\n{}'.format(self.quiver_mutations)
        print '\nHomology class permutations\n{}\n\n'.format(
            self.homology_permutations
        )


def intersection_pairing(graph, gamma_1, gamma_2):
    """
    Takes a graph, and two labels of homology classes.
    Returns their intersection pairing <gamma_1, gamma_2>
    """
    if gamma_1 == gamma_2:
        return 0

    hom_dic = graph.homology_dictionary()

    streets_1 = [graph.streets[s] for s in hom_dic[gamma_1]]
    streets_2 = [graph.streets[s] for s in hom_dic[gamma_2]]

    pairing = 0
    for s1 in streets_1:
        for s2 in streets_2:
            for ep1 in s1.endpoints:
                for ep2 in s2.endpoints:
                    if (
                        ep1.end_point.label == ep2.end_point.label and
                        ep1.end_point.__class__.__name__ == 'BranchPoint'
                    ): 
                        if (ep1.slot % 3) == ((ep2.slot + 1) % 3):
                            # street of gamma_1 sits CCW of street of gamma_2
                            # then this contributes +1
                            pairing += 1
                        elif (ep1.slot % 3) == ((ep2.slot - 1) % 3):
                            # street of gamma_1 sits CW of street of gamma_2
                            # then this contributes -1
                            pairing += -1
                        else:
                            raise Exception(
                                '{} and {} cannot belong to two different '
                                'homology classes and end on the same slot '
                                'of {}'.format(
                                    s1.label, s2.label, ep1.end_point.label
                                )
                            )
    return pairing


def build_face(end_pt):
    """
    build a face associated with a street's endpoint
    where 'endpoint' really means an instance of the class
    'StreetEndpoint'
    By convention we start from the street to whom the end point 
    belongs, we go INTO the end point, and then we follow 
    the face counterclockwise.
    We keep track of all nodes and streets along the face boundary.
    """
    node_in = end_pt.end_point
    slot_in = end_pt.slot
    street_in = end_pt.street
    initial_point = [node_in, slot_in]
    # record both nodes and streets, in ordered sequences
    face_seq = [street_in, end_pt]
    face_seq_nodes = [node_in.label]
    face_seq_streets = [street_in.label]

    closed = False

    while closed is False:
        if node_in.__class__.__name__ == 'BranchPoint':
            n_slots = 3
            slot_out = (slot_in - 1) % n_slots
        elif node_in.__class__.__name__ == 'Joint':
            n_slots = 6
            # at a joint, it can be that the next slot is 
            # empty, for example in a tri-valent joint
            # then take the next one.
            # Avoid 0 (that would be the 'in' slot), start from 1
            for k in range(1, 5):
                if node_in.streets[slot_in - k] is not None:
                    slot_out = (slot_in - k) % n_slots
                    break

        street_out = node_in.streets[slot_out]

        # determine which of the new street's endpoints will be the next one
        [new_ep0, new_ep1] = street_out.endpoints
        if new_ep0.end_point == node_in and new_ep0.slot == slot_out:
            ep_out = new_ep1
            node_out = new_ep1.end_point
            slot_out = new_ep1.slot
        elif new_ep1.end_point == node_in and new_ep1.slot == slot_out:
            ep_out = new_ep0
            node_out = new_ep0.end_point
            slot_out = new_ep0.slot
        else:
            raise Exception(
                'Cannot determine next step of the face path. '
                'Maybe an open face? (irregular singularity)'
            )
        if [node_out, slot_out] == initial_point:
            closed = True
            break
        else:
            # record progress
            face_seq += [street_out, ep_out]
            face_seq_nodes += [node_out.label]
            face_seq_streets += [street_out.label]

            # update 'in' data for next loop
            node_in = node_out
            slot_in = slot_out

    return [face_seq, face_seq_nodes, face_seq_streets]


def determine_face_label(endpoint, graph, face_data):
    """
    Determine how to label a certain face of a BPS graph,
    based on the given 'face_data'.

    The face is identified by 'endpoint': this is really an instance 
    of the class 'StreetEndpoint'
    By convention the face is built by starting from the street 
    to whom the end point belongs, going INTO the end point, and following
    the face counterclockwise.
    We keep track of all nodes and streets along the face boundary.

    The variable 'face_data' is a dictionary like
    {..., 'f_i' : [p_j, n_i, p_k,...],...}
    with face labels as keys and the 'STREET' sequence of edges  
    bounding the face as values.
    The ordering of the edges matters, up to cyclic permutations.
    """

    # the last entry is the street sequence
    [face_seq, face_seq_nodes, face_seq_streets] = build_face(endpoint)
    found_match = False
    for f_k, f_v in face_data.iteritems():
        if are_same_cycle(f_v, face_seq_streets):
            found_match = True
            face_label = f_k
            break

    if found_match is True:
        return face_label
    else:
        raise Exception('Cannot determine face label, from given face data.')
 

def shift_number_into_interval_01(x):
    # NOTE: crucial to use numpy.fmod, because of issues described e.g. here
    # http://stackoverflow.com/questions/14763722/python-modulo-on-floats

    # this returns a number between -1.0 and +1.0
    ans = numpy.fmod(x, 1.0)
    # Then, we shift it between 0.0 and 1.0, if necessary
    if ans < 0.0:
        ans = ans + 1.0
    return ans

def are_within_range(xy_A, xy_B, eps):
    """
    Given xy_A = [x_0, y_0] and similarly xy_B,
    determine whether, taken modulo 1, their distance is less than epsilon.
    For example, if xy_A = [0,0] and xy_B = [1,1] then they are 
    effectively coincident. (The torus is square with side 1).
    """
    xy_A_mod = [shift_number_into_interval_01(z) for z in xy_A]
    xy_B_mod = [shift_number_into_interval_01(z) for z in xy_B]

    # now all the numbers are shifted between [0.0 , 1.0]
    # however there is still one problem: the above shift would take
    # -0.00001 --> 0.999999
    # but also 
    # 1.0 --> 0.0
    # so the two would look far away.
    # So we compute the distance and we take it modulo 1 and finally 
    # take its absolute value. If this is small enough, 
    # the two points are close.
    # MORE PRECISELY, we must take the difference of x-coordinates mod 1
    # and the difference of y-coordinates mod 1, SEPARATELY!
    # Then we build a 'modded' distance vector with those remainders,
    # and check if its norm is smaller than epsilon.
    delta_x = min(
        abs(shift_number_into_interval_01(xy_A_mod[0] - xy_B_mod[0])),
        abs(shift_number_into_interval_01(xy_B_mod[0] - xy_A_mod[0]))
    )
    delta_y = min(
        abs(shift_number_into_interval_01(xy_A_mod[1] - xy_B_mod[1])),
        abs(shift_number_into_interval_01(xy_B_mod[1] - xy_A_mod[1]))
    )
    d = numpy.linalg.norm(numpy.array([delta_x, delta_y]))
    if d < eps:
        return True
    else:
        return False


def other_edge_coord(edge_coordinates, xy):
    """
    Given a list of two coordinates [[x_0, y_0], [x_1, y_1]] and
    xy = [x, y], it returns [x_0, y_0] if [x, y] = [x_1, y_1], and
    viceversa. Everything up to shifts by [n, m] (square torus)
    """
    if (
        are_within_range(xy, edge_coordinates[0], EPSILON) is True and 
        are_within_range(xy, edge_coordinates[1], EPSILON) is False
    ):
        return edge_coordinates[1]
    elif (
        are_within_range(xy, edge_coordinates[0], EPSILON) is False and 
        are_within_range(xy, edge_coordinates[1], EPSILON) is True
    ):
        return edge_coordinates[0]
    elif (
        are_within_range(xy, edge_coordinates[0], EPSILON) is True and 
        are_within_range(xy, edge_coordinates[1], EPSILON) is True
    ):
        raise Exception()
        
    else:
        print 'xy = {}'.format(xy)
        print 'edge_coordinates = {}'.format(edge_coordinates)
        raise Exception()


def determine_edge_shift(edge_coordinates, xy):
    """
    Given a list of two coordinates [[x_0, y_0], [x_1, y_1]] and
    xy = [x, y], it computes whether [x, y] = [x_0, y_0], or the opposite.
    Then, it computes the shift of fundamental domain e.g. 
    [n, m] = [x_0 - x, y_0 - y]
    """
    if (
        are_within_range(xy, edge_coordinates[0], EPSILON) is True and 
        are_within_range(xy, edge_coordinates[1], EPSILON) is False
    ):
        return list( numpy.array(edge_coordinates[0]) - numpy.array(xy) ) 
    elif (
        are_within_range(xy, edge_coordinates[0], EPSILON) is False and 
        are_within_range(xy, edge_coordinates[1], EPSILON) is True
    ):
        return list( numpy.array(edge_coordinates[1]) - numpy.array(xy) )
    elif (
        are_within_range(xy, edge_coordinates[0], EPSILON) is True and 
        are_within_range(xy, edge_coordinates[1], EPSILON) is True
    ):
        raise Exception(
            'The endpoints of some edges appear to coincide. '
            'Try changing the initial positions of edges by a bit.'
        )
    else:
        print 'xy = {}'.format(xy)
        print 'edge_coordinates = {}'.format(edge_coordinates)
        raise Exception()


def flip_edge(graph, edge):
    """
    Perform a flip move on a BPS graph, 
    at the edge whose label is given as an argument.
    NOTE: this function can also be used for an H-flip,
    that is, it can be used to flip and edge within an H-web,
    which ends on joints at both endpoints.
    """
    if edge not in graph.streets.keys():
        raise Exception('This edge is not part of the graph.')
    else:
        # the street
        p = graph.streets[edge]

    ep0, ep1 = p.endpoints

    # check that we can perform a flip on this edge: must be an I-web
    if can_be_flipped(p) is False and can_be_H_flipped(p) is False:
        raise Exception('Cannot flip / H-flip on edge {}.'.format(edge))

    # check whether it's an I-web or an H-web
    if (
        ep0.end_point.__class__.__name__ == 'BranchPoint' and 
        ep1.end_point.__class__.__name__ == 'BranchPoint'
    ):
        edge_type = 'I'
        slot_shift = 1
    elif (
        ep0.end_point.__class__.__name__ == 'Joint' and 
        ep1.end_point.__class__.__name__ == 'Joint'
    ):
        edge_type = 'H'
        slot_shift = 2
        ### DEBUG BEGINS
        # print '\n\n*** H-flipping an edge***\n\n'
        ### DEBUG ENDS
    else:
        raise Exception('Cannot flip this edge: its ends are {}'.format(
            [
                ep0.end_point.__class__.__name__, 
                ep1.end_point.__class__.__name__
            ]
        ))


    new_streets = [k for k in graph.streets.keys()]
    new_branch_points = {
        b.label : [s.label for s in b.streets] 
        for b in graph.branch_points.values()
    }
    new_joints = {
        j.label : [get_label(s) for s in j.streets] 
        for j in graph.joints.values()
    }
    new_homology_classes = {
        h.label : [s.label for s in h.streets] 
        for h in graph.basis_homology_classes.values()
    }

    l_0 = ep0.end_point.label
    slot_0 = ep0.slot
    
    l_1 = ep1.end_point.label
    slot_1 = ep1.slot

    if edge_type == 'I':
        old_node_xy_0 = graph.bp_coordinates[l_0]
        old_node_xy_1 = graph.bp_coordinates[l_1]
        del new_branch_points[l_0]
        del new_branch_points[l_1]
    elif edge_type == 'H':
        old_node_xy_0 = graph.j_coordinates[l_0]
        old_node_xy_1 = graph.j_coordinates[l_1]
        del new_joints[l_0]
        del new_joints[l_1]

    # now I label the streets meeting at the two endpoints as follows:
    # at endpoint ep0 I have in ccw order: {p, s00, s01}
    # at endpoint ep1 I have in ccw order: {p, s10, s11}
    # then I must make new end points (branch points for I-web, or 
    # joints for an H-web) where I have, in ccw order
    # {p, s11, s00} at l_0 and {p, s01, s10} at l_1
    # let's get the labels fo these streets then:
    s00 = ep0.end_point.streets[
        (slot_0 + 1 * slot_shift) % (3 * slot_shift)
    ].label
    s01 = ep0.end_point.streets[
        (slot_0 + 2 * slot_shift) % (3 * slot_shift)
    ].label
    s10 = ep1.end_point.streets[
        (slot_1 + 1 * slot_shift) % (3 * slot_shift)
    ].label
    s11 = ep1.end_point.streets[
        (slot_1 + 2 * slot_shift) % (3 * slot_shift)
    ].label

    # now add the new keys to the dictionary of branch points / joints
    if edge_type == 'I':
        new_branch_points[l_0] = [p.label, s11, s00]
        new_branch_points[l_1] = [p.label, s01, s10]
    elif edge_type == 'H':
        new_joints[l_0] = [p.label, None, s11, None, s00, None]
        new_joints[l_1] = [p.label, None, s01, None, s10, None]

    # Now update the coordinates for endpoints of edges of the graph.
    # Recall that these are the same as the coordinates of the endpoints,
    # defined above as old_node_xy_0/1, but only up to
    # the ambiguity by integer shifts around cycles of the torus.
    # Therefore they must be computed differently, using directy the 
    # information stored about the edges.
    # First of all, note that we don't know if the endpoints (nodes) 
    # of the edge that we're flipping match with those of its coordinate data, 
    # so we first match them.
    if (
        are_within_range(
            old_node_xy_0, graph.e_coordinates[edge][0], EPSILON
        ) and are_within_range(
            old_node_xy_1, graph.e_coordinates[edge][1], EPSILON
        )
    ):
        old_edge_xy_0 = graph.e_coordinates[edge][0]
        old_edge_xy_1 = graph.e_coordinates[edge][1]
    elif (
        are_within_range(
            old_node_xy_0, graph.e_coordinates[edge][1], EPSILON
        ) and are_within_range(
            old_node_xy_1, graph.e_coordinates[edge][0], EPSILON
        )
    ):
        old_edge_xy_0 = graph.e_coordinates[edge][1]
        old_edge_xy_1 = graph.e_coordinates[edge][0]
    else: 
        print 'old coordinate data of nodes {}'.format(old_node_xy_0)
        print 'old coordinate data of edges {}'.format(graph.e_coordinates[edge])
        raise Exception()
    # now at this point we are guaranteed that old_edge_xy_0 (resp 1)
    # correspond to the coordinate of l_0 (resp 1), before flipping.
    edge_center = (numpy.array(old_edge_xy_0) + numpy.array(old_edge_xy_1)) / 2
    edge_sep = numpy.array(old_edge_xy_0) - numpy.array(old_edge_xy_1) 
    rot = numpy.array([[0, 1], [-1, 0]])
    new_edge_xy_0 = list(edge_center + (rot.dot(edge_sep) / 2))
    new_edge_xy_1 = list(edge_center - (rot.dot(edge_sep) / 2))
    new_e_coordinates = graph.e_coordinates.copy()
    new_e_coordinates[edge] = [new_edge_xy_0, new_edge_xy_1]
    # Then we update the cooridnates of those four edges that ended on the 
    # flipping edge. But we must be careful: each of them may belong 
    # to a different 'sector' of the covering!
    # So we must determine the suitable shifts for their own endpoints,
    # by comparing the position of the old branch point with the position 
    # of the old edge's endpoint, for each edge.
    s00_shift = determine_edge_shift(
        graph.e_coordinates[s00], old_edge_xy_0
    )
    s00_other_endpt = other_edge_coord(
        graph.e_coordinates[s00], old_edge_xy_0
    )
    new_e_coordinates[s00] = [
        list(numpy.array(new_edge_xy_0) + numpy.array(s00_shift)), 
        s00_other_endpt
    ]

    s01_shift = determine_edge_shift(
        graph.e_coordinates[s01], old_edge_xy_0
    )
    s01_other_endpt = other_edge_coord(
        graph.e_coordinates[s01], old_edge_xy_0
    )
    new_e_coordinates[s01] = [
        list(numpy.array(new_edge_xy_1) + numpy.array(s01_shift)), 
        s01_other_endpt
    ]

    s10_shift = determine_edge_shift(
        graph.e_coordinates[s10], old_edge_xy_1
    )
    s10_other_endpt = other_edge_coord(
        graph.e_coordinates[s10], old_edge_xy_1
    )
    new_e_coordinates[s10] = [
        list(numpy.array(new_edge_xy_1) + numpy.array(s10_shift)), 
        s10_other_endpt
    ]

    s11_shift = determine_edge_shift(
        graph.e_coordinates[s11], old_edge_xy_1
    )
    s11_other_endpt = other_edge_coord(
        graph.e_coordinates[s11], old_edge_xy_1
    )
    new_e_coordinates[s11] = [
        list(numpy.array(new_edge_xy_0) + numpy.array(s11_shift)), 
        s11_other_endpt
    ]

    # Also update positions of the branch points / joints, now using
    # the separation induced by the edge coordinates, though!
    # TODO: get rid of this step, by automatically computing positions
    # of branch points and joints from the edge coordinates.
    # TODO: DECIDE whether keep all nodes within the unit square, 
    # or let them go around. Probably it's ok to have either choice.
    new_node_xy_0 = new_edge_xy_0
    new_node_xy_1 = new_edge_xy_1

    if edge_type == 'I':
        new_bp_coordinates = graph.bp_coordinates.copy()
        new_bp_coordinates[l_0] = new_node_xy_0
        new_bp_coordinates[l_1] = new_node_xy_1
        new_j_coordinates = graph.j_coordinates.copy()
    elif edge_type == 'H':
        new_j_coordinates = graph.j_coordinates.copy()
        new_j_coordinates[l_0] = new_node_xy_0
        new_j_coordinates[l_1] = new_node_xy_1
        new_bp_coordinates = graph.bp_coordinates.copy()

    # #### DEBUGGING BEGINS ####   
    # print '\n\nflipped edge: {}'.format(edge)
    # print 'old edge coordinates = \n{}\n{}'.format(old_edge_xy_0, old_edge_xy_1)
    # print 'new edge coordinates = \n{}\n{}'.format(new_edge_xy_0, new_edge_xy_1)
    # print 'new l_0: {}= {} at {}'.format(l_0, new_branch_points[l_0], new_bp_coordinates[l_0])
    # print 'new l_1: {}= {} at {}'.format(l_1, new_branch_points[l_1], new_bp_coordinates[l_1])
    # print 'other edges involved: {}'.format([s00, s01, s10, s11])
    # print 'these are the shifts: {}'.format([s00_shift, s01_shift, s10_shift, s11_shift])
    # for s_i in [s00, s01, s10, s11]:
    #     print '{} goes from \n{}\nto\n{} '.format(
    #         s_i, graph.e_coordinates[s_i], new_e_coordinates[s_i]
    #     )
    # #### DEBUGGING ENDS ####

    # must mutate the seed if the flip was done on an I-web
    if edge_type == 'I':
        # determine homology class of the flipped edge
        for hc_k, hc_v in graph.homology_dictionary().iteritems():
            if len(hc_v) == 1 and hc_v[0] == edge:
                gamma = hc_k
                break
        
        new_seed = graph.seed.copy()
        # temporarily disable
        global SEED_WARNING_MESSAGE
        if SEED_WARNING_MESSAGE is False:
            print 'WARNING: disabled seed mutation'
            SEED_WARNING_MESSAGE = True
        # new_seed.mutate(gamma)
    else:
        new_seed = graph.seed.copy()

    # must determine the new face structure, to pass the correct labeling
    # of faces to the next BPS graph created after the flip
    old_face_full_data = (
        {f_k : f_v.full_sequence for f_k, f_v in graph.faces.iteritems()}
    )
    old_face_data = (
        {f_k : f_v.street_sequence for f_k, f_v in graph.faces.iteritems()}
    )
    # determine which are the four faces involved (before the flip)
    # call them f_N, f_S, f_W, f_E (North, South, etc.)
    # For example f_W must have [..., s01, s00, ...] in its street 
    # sequence (which goes CCW). 
    # Likewise f_N must have [..., s10, p, s01 ,...], and so on.
    # IMPORTANT: in some systems the North/South/West/East face may be the 
    # SAME FACE! Then, we must check for all these conditions on each face,
    # and update sequentially.
    # ALSO IMPORTANT: some edges may appear more than once within a given 
    # face. In that case one must figure out where to insert/remove the new
    # edge by looking at ALL positions!
    # So, to do this properly we take each face, then go around it by cycling 
    # through the list of its boundary edges.
    # At each step we study the condition whether we are somewhere in the 
    # neighborhood of the flipping edge, and if that's the case, we 
    # modify the list of boundary edges accordingly.
    # IMPORTANT also is to realize that a pair of edges being sequential 
    # like s01 and s00 does not imply that we must insert the edge there.
    # Because those edges may have pairing 2 and the flipping edge
    # may only be attached to one of their common nodes.
    # So, must really use the full face data, and check also that the node
    # between those edges contains the flipping edge.

    ### TO DO: detach this algorithm into a separate function, it's too 
    ### complex to stay here.
    new_face_data = {}
    for f_k, f_v in old_face_full_data.iteritems():
        perimeter_length = len(f_v)

        # only build the sequence of streets for the new face,
        # no info about branch points or joints
        new_face_data[f_k] = []

        # scan through all the elements of the boundary f_v
        for el_i, el in enumerate(f_v):
            if el.__class__.__name__ != 'Street':
                # do not keep track of joints or branch points
                # just record the street sequence
                pass

            else:
                # first look into whether this face is the west-face:
                if (
                    el.label == s01 and 
                    f_v[(el_i + 2) % perimeter_length].label == s00 and (
                        edge in 
                        [s.label for s in 
                        f_v[(el_i + 1) % perimeter_length].end_point.streets 
                        if s is not None]
                    )
                ):
                    # print 'found W!'
                    new_face_data[f_k].append(el.label)
                    new_face_data[f_k].append(edge)
                # then look into whether this face is the east-face:
                elif (
                    el.label == s11 and 
                    f_v[(el_i + 2) % perimeter_length].label == s10 and (
                        edge in 
                        [s.label for s in 
                        f_v[(el_i + 1) % perimeter_length].end_point.streets 
                        if s is not None]
                    )
                ):
                    # print 'found E!'
                    new_face_data[f_k].append(el.label)
                    new_face_data[f_k].append(edge)
                # then look into whether this face is the north-face:
                elif (
                    el.label == edge and 
                    f_v[(el_i + 2) % perimeter_length].label == s11 and 
                    f_v[(el_i - 2) % perimeter_length].label == s00
                ):  
                    # print 'found N!'
                    # don't record this edge then
                    pass
                # then look into whether this face is the south-face:
                elif (
                    el.label == edge and 
                    f_v[(el_i + 2) % perimeter_length].label == s01 and 
                    f_v[(el_i - 2) % perimeter_length].label == s10
                ):  
                    # print 'found S!'
                    # don't record this edge then
                    pass
                
                else:
                    new_face_data[f_k].append(el.label)
            
        ### DEBUG
        # print 'turned face \n{}\ninto\n{}'.format(old_face_data[f_k], new_face_data[f_k])
        ### DEBUG END

    # finally create the new BPS graph
    return BPSgraph(
        branch_points=new_branch_points, 
        streets=new_streets, 
        joints=new_joints, 
        homology_classes=new_homology_classes,
        faces=new_face_data,
        bp_coordinates=new_bp_coordinates,
        j_coordinates=new_j_coordinates,
        e_coordinates=new_e_coordinates,
        seed=new_seed,
    )


def cootie_face(graph, face):
    """
    Perform a cootie move on a BPS graph, 
    at the face whose label is given as an argument
    """
    if face not in graph.faces.keys():
        raise Exception('This face is not part of the graph.')
    else:
        # the street
        f = graph.faces[face]

    # check that we can perform a cootie on this face:
    if can_cootie(graph, f) is False:
        print '\n\nFACE INFO'
        f.print_info()
        raise Exception('Cannot perform a cootie on face {}.'.format(face))

    nodes = f.node_sequence
    streets = f.street_sequence

    new_streets = [k for k in graph.streets.keys()]
    new_branch_points = {
        b.label : [s.label for s in b.streets] 
        for b in graph.branch_points.values()
    }
    new_joints = {
        j.label : [get_label(s) for s in j.streets] 
        for j in graph.joints.values()
    }
    new_homology_classes = {
        h.label : [s.label for s in h.streets] 
        for h in graph.basis_homology_classes.values()
    }
    new_seed = graph.seed.copy()

    # identify the branch points and the joints involved in the cootie move
    # NOTE the convention: j_0 comes AFTER b_0 counterclockwise!
    # also note: at this point we are just collecting the labels 
    # of joints and branch points
    if f.face_type == ['j', 'b', 'j', 'b']:
        b_0 = f.node_sequence[1]
        b_1 = f.node_sequence[3]
        j_0 = f.node_sequence[2]
        j_1 = f.node_sequence[0]
    elif f.face_type == ['b', 'j', 'b', 'j']:
        b_0 = f.node_sequence[0]
        b_1 = f.node_sequence[2]
        j_0 = f.node_sequence[1]
        j_1 = f.node_sequence[3]

    # Identify streets uniquely: for example the external street
    # attached to b_0 will be s_b0, while the boundary street
    # and running b_0 and j_1 will be s_b0_j1, and so on

    b0_streets = [s.label for s in graph.branch_points[b_0].streets]
    b1_streets = [s.label for s in graph.branch_points[b_1].streets]
    j0_streets = [s.label for s in graph.joints[j_0].streets if s is not None]
    j1_streets = [s.label for s in graph.joints[j_1].streets if s is not None]

    # Start with the four external edges
    for s in b0_streets:
        if s not in f.street_sequence:
            s_b0 = s
            break
    for s in b1_streets:
        if s not in f.street_sequence:
            s_b1 = s
            break
    for s in j0_streets:
        if s not in f.street_sequence:
            s_j0 = s
            break
    for s in j1_streets:
        if s not in f.street_sequence:
            s_j1 = s
            break
    
    # Then the four boundary edges
    for s in f.street_sequence:
        if (
            s in b0_streets and 
            s in j0_streets
        ):
            s_b0_j0 = s
        elif (
            s in b0_streets and 
            s in j1_streets
        ):
            s_b0_j1 = s
        elif (
            s in b1_streets and 
            s in j0_streets
        ):
            s_b1_j0 = s
        elif (
            s in b1_streets and 
            s in j1_streets
        ):
            s_b1_j1 = s
        else:
            # ### DEBUG
            # print 'streets on the boundary: {}'.format(f.street_sequence)
            # print 'streets of b0 {}'.format(b0_streets)
            # print 'streets of b1 {}'.format(b1_streets)
            # print 'streets of j0 {}'.format(j0_streets)
            # print 'streets of j1 {}'.format(j1_streets)
            # ### END DEBUG
            raise Exception('{} is not a boundary edge for face {}'.format(
                s, face)
        )

    # Delete old branch points and joints that will be replaced.
    del new_branch_points[b_0]
    del new_branch_points[b_1]
    del new_joints[j_0]
    del new_joints[j_1]

    # now I add a new pair of joints and a new pair of branch points
    # this is actually very easy, because I just need to replace labels
    # for example from
    # {..., b_0 : [p_1, p_2, p_3], ...}
    # {..., j_0 : [p_1, None, p_4, None, p_5, None], ...}
    # to 
    # {..., b_0 : [p_1, p_4, p_5], ...}
    # {..., j_0 : [p_1, None, p_2, None, p_3, None], ...}
    # (Alternatively, I could have used the explicit identification 
    # of all streets obtained above.)
    streets_b_0 = [
        s.label for s in graph.branch_points[b_0].streets if s is not None
    ]
    streets_b_1 = [
        s.label for s in graph.branch_points[b_1].streets if s is not None
    ]
    streets_j_0 = [
        s.label for s in graph.joints[j_0].streets if s is not None
    ]
    streets_j_1 = [
        s.label for s in graph.joints[j_1].streets if s is not None
    ]

    # now add the new keys to the dictionaries of branch points and joints
    new_branch_points[b_0] = streets_j_0
    new_branch_points[b_1] = streets_j_1
    new_joints[j_0] = [
        streets_b_1[0], None, streets_b_1[1], None, streets_b_1[2], None
    ]
    new_joints[j_1] = [
        streets_b_0[0], None, streets_b_0[1], None, streets_b_0[2], None
    ]

    # update coordinates of branch points and joints
    old_xy_b_0 = graph.bp_coordinates[b_0]
    old_xy_b_1 = graph.bp_coordinates[b_1]
    old_xy_j_0 = graph.j_coordinates[j_0]
    old_xy_j_1 = graph.j_coordinates[j_1]

    new_bp_coordinates = graph.bp_coordinates.copy()
    new_j_coordinates = graph.j_coordinates.copy()
    new_e_coordinates = graph.e_coordinates.copy()

    new_bp_coordinates[b_0] = old_xy_j_0
    new_bp_coordinates[b_1] = old_xy_j_1
    new_j_coordinates[j_0] = old_xy_b_1
    new_j_coordinates[j_1] = old_xy_b_0

    # Update the homology classes of the BPS graph.
    # They are denoted by gamma_b0, gamma_b1, gamma_j0, gamma_j1
    # for example gamma_b0 should include s_b0 while 
    # gamma_j0 should include s_j0, s_b0_j0, s_b1_j0
    # and so on.
    
    for h_k, h_v in new_homology_classes.iteritems():
        # To gamma_b0 we must add streets s_b0_j0, s_b0_j1 
        if s_b0 in h_v:
            gamma_b0 = h_k
            new_streets_gamma_b0 = [s for s in h_v] + [s_b0_j0, s_b0_j1]
        # To gamma_b1 we must add streets s_b1_j0, s_b1_j1 
        elif s_b1 in h_v:
            gamma_b1 = h_k
            new_streets_gamma_b1 = [s for s in h_v] + [s_b1_j0, s_b1_j1]
        # To gamma_j0 we must remove streets s_b0_j0, s_b1_j0
        if s_j0 in h_v:
            gamma_j0 = h_k
            new_streets_gamma_j0 = [
                s for s in h_v if (s != s_b0_j0 and s != s_b1_j0)
            ]
        # To gamma_j1 we must remove streets s_b0_j1, s_b1_j1
        if s_j1 in h_v:
            gamma_j1 = h_k
            new_streets_gamma_j1 = [
                s for s in h_v if (s != s_b0_j1 and s != s_b1_j1)
            ]
    # Now update the new homology classes:
    new_homology_classes[gamma_b0] = new_streets_gamma_b0
    new_homology_classes[gamma_b1] = new_streets_gamma_b1
    new_homology_classes[gamma_j0] = new_streets_gamma_j0
    new_homology_classes[gamma_j1] = new_streets_gamma_j1

    # the face data remains unchanged, because it identifies a face 
    # only based on the (ordered sequence of) streets that circle around it 
    # counter-clockwise. The streets are unchanged by the cootie.
    new_face_data = (
        {f_k : f_v.street_sequence for f_k, f_v in graph.faces.iteritems()}
    )

    return BPSgraph(
        branch_points=new_branch_points, 
        streets=new_streets, 
        joints=new_joints, 
        homology_classes=new_homology_classes,
        faces=new_face_data,
        bp_coordinates=new_bp_coordinates,
        j_coordinates=new_j_coordinates,
        e_coordinates=new_e_coordinates,
        seed=new_seed,
    )


def get_label(obj):
    """
    if the object is a street, it returns its label
    otherwise returns the object
    """
    if obj.__class__.__name__ == 'Street': 
        return obj.label
    else:
        return obj


def can_be_flipped(street):
    """
    Determines whether an edge can be flipped.
    Minimal requirements are that it is bounded by two branch points,
    and that each branch point is trivalent.
    (the second can be relaxed, but for now we enforce this)
    Also require that an edge does not end twice on the same branch point.
    (cannot flip the internal edge of a self-folded triangle)
    """
    ep0, ep1 = street.endpoints

    # check that both endpoints are tri-valent branch points
    if (
        ep0.end_point.type != 'type_3_branch_point' or
        ep1.end_point.type != 'type_3_branch_point'
    ):
        return False
    elif ep0.end_point == ep1.end_point:
        return False
    else:
        return True


def can_cootie(graph, face):
    """
    Determines whether a cootie can be performed on a face.
    Minimal requirements are that it is bounded by two trivalent 
    branch points, and two trivalent joints.
    (the second can be relaxed, but for now we enforce this)
    """
    
    # check that the face is bounded by two branch points 
    # and two joints which alternate
    if (
        face.face_type != ['j', 'b', 'j', 'b'] and
        face.face_type != ['b', 'j', 'b', 'j']
    ):
        return False
    else:
        pass

    # identify the branch points and the joints involved in the cootie move
    if face.face_type == ['j', 'b', 'j', 'b']:
        b_0 = face.node_sequence[1]
        b_1 = face.node_sequence[3]
        j_0 = face.node_sequence[2]
        j_1 = face.node_sequence[0]
    elif face.face_type == ['b', 'j', 'b', 'j']:
        b_0 = face.node_sequence[0]
        b_1 = face.node_sequence[2]
        j_0 = face.node_sequence[1]
        j_1 = face.node_sequence[3]

    # now I check that the joints are both tri-valent, 
    # and that also the branch points are tri-valent,
    # otherwise we cannot perform a cootie move
    if (
        graph.joints[j_0].type != 'type_3_joint' or
        graph.joints[j_1].type != 'type_3_joint' or
        graph.branch_points[b_0].type != 'type_3_branch_point' or
        graph.branch_points[b_1].type != 'type_3_branch_point'
    ):
        return False
    else:
        return True


    # check that both endpoints are tri-valent branch points
    if (
        ep0.end_point.type != 'type_3_branch_point' or
        ep1.end_point.type != 'type_3_branch_point'
    ):
        return False
    else:
        return True


def can_be_H_flipped(street):
    """
    Determines whether an edge can be 'H-swapped'.
    This is the transformation s-channel --> t-channel for an edge
    that it is bounded by two joints, each of which must be trivalent.
    """
    ep0, ep1 = street.endpoints

    # check that both endpoints are tri-valent branch points
    if (
        ep0.end_point.type != 'type_3_joint' or
        ep1.end_point.type != 'type_3_joint'
    ):
        return False
    elif ep0.end_point == ep1.end_point:
        return False
    else:
        return True


def are_same_cycle(c_1, c_2):
    """
    determine whether c_1 and c_2 are related by a cyclic permutation
    """
    l_1 = len(c_1) 
    l_2 = len(c_2) 
    
    if l_1 != l_2:
        return False
    
    for i in range(l_1):
        c_1_prime = [c_1[(j + i) % l_1] for j in range(l_1)]
        if c_1_prime == c_2:
            return True

    return False

def have_same_face_types(graph_1, graph_2):
    """
    Determine whether two graphs have faces of the 
    same types, in terms of their 'type' as determined
    by the ordered sequence of nodes bounding a face
    e.g. ['b', 'b', 'j', 'j']
    """
    face_types_1 = [f.face_type for f in graph_1.faces.values()]
    face_types_2 = [f.face_type for f in graph_2.faces.values()]
    # for each element in the first list, 
    # sort through the elements of the second list.
    # If a match is found, remove that element from the second 
    # list, then keep checking.
    for i_1, ft_1 in enumerate(face_types_1):
        for i_2, ft_2 in enumerate(face_types_2):
            if are_same_cycle(ft_1, ft_2):
                face_types_2.pop(i_2)
                break
    if len(face_types_2) > 0:
        return False
    else:
        return True


def validate_edge_perm(graph_1, graph_2, edge_perm_dic):
    """
    check if the permutation of edges encoded by a dictionary 
    of the form {...,  edge_of_g1 : corresp_edge_of_g2, ...}
    is valid, based on joint and branch point structure
    """
    s_1 = graph_1.streets.keys()
    bp_1 = {
        b.label : [s.label for s in b.streets] 
        for b in graph_1.branch_points.values()
    }
    j_1 = {
        j.label : [get_label(s) for s in j.streets] 
        for j in graph_1.joints.values()
    }

    s_2 = graph_2.streets.keys()
    bp_2 = {
        b.label : [s.label for s in b.streets] 
        for b in graph_2.branch_points.values()
    }
    j_2 = {
        j.label : [get_label(s) for s in j.streets] 
        for j in graph_2.joints.values()
    }

    # Now apply the permutation to the set of streets s_1
    # and recast branch points bp_1 and joints j_1
    # in the new labels
    # If they both coincide with bp_2 and j_2 respectively,
    # the two graphs are equivalent!

    # list of permuted edges, in the order corresponding to s_1:
    p = [edge_perm_dic[s] for s in s_1]

    # ### DEBUG
    # print '\n\npermuting edges \n{}\ninto\n{}\n'.format(s_1, p)
    # ###

    # start by checking equality of branch points
    bp_1_perm = replace(s_1, p, bp_1)
    # ### DEBUG
    # print 'permuting BP \n{}\ninto\n{}\ncompare\n{}'.format(bp_1, bp_1_perm, bp_2)
    # ###
    if equivalent_as_dictionaries(bp_1_perm, bp_2) is True:
        # then also check equality of joints
        j_1_perm = replace(s_1, p, j_1)
        # ### DEBUG
        # print 'permuting J \n{}\ninto\n{}\ncompare\n{}'.format(j_1, j_1_perm, j_2)
        # ###
        if equivalent_as_dictionaries(j_1_perm, j_2) is True:
            return p
        else:
            return None
    else:
        return None



def are_equivalent_by_edge_perm(graph_1, graph_2):
    """
    Determine whether two graphs are equivalent, in the 
    sense that there exists a permutation of the edges
    which makes them identical.
    If they are inequivalent, returns None.
    If they are equivalent, returns all the permutations of edges
    which make those graphs equivalent.
    """

    s_1 = graph_1.streets.keys()

    # create all permutations of the list of streets
    all_perms = map(list, list(itertools.permutations(
        s_1
    )))

    ok_perms = []

    for p in all_perms:
        # construct the permutation dictionary
        edge_perm_dic = {s_1[i] : p[i] for i in range(len(s_1))}
        # check if the permutation is valid
        check = validate_edge_perm(graph_1, graph_2, edge_perm_dic)
        if check is not None:
            ok_perms.append([s_1, p])

    if len(ok_perms) > 0:
        return ok_perms
    else:
        return None


def equivalent_as_dictionaries(dic_1, dic_2):
    """
    Determine whether two dictionaries (of branch points or joints)
    are equivalent.
    Include the possibility of relabeling of branch points/joints
    as well as the possibility of cyclic shifts among the streets 
    ending at a node.
    In practice, we just compare dic.values() and see if the two sets 
    are equivalent (up to cyclic permutations of the entries)
    """
    v_1 = dic_1.values()
    v_2 = dic_2.values()
    # for each element in the first list, 
    # sort through the elements of the second list.
    # If a match is found, remove that element from the second 
    # list, then keep checking.
    for i_1, el_1 in enumerate(v_1):
        for i_2, el_2 in enumerate(v_2):
            if are_same_cycle(el_1, el_2):
                v_2.pop(i_2)
                break
    if len(v_2) > 0:
        return False
    else:
        return True

def replace(old_vars, new_vars, dic):
    """
    Replace the variables appearing in the values of a dictionary
    according to a permutation.
    IMPORTANT: this is a specialized functoin, the dictionary is
    assumed to be structured like this:
    {..., key : [entry1, entry2, ..., entryk], ...}
    that is, the values of the dictionary are lists of entries
    """
    if len(old_vars) != len(new_vars):
        raise Exception('Replacement impossible')

    # then, let's perform the subs, first create a replacement dictionary
    rep = {old_vars[i] : new_vars[i] for i in range(len(old_vars))}
    
    new_dic_str = {}
    for key, val in dic.iteritems():
        # build the new value obtained after substitution
        # for each key of the dictionary
        new_val = []
        for entry in val:
            # in Joints, the entries of val = [entry1, entry2,...]
            # may happen to be 'None', so be careful to handle 
            # this separately
            if entry is not None:
                new_val.append(rep[entry])
            else:
                new_val.append(None)
        new_dic_str[key] = new_val

    return new_dic_str


### OLD -- was finding sequences using both face permutations and
### using edge permutations
###
# def find_invariant_sequences(
#     graph, depth, level=None, ref_graph=None, sequence=None,
#     face_cycle_min_length=None, edge_cycle_min_length=None,
# ):
#     """
#     Find all sequences of flips or cootie moves, up to length 
#     'depth' which take a BPS graph back to itself, up to an automorphism.
#     The criterion for the automorphism is that the two graphs must 
#     coincide up to a permutation of the streets.
#     Since the permutation of streets can take much time, 
#     to speed up the comparison we first make weaker checks.
#     For instance we check that two graphs have the same types of faces.
#     Only retain 
#     """
#     if face_cycle_min_length is None:
#         face_cycle_min_length = 0
#     if edge_cycle_min_length is None:
#         edge_cycle_min_length = 0

#     if level is None:
#         level = 0

#     if level == depth:
#         # print 'max depth reached'
#         return []

#     if sequence is None:
#         sequence = []

#     # print '\nStudying sequence {}'.format(sequence)
    
#     if ref_graph is None:
#         ref_graph = graph

#     self_similar_graphs = []

#     for f in graph.mutable_faces:
#         # avoid repeating the last move in the sequence
#         # that would be a trivial result (it's involutive)
#         if len(sequence) > 0:
#             if f == sequence[-1]:
#                 # print 'cootie on {} is the same as last move'.format(f)
#                 continue

#         g_new = cootie_face(graph, f)
#         new_sequence = sequence + [f]
#         if have_same_face_types(ref_graph, g_new) is True:
#             perms = are_equivalent(ref_graph, g_new, faces=True, edges=False)
#             if perms is not None:
#                 # FIXME: this only works if 'perms' is a honest list of 
#                 # permutations
#                 # it won't be true if both edge-permtuations and 
#                 # face-permutations are turned on.
                
#                 # Now among all the permutations which would 
#                 # make two graphs coincide, we take the minimal one
#                 # only. Because we want to discard two graphs which 
#                 # can be related to each other by 'small permutations'
#                 min_perm = []
#                 n_faces = len(graph.faces.keys())
#                 min_l = n_faces
#                 # find the minimal length of all the maximal 
#                 # cycles within each permutation
#                 for p in perms:
#                     cycles = determine_perm_cycles(p[0], p[1])
#                     if max(cycles) < min_l:
#                         min_l = max(cycles)
#                         min_perm = p
#                 # now let's see if there is a permutation with 
#                 # all cycles below the threshold cycle_min_length
#                 # in that case, we don't keep this graph.
#                 if min_l < face_cycle_min_length:
#                         continue
#                 else:
#                     print (
#                         'This is a good sequence: {}'.format(new_sequence)
#                     )
#                     self_similar_graphs.append([new_sequence, p])
#         else:
#             # print 'Will try going deeper with: {}'.format(new_sequence)
#             deeper_sequences = find_invariant_sequences(
#                 g_new,
#                 depth,
#                 level=level+1,
#                 ref_graph=ref_graph,
#                 sequence=new_sequence,
#                 face_cycle_min_length=face_cycle_min_length,
#                 edge_cycle_min_length=edge_cycle_min_length,
#             )
#             self_similar_graphs += deeper_sequences

#     for e in graph.mutable_edges:
#         # avoid repeating the last move in the sequence
#         # that would be a trivial result (it's involutive)
#         # print 'studying mutation on {}'.format(e)
#         if len(sequence) > 0:
#             if e == sequence[-1]:
#                 # print 'flipping {} is the same as last move'.format(e)
#                 continue

#         g_new = flip_edge(graph, e)
#         new_sequence = sequence + [e]
#         if have_same_face_types(ref_graph, g_new) is True:
#             perms = are_equivalent(ref_graph, g_new, faces=True, edges=False)
#             if perms is not None:
#                 # FIXME: this only works if 'perms' is a honest list of 
#                 # permutations
#                 # it won't be true if both edge-permtuations and 
#                 # face-permutations are turned on.
                
#                 # Now among all the permutations which would 
#                 # make two graphs coincide, we take the minimal one
#                 # only. Because we want to discard two graphs which 
#                 # can be related to each other by 'small permutations'
#                 min_perm = []
#                 n_faces = len(graph.faces.keys())
#                 min_l = n_faces
#                 # find the minimal length of all the maximal 
#                 # cycles within each permutation
#                 for p in perms:
#                     cycles = determine_perm_cycles(p[0], p[1])
#                     if max(cycles) < min_l:
#                         min_l = max(cycles)
#                         min_perm = p
#                 # now let's see if there is a permutation with 
#                 # all cycles below the threshold cycle_min_length
#                 # in that case, we don't keep this graph.
#                 if min_l < face_cycle_min_length:
#                         continue
#                 else:
#                     print (
#                         'This is a good sequence: {}'.format(new_sequence)
#                     )
#                     self_similar_graphs.append([new_sequence, p])
                    
#         else:
#             # print 'Will try going deeper with: {}'.format(new_sequence)
#             deeper_sequences = find_invariant_sequences(
#                 g_new,
#                 depth,
#                 level=level+1,
#                 ref_graph=ref_graph,
#                 sequence=new_sequence,
#                 face_cycle_min_length=face_cycle_min_length,
#                 edge_cycle_min_length=edge_cycle_min_length,
#             )
#             self_similar_graphs += deeper_sequences

#     return [ssg for ssg in self_similar_graphs if len(ssg) > 0]


def move_count(sequence, labels_list, diff=None):
    """
    Given a sequence, like 
    ['p_8', 'p_4', 'f_2', 'p_8', 'f_2', 'p_4']
    counts how many occurrences of labels of a certain type there are.
    For example if labels_list is a list of labels of faces, then it counts 
    the number of cootie moves in the sequence.
    Likewise, supplying the list of edges will count flips.

    If diff is True, then only count each occurrence of an edge/face once.
    """
    if diff is None:
        diff=False

    tot = 0
    if diff is False:
        for f in labels_list:
            tot += sequence.count(f)
    elif diff is True:
        for f in labels_list:
            if f in sequence:
                tot += 1
    
    return tot


# def overall_graph_displacement(graph_1, graph_2):
#     """
#     Computes the center of mass of all edges for each graph,
#     returns the modulus of the difference of those vectors.
#     """
#     n_edges = float(len(graph_1.e_coordinates.keys()))
#     # recall that graph.e_coordinates is a dictionary like
#     # {..., p_i : [[x_0, y_0] , [x_1, y_1]], ...}
#     center_1 = sum(
#         [(numpy.array(v_0) + numpy.array(v_1))/ 2 
#         for [v_0 , v_1] in graph_1.e_coordinates.values()]
#     ) / n_edges
#     center_2 = sum(
#         [(numpy.array(v_0) + numpy.array(v_1))/ 2 
#         for [v_0 , v_1] in graph_2.e_coordinates.values()]
#     ) / n_edges

#     return numpy.linalg.norm(center_2 - center_1)


def find_graph_corner(graph, direction=None):
    """
    Find the corner of a given BPS graph, in the direction specified by 
    'direction' which can be: NW, NE, SE, SW.
    Considers only the edges in the fundamental region, 
    as defined by those whose coordinates are given explicitly 
    in the attribute 'graph.e_coordinates'.
    Also determine which of its two nodes is the bottom-right-most.
    """
    if direction == 'NW':
        v = [-1, 1]
    elif direction == 'NE':
        v = [1, 1]
    elif direction == 'SE':
        v = [1, -1]
    elif direction == 'SW':
        v = [-1, -1]

    # set a starting point for the search
    # recall that graph.e_coordinates[p_i] is of the form 
    # [[x_0, y_0], [x_1, y_1]]
    edge_label = graph.e_coordinates.keys()[0]
    node_coord = graph.e_coordinates[edge_label][0]
    # compute the 'height' of the node, in the bottom-right direction,
    # by projecting onto the vector [1, -1]
    node_height = v[0] * node_coord[0] + v[1] * node_coord[1]

    for e_k, e_v in graph.e_coordinates.iteritems():
        # recall that e_v is of the form [[x_0, y_0], [x_1, y_1]]
        # compute the 'height' of each node, in the bottom-right direction,
        # by projecting onto the vector [1, -1]
        h_0 = v[0] * e_v[0][0] + v[1] * e_v[0][1]
        h_1 = v[0] * e_v[1][0] + v[1] * e_v[1][1]
        if h_0 > h_1:
            ind = 0
            max_h = h_0
        else:
            ind = 1
            max_h = h_1

        # if this node sits 'higher', retain it as the candidate
        if max_h > node_height:
            edge_label = e_k
            node_coord = e_v[ind]
            node_height = max_h

    # now determine the label of the node
    found=False
    for bp_k, bp_v in graph.branch_points.iteritems():
        if are_within_range(graph.bp_coordinates[bp_k], node_coord, EPSILON):
            if found is False:
                found = True
                node_label = bp_k
            else:
                raise Exception(
                    'Two candidate nodes seem to sit at {}'.format(node_coord)
                )
    for j_k, j_v in graph.joints.iteritems():
        if are_within_range(graph.j_coordinates[j_k], node_coord, EPSILON):
            if found is False:
                found = True
                node_label = j_k
            else:
                raise Exception(
                    'Two candidate nodes seem to sit at {}'.format(node_coord)
                )
    if found is False:
        raise Exception(
            'Cannot find a branch point or joint at {}'.format(node_coord)
        )

    return [edge_label, node_label, node_coord]



def sum_up_edges(edge_sequence, graph, int_round=None):
    """
    Given an ordered sequence of edges [p0,p2,...,pk] 
    checks that they are actually concatenated, and computes 
    the overall displacement from beg(p0) to end(pk).

    The option int_round, if set to True, rounds them to the nearest
    integer.
    """
    p0 = edge_sequence[0]
    p1 = edge_sequence[1]
    # recall that graph.e_coordinates is a list of two vectors for each edge
    # [..., [[x_0, y_0], [x_1, y_1]] , ...]
    v_00 = graph.e_coordinates[p0][0]
    v_01 = graph.e_coordinates[p0][1]
    v_10 = graph.e_coordinates[p1][0]
    v_11 = graph.e_coordinates[p1][1]
    if (
        are_within_range(v_00, v_10, EPSILON) or 
        are_within_range(v_00, v_11, EPSILON)
    ):
        start = v_01
        next = v_00
    elif (
        are_within_range(v_01, v_10, EPSILON) or 
        are_within_range(v_01, v_11, EPSILON)
    ):
        start = v_00
        next = v_01
    else:
        raise Exception(
            'The first two edges in the sequence {} appear not to be '
            'concatenated'.format(edge_sequence)
        )

    delta_x = next[0] - start[0]
    delta_y = next[1] - start[1]

    for i, p in enumerate(edge_sequence[1:]):
        v_p0 = graph.e_coordinates[p][0]
        v_p1 = graph.e_coordinates[p][1]
        if are_within_range(v_p0, next, EPSILON):
            start = v_p0
            next = v_p1
        elif are_within_range(v_p1, next, EPSILON):
            start = v_p1
            next = v_p0
        else:
            raise Exception(
                'The two edges {}, {} in the sequence {} appear not to be '
                'concatenated'.format(
                    edge_sequence[i-1], edge_sequence[i], edge_sequence
                )
            )
        delta_x += next[0] - start[0]
        delta_y += next[1] - start[1]
    
    if int_round is True:
        return [int(round(delta_x)), int(round(delta_y))]
    else:
        return [delta_x, delta_y]


def compute_modular_parameter(one, tau, graph):
    """
    Takes two sequences of streets 'one' and 'tau', sums up their 
    displacements, converts them into complex numbers a, b
    and finally returns the ratio a/b.
    """
    a = sum_up_edges(tau, graph)
    b = sum_up_edges(one, graph)
    return (a[0] + 1j * a[1]) / (b[0] + 1j * b[1])


def quiver_node_permutation(graph_1, graph_2, edge_permutation):
    """
    Given an initial graph, graph1, and a final graph, graph2, 
    related by a permutation of the edges, determine the corresponding
    permutation of the homology classes / quiver nodes.

    More precisely:
    - 'edge_permutation' is a dictionary which, applied to the edges 
    of graph_1 gives graph_2 {..., p : p', ...}
    """
    # let's build the permutation dictionary for homology 
    # classes, by scanning through the edges of each graph
    hc_dict = {}
    for e_1 in graph_1.streets.keys():
        gamma_1 = (
            [key for key, value in graph_1.homology_dictionary().iteritems() 
            if e_1 in value][0]
        )
        
        e_2 = edge_permutation[e_1]
        gamma_2 = (
            [key for key, value in graph_2.homology_dictionary().iteritems() 
            if e_2 in value][0]
        )
        # check if we already found this mapping of homology class,
        # and check that it is compatible
        if gamma_1 in hc_dict.keys():
            if hc_dict[gamma_1] == gamma_2:
                pass
            else:
                raise Exception(
                    'Conflict in assignment of homology class dictionary.'
                )
        elif gamma_2 in hc_dict.values():
            for g1, g2 in hc_dict.iteritems():
                if g2 == gamma_2:
                    if g1 == gamma_1:
                        pass
                    else:
                        raise Exception(
                        'Conflict in assignment of homology class dictionary.'
                    )
        else:
            hc_dict[gamma_1] = gamma_2


    return hc_dict


def find_invariant_sequences(
    prev_graph, depth, level=None, ref_graph=None, sequence=None,
    mutation_sequence=None,
    edge_cycle_min_length=None, min_n_cooties=None,
    fundamental_region=None, drop_if_trivial_perm=None,
    avoid_last_n_moves=None, avoid_last_n_mutations=None,
):
    """
    Find all sequences of flips or cootie moves, up to length 
    'depth' which take a BPS graph back to itself, up to an automorphism.
    The criterion for the automorphism is that the two graphs must 
    coincide up to a permutation of the streets.
    We try to match streets on both graphs by following the node structure.
    Only retain sequences in which streets are permuted with cycles 
    of a given minimal length.

    We avoid sequences whose corresponding mutations include two identical
    consecutive elements.

    'level' is an auxiliary variable, it is used to keep track of recursion.
    'ref_graph' is the original graph, a sequence must eventually reproduce it.
    'sequence' is also an auxiliary variable, it keeps track of how many moves
    have been done, in fact level=len(sequence)
    'mutation_sequence' contains the same data as 'sequence', except that
    it translates it into actual quiver mutations, that is: it doesn't
    record cooties, and it translates flips on streets into mutations on 
    quvier nodes.
    'edge_cycle_min_length' determines whether to retain a valid sequence, 
    based on how the edges of the graph are permuted: if it's a cyclic 
    permutation of length greater than a certain minimal value
    'min_n_cooties' determines whether to retain a valid sequence based on the
    number of cootie moves it contains
    'fundamental_region' is a list [one, tau] which contains information about 
    the fundamental region of the torus where the reference BPS graph sits.
    For example, in the [2,1] torus we would have
    one = ['p_1', 'p_8', 'p_9', 'p_2']
    tau = ['p_2', 'p_5', 'p_4', 'p_1']
    'drop_if_trivial_perm' drops a sequence of moves if there exists a trivial 
    permutation which relates its outcome to the original graph.

    """
    # ### DEBUG TEMP
    # temp_perms = are_equivalent_as_graphs(prev_graph, prev_graph)
    # if temp_perms is not None:
    #     if len(temp_perms) > 2:
    #         print '\n\n\nTHIS GRAPHS HAS A Z{} SYMMETRY\n{}\n\n\n'.format(
    #             len(temp_perms), sequence
    #             )
    # #### END DEBUG
    if edge_cycle_min_length is None:
        edge_cycle_min_length = 0

    if min_n_cooties is None:
        min_n_cooties = 0

    if level is None:
        level = 0

    if drop_if_trivial_perm is None:
        drop_if_trivial_perm = False

    if level == depth:
        # print 'max depth reached'
        return []

    if sequence is None:
        sequence = []

    if mutation_sequence is None:
        mutation_sequence = []
    
    if ref_graph is None:
        ref_graph = prev_graph

    # the graph (as an object) may change during successive iterations
    # hence record now the important data
    mutable_faces = prev_graph.mutable_faces
    mutable_edges = prev_graph.mutable_edges

    self_similar_graphs = []

    for f in mutable_faces:
        # avoid repeating the last move in the sequence
        # that would be a trivial result (it's involutive)
        if len(sequence) > 0:
            if f == sequence[-1]:
                # print 'cootie on {} is the same as last move'.format(f)
                continue

        # avoid repeating any move in the last k ones, 
        # where k is given as an argument of this function
        if avoid_last_n_moves is not None:            
            prev_n = min(len(sequence), avoid_last_n_moves)
            if f in sequence[-prev_n:]:
                continue

        new_graph = cootie_face(prev_graph, f)
        ### DEBUG
        # if sequence == ['p_6', 'p_2'] and f == 'f_1':
        #     print '\n\n\nlevel {}: cootie on face {}'.format(level, f)
        #     print '\n\n\tNEW GRAPH\n\n'
        #     new_graph.print_info()
        ### DEBUG END
        new_sequence = [m for m in sequence] + [f]
        # Cooties are not recorded as mutations
        new_mutation_sequence = [m for m in mutation_sequence]

        if have_same_face_types(ref_graph, new_graph) is True:
            perms = are_equivalent_as_graphs(ref_graph, new_graph)
            if perms is not None:
                # get ALL valid permutation dictionaries of edges
                
                # Now among all the permutations which would 
                # make two graphs coincide, we take only those whose cycles
                # have at least a certain length.
                # Because we want to discard two graphs which 
                # can be related to each other by 'small permutations'
                
                # Due to discrete symmetries of the BPS graph
                # it is typically the case that there is more than 
                # one permutation that establishes equivalence between 
                # two graphs.
                # Among them, we keep only those whose cycles 
                # contain at least one of length 'edge_cycle_min_length'
                selected_perms = []
                for p in perms:
                    # find the minimal length of all the maximal 
                    # cycles within each permutation.
                    # Recall each permutation is a dictionary
                    p_k = p.keys()
                    p_v = [p[key] for key in p_k]
                    cycles = determine_perm_cycles(p_k, p_v)
                    if max(cycles) < edge_cycle_min_length:
                            continue
                    else:
                        selected_perms.append(p)
            else:
                selected_perms = []

            if len(selected_perms) > 0:
                print (
                    'This is a good sequence: {}'.format(new_sequence)
                )
                print 'With these permutations:'
                for p in selected_perms:
                    print p

                # ### DEBUG : Disable for now
                homology_permutations = ([
                    quiver_node_permutation(ref_graph, new_graph, p) 
                    for p in selected_perms
                ])
                # ### DEBUG END
                self_similar_graphs.append(
                    MoveSequence(
                        moves=new_sequence, 
                        edge_permutations=selected_perms,
                        quiver_mutations=new_mutation_sequence, 
                        homology_permutations=homology_permutations,
                    )
                )
        
        # If the graphs had the same face types, but no permutation 
        # related them, keep going deeper.
        # Also if they didn't have the same face types, keep going deeper.
        
        # print 'Will try going deeper with: {}'.format(new_sequence)
        deeper_sequences = find_invariant_sequences(
            new_graph,
            depth,
            level=level+1,
            ref_graph=ref_graph,
            sequence=new_sequence,
            mutation_sequence=new_mutation_sequence,
            edge_cycle_min_length=edge_cycle_min_length,
            min_n_cooties=min_n_cooties,
            avoid_last_n_moves=avoid_last_n_moves,
            avoid_last_n_mutations=avoid_last_n_mutations,
        )
        self_similar_graphs += deeper_sequences

    for e in mutable_edges:
        # Note: this includes both flips on I-webs
        # and 'H-flips' on the middle edges of H-webs.

        # avoid repeating the last move in the sequence
        # that would be a trivial result (it's involutive)
        # print 'studying mutation on {}'.format(e)
        if len(sequence) > 0:
            if e == sequence[-1]:
                # print 'flipping {} is the same as last move'.format(e)
                continue

        # avoid repeating any move in the last k ones, 
        # where k is given as an argument of this function
        if avoid_last_n_moves is not None:
            prev_n = min(len(sequence), avoid_last_n_moves)
            if e in sequence[-prev_n:]:
                continue

        new_graph = flip_edge(prev_graph, e)
        ### DEBUG
        # if sequence == ['p_6'] and e == 'p_2':
        #     print '\n\n\nlevel {}: flip on face {}'.format(level, e)
        #     print '\n\n\tNEW GRAPH\n\n'
        #     new_graph.print_info()
        ### DEBUG END
        new_sequence = [m for m in sequence] + [e]

        # if the edge is an I-web (as opposed to the middle of an H-web)
        # record this into the mutation sequence as well.
        p = prev_graph.streets[e]
        ep0, ep1 = p.endpoints
        if (
            ep0.end_point.__class__.__name__ == 'BranchPoint' and 
            ep1.end_point.__class__.__name__ == 'BranchPoint'
        ):
            gamma = (
                [key for key, value 
                in prev_graph.homology_dictionary().iteritems() 
                if e in value][0]
            )
            # if gamma is the same as the last entry of the mutation sequence, 
            # stop exploring these sequence and its descendants.
            if len(mutation_sequence) > 0:
                if gamma == mutation_sequence[-1]:
                    continue
                # additionally check if we are repeating one of the last n 
                # quiver mutations
                if avoid_last_n_mutations is not None:
                    prev_n = min(
                        len(mutation_sequence), avoid_last_n_mutations
                    )
                    if gamma in mutation_sequence[-prev_n:]:
                        continue
                    else:
                        pass
                else:
                    pass

            new_mutation_sequence = [m for m in mutation_sequence] + [gamma]
        # otherwise, if it's a H-web, just keep the mutation sequence as it was
        elif (
            ep0.end_point.__class__.__name__ == 'Joint' and 
            ep1.end_point.__class__.__name__ == 'Joint'
        ):
            new_mutation_sequence = [m for m in mutation_sequence] 

        # ### DEBUGGING BEGINS
        # # print out the problematic sequence and save the graphs in plots
        # if new_graph.wrong_coordinates is True:
        #     print '\n\n\nPROBLEM with sequence {}'.format(new_sequence)
        #     mydir = os.path.join(
        #         os.getcwd(), 'mcg_moves', 
        #         (datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S') + '_DEBUG')
        #     )
        #     os.makedirs(mydir)
        #     print 'will save plots for debug, in folder {}'.format(mydir)
        #     new_graph = apply_sequence(ref_graph, new_sequence, save_plot=mydir)
        # ### DEBUGGING ENDS
        
        if have_same_face_types(ref_graph, new_graph) is True:
            perms = are_equivalent_as_graphs(ref_graph, new_graph)
            selected_perms = []
            if perms is not None:
                # get ALL valid permutation dictionaries of edges
                
                # Now among all the permutations which would 
                # make two graphs coincide, we take only those whose cycles
                # have at least a certain length.
                # Because we want to discard two graphs which 
                # can be related to each other by 'small permutations'
                
                # Due to discrete symmetries of the BPS graph
                # it is typically the case that there is more than 
                # one permutation that establishes equivalence between 
                # two graphs.
                # Among them, we keep only those whose cycles 
                # contain at least one of length 'edge_cycle_min_length'
                
                for p in perms:
                    # find the minimal length of all the maximal 
                    # cycles within each permutation.
                    # Recall each permutation is a dictionary
                    p_k = p.keys()
                    p_v = [p[key] for key in p_k]
                    cycles = determine_perm_cycles(p_k, p_v)
                    if max(cycles) < edge_cycle_min_length:
                            continue
                    else:
                        selected_perms.append(p)
            if len(selected_perms) > 0:
                print (
                    'This is a candidate sequence: {}'.format(new_sequence)
                )
                # ### DEBUG : Disable for now
                homology_permutations = ([
                    quiver_node_permutation(ref_graph, new_graph, p) 
                    for p in selected_perms
                ])
                # ### DEBUG : END
                self_similar_graphs.append(
                    MoveSequence(
                        moves=new_sequence, 
                        edge_permutations=selected_perms,
                        quiver_mutations=new_mutation_sequence, 
                        homology_permutations=homology_permutations,
                    )
                )

        # If the graphs had the same face types, but no permutation 
        # related them, keep going deeper.
        # Also if they didn't have the same face types, keep going deeper.
        
        # print 'Will try going deeper with: {}'.format(new_sequence)
        deeper_sequences = find_invariant_sequences(
            new_graph,
            depth,
            level=level+1,
            ref_graph=ref_graph,
            sequence=new_sequence,
            mutation_sequence=new_mutation_sequence,
            edge_cycle_min_length=edge_cycle_min_length,
            min_n_cooties=min_n_cooties,
            avoid_last_n_moves=avoid_last_n_moves,
            avoid_last_n_mutations=avoid_last_n_mutations,
        )
        self_similar_graphs += deeper_sequences

    # drop empty sequences
    valid_ssg = (
        [ssg for ssg in self_similar_graphs 
        if len(ssg.moves_sequence) > 0]
    )

    # determine which sequences to retain, based on several factors
    # But only at the very end (if the 'level' is 0)
    if level == 0:
        faces = ref_graph.faces.keys()
        edges = ref_graph.streets.keys()
        
        # require a certain number of cootie moves;
        # also drop a sequence if there exists a trivial permutation 
        # which relates it to the original graph
        trivial_perm = {s : s for s in edges}
        selected_ssg = []
        for ssg in valid_ssg:
            if move_count(ssg.moves_sequence, faces) >= min_n_cooties :
                perms = ssg.edge_permutations
                if trivial_perm not in perms or drop_if_trivial_perm is False:
                    selected_ssg.append(ssg)
                else:
                    print (
                        'Dropping sequence because it is related to the '
                        'original graph by a trivial permutation'
                    )
                    pass
            else:
                pass
        valid_ssg = selected_ssg

    # Now study the new modular parameter for each sequence.
    # But only at the very end (if the 'level' is 0)
    # And ONLY if the final graph is equivalent to the original one
    # (NOTE: being a bit sloppy, we just require that they are EXACTLY the 
    # same object, we don't actually check equivalence)
    if level == 0 and ref_graph == prev_graph:
        old_modular_parameter = compute_modular_parameter(one, tau, ref_graph)
        for i, s in enumerate(valid_ssg):
            # Each sequence of moves comes with a set of permutations
            # that relate the final graph to the original one.
            # Must consider each of them
            for j, p in enumerate(s.edge_permutations):
                # use the infor about the funadmental region to compute the 
                # new modular parameter, by tracking where the edges 
                # ended up.
                # Before the moves, the two parameteters [1, \tau] were 
                # specified by a sequence of edges each.
                # Now tracking where those edges ended up we are
                # able to compute the new [1', \tau'].
                # FIrst of all, we translate the old sequence of edges
                # into the new sequence, by applying the permutation 
                # dictionary
                new_one = [p[edge] for edge in one]
                new_tau = [p[edge] for edge in tau]
                new_graph = apply_sequence(ref_graph, s.moves_sequence)
                new_modular_parameter = compute_modular_parameter(
                    new_one, new_tau, new_graph
                )
                print (
                    'Modular parameter of sequence #{} '
                    'with permutation #{} : {}'
                    .format(i, j, new_modular_parameter)
                )
                print (
                    'The new one : {} = {}\nThe new tau : {} = {}'.format(
                        sum_up_edges(new_one, new_graph, int_round=True), new_one,  
                        sum_up_edges(new_tau, new_graph, int_round=True), new_tau
                    )
                )

    return valid_ssg


def apply_sequence(graph, seq, save_plot=None, include_mutation_sequence=None):
    """
    Perform a sequence of moves.
    Gives an option to save the plots of all steps, 
    takes the directory path as argument for this.

    An option to return also the mutation sequence for homology classes 
    is provided.
    """
    if save_plot is not None:
        plot = plot_BPS_graph(graph)
        plt.savefig(os.path.join(save_plot, '0.png'))
        plt.clf()  # Clear the figure for the next loop

    # collect all names of edges and faces
    # although the graph will change at each step,
    # the set of names will remain invariant.
    faces = [k for k in graph.faces.keys()]
    edges = [k for k in graph.streets.keys()]

    graph_sequence = [graph]
    hc_sequence = []
    for i, el in enumerate(seq):
        if el in faces:
            # print '\n\ncootie on {}'.format(el)
            graph_sequence.append(cootie_face(graph_sequence[-1], el))
        elif el in edges:
            # print '\n\nflip on {}'.format(el)
            # compute which homology class corresponds to the flipped edge
            prev_graph = graph_sequence[-1]
            p = prev_graph.streets[el]
            ep0, ep1 = p.endpoints
            if (
                ep0.end_point.__class__.__name__ == 'BranchPoint' and 
                ep1.end_point.__class__.__name__ == 'BranchPoint'
            ):
                gamma = (
                    [key for key, value 
                    in prev_graph.homology_dictionary().iteritems() 
                    if el in value][0]
                )
                hc_sequence.append(gamma)
            #
            graph_sequence.append(flip_edge(graph_sequence[-1], el))
        else:
            print 'got this command: {}'.format(el)
            raise Exception('Unknown step.')
        
        if save_plot is not None:
            plot = plot_BPS_graph(graph_sequence[-1])
            plt.savefig(os.path.join(save_plot, str(i + 1)+'.png'))
            plt.clf()  # Clear the figure for the next loop
        else:
            pass

        pass   

    if include_mutation_sequence is True:
        return [graph_sequence[-1], hc_sequence]
    else:
        return graph_sequence[-1]



def check_sequence(graph, seq):
    """
    Perform a sequence of moves, then check if the final graph is
    equivalent to the initial one.
    """
    new_graph = apply_sequence(graph, seq)
    perm = are_equivalent(graph, new_graph, faces=True, edges=False)
    if perm is not None:
        print 'The sequence produces an equivalent graph.'
        print 'The permutation is: {}'.format(perm)
        return True
    else:
        print 'The sequence produces an inequivalent graph.'
        return False

def faces_across_edge(graph, street):
    """
    Returns the faces that bound an edge.
    If an edge bounds the same face twice it returns that twice.
    """
    faces = []
    for f_k, f_v in graph.faces.iteritems():
        if street in f_v.street_sequence:
            faces.append(f_k)

    # check if the edge is included in two faces, or only in one
    if len(faces) == 2:
        return faces
    elif len(faces) == 1:
        return faces + faces
    else:
        raise Exception(
            'Edge {} doesnt delimit faces correctly'.format(street)
        )

def other_face(graph, street, face):
    """
    Returns the other face across an edge which bounds it.
    If an edge bounds the same face twice it returns the original face.
    """
    faces = faces_across_edge(graph, street)

    if faces.count(face) == 0:
        raise Exception('Edge does not bound face.')        
    elif faces.count(face) == 1:
        if faces[0] == face:
            return faces[1]
        else:
            return faces[0]
    elif faces.count(face) == 2:
        return face

def face_neighbors(graph, face):
    """
    Return the neighboring faces around a given face,
    ordered counter-clockwise.
    """ 
    boundary_edges = graph.faces[face].street_sequence
    neighbors = [other_face(graph, s, face) for s in boundary_edges]

    return neighbors


def are_equivalent_by_face_perm(graph_1, graph_2):
    """
    Determine whether two graphs are equivalent, in the 
    sense that there exists a permutation of the faces
    which makes them identical.
    If they are inequivalent, returns None.
    If they are equivalent, returns all the permutations of faces
    which make those graphs equivalent.
    """
    f_1 = graph_1.faces.keys()
    t_1 = graph_1.tessellation_data

    f_2 = graph_2.faces.keys()
    t_2 = graph_2.tessellation_data

    # Now apply a permutation to the set of faces f_1
    # and recast the tessellation data t_1
    # in the new labels.
    # If it coincides with t_2 the two graphs are equivalent!

    # create all permutations of the list of streets
    all_perms = map(list, list(itertools.permutations(
        f_1
    )))

    ok_perms = []

    for p in all_perms:
        t_1_perm = replace(f_1, p, t_1)
        if equivalent_as_dictionaries(t_1_perm, t_2) is True:
            ok_perms.append([f_1, p])

    if len(ok_perms) > 0:
        return ok_perms
    else:
        return None

def grow_dict(
    edge_dic, next_node_1, next_node_2
):
    """
    takes a dictionary of edges between two graphs g1, g2
    and keeps growing it from a specified pair of end_points
    on each graph. (IMPORTANT: these are objects of the type
    StreetEndPoint)
        *** IMPORTANT: ASSUMING TRI-VALENT BP AND JOINTS ***
    """
    new_dic = edge_dic.copy()

    # distinguish between branch points and joints
    if next_node_1.end_point.__class__.__name__ == 'BranchPoint':
        q_1 = 1
    elif next_node_1.end_point.__class__.__name__ == 'Joint':
        q_1 = 2
    if next_node_2.end_point.__class__.__name__ == 'BranchPoint':
        q_2 = 1
    elif next_node_2.end_point.__class__.__name__ == 'Joint':
        q_2 = 2

    if q_1 != q_2:
        # interrupt this cycle because the two streets cannot be 
        # compatible with each other if their endpoints are different
        # print 'incompatibility 1'
        return None

    # now that we identified the corresponding endpoints on each
    # graph, let's use this data to identify the new edges
    slot_1 = next_node_1.slot
    next_edge_1_R = next_node_1.end_point.streets[
        (slot_1 + q_1) % (3 * q_1)
    ]
    next_edge_1_L = next_node_1.end_point.streets[
        (slot_1 + 2 * q_1) % (3 * q_1)
    ]
    slot_2 = next_node_2.slot
    next_edge_2_R = next_node_2.end_point.streets[
        (slot_2 + q_2) % (3 * q_2)
    ]
    next_edge_2_L = next_node_2.end_point.streets[
        (slot_2 + 2 * q_2) % (3 * q_2)
    ]
    # add new entries to the dictionary
    if have_compatible_endpoints(
        next_edge_1_L, next_edge_2_L
    ):
        if (
            (next_edge_1_L.label not in new_dic.keys()) and
            (next_edge_2_L.label not in new_dic.values())
        ): 
            # the pair is not yet in the dictionary, add it
            new_dic[next_edge_1_L.label] = next_edge_2_L.label
        elif (
            (next_edge_1_L.label in new_dic.keys()) and
            (next_edge_2_L.label in new_dic.values())
        ): 
            # the pair is already in the dictionary
            if new_dic[next_edge_1_L.label] == next_edge_2_L.label:
                pass
            else:
                # pair is already there but there is a mismatch
                # print 'incompatibility 2'
                return None
        else:
            # only part of the pair is there, so there is an incompatibility
            # print 'incompatibility 3'
            return None
    else:
        # interrupt this cycle because the two streets cannot be 
        # compatible with each other if their endpoints are different
        # print 'incompatibility 4'
        return None

    if have_compatible_endpoints(
        next_edge_1_R, next_edge_2_R
    ):       
        if (
            (next_edge_1_R.label not in new_dic.keys()) and
            (next_edge_2_R.label not in new_dic.values())
        ): 
            # the pair is not yet in the dictionary, add it
            new_dic[next_edge_1_R.label] = next_edge_2_R.label
        elif (
            (next_edge_1_R.label in new_dic.keys()) and
            (next_edge_2_R.label in new_dic.values())
        ): 
            # the pair is already in the dictionary
            if new_dic[next_edge_1_R.label] == next_edge_2_R.label:
                pass
            else:
                # pair is already there but there is a mismatch
                # print 'incompatibility 5'
                return None
        else:
            # only part of the pair is there, so there is an incompatibility
            # print 'incompatibility 6'
            return None
    else:
        # interrupt this cycle because the two streets cannot be 
        # compatible with each other if their endpoints are different
        # print 'incompatibility 7'
        return None

    # move on to the next step, next generation that is.
    # But only if we changed the dictionary at this step, 
    # otherwise we are done.
    if new_dic == edge_dic:
        return new_dic
    else:
        pass

    next_end_point_1_L = other_node(next_edge_1_L, next_node_1.end_point)
    next_end_point_2_L = other_node(next_edge_2_L, next_node_2.end_point)

    next_end_point_1_R = other_node(next_edge_1_R, next_node_1.end_point)
    next_end_point_2_R = other_node(next_edge_2_R, next_node_2.end_point)
    
    new_entries_L = grow_dict(
        new_dic, next_end_point_1_L, next_end_point_2_L
    )
    if new_entries_L == None:
        # there is an incompatibility, stop here
        # print 'incompatibility 8'
        return None
    else:
        # merge with previous dictionary
        updated_new_dic = dict(new_entries_L.items() + new_dic.items())
        new_dic = updated_new_dic.copy()

    new_entries_R = grow_dict(
        new_dic, next_end_point_1_R, next_end_point_2_R
    )
    if new_entries_R == None:
        # there is an incompatibility, stop here
        # print 'incompatibility 9'
        return None
    else:
        # merge with previous dictionary
        updated_new_dic = dict(new_entries_R.items() + new_dic.items())
        new_dic = updated_new_dic

    return new_dic


def other_node(street, node):
    """
    given an actual branch point or joint, returns the opposite StreetEndPoint
    """
    pts = street.endpoints
    if node == pts[0].end_point:
        return pts[1]
    elif node == pts[1].end_point:
        return pts[0]
    else:
        raise Exception()

def node_coordinates(node_label, graph):
    """
    Given a label of a branch point, or a joint, returns the corresponding 
    coordinates as retrieved from the data of the BPSgraph instance 
    given as 'graph'
    """
    if (
        node_label in graph.branch_points.keys() and 
        node_label not in graph.joints.keys()
    ):
        return graph.bp_coordinates[node_label]
    elif (
        node_label not in graph.branch_points.keys() and 
        node_label in graph.joints.keys()
    ):
        return graph.j_coordinates[node_label]
    elif (
        node_label in graph.branch_points.keys() and 
        node_label in graph.joints.keys()
    ):
        raise Exception(
            'Node {} appears to be both branch point and joint'.format(
                node_label
            ))
    else:
        raise Exception(
            'Node {} does not appear among branch points nor joints'.format(
                node_label
            ))


def match_graphs_from_starting_point(
    g1, g2, start_2, orientation_2
):
    """
    Try to match graphs g1, g2 by taking the first edge of graph g1 
    with canonical orientation, and match it with street 'start_2' of graph 
    g2, with orientation specified by 'orientation_2'. 
    Then follow through all the nodes of the graph and see if it works out.
    """
    # ### DEBUG BEGIN
    # print '\n\nwill try to match this graph {}'.format({j : [s.label for s in g1.joints[j].streets if s is not None] for j in g1.joints.keys()})
    # print 'to this one {}'.format({j : [s.label for s in g2.joints[j].streets if s is not None] for j in g1.joints.keys()})
    # ### DEBUG END
    edges_1 = g1.streets.keys()
    edges_2 = g2.streets.keys()

    # fix the initial edge
    if have_compatible_endpoints(
        g1.streets[edges_1[0]], g2.streets[start_2]
    ):
        edge_dic = {edges_1[0] : start_2}
    else:
        # interrupt this because the two streets cannot be 
        # compatible with each other if their endpoints are different
        return None

    # then use the orientation to move forward along the graph
    if orientation_2 == +1:
        next_node_1 = g1.streets[edges_1[0]].initial_point()
    elif orientation_2 == -1:
        next_node_1 = g1.streets[edges_1[0]].final_point()
    next_node_2 = g2.streets[start_2].final_point()
    
    # distinguish between branch points and joints
    if next_node_1.end_point.__class__.__name__ == 'BranchPoint':
        q_1 = 1
    elif next_node_1.end_point.__class__.__name__ == 'Joint':
        q_1 = 2
    if next_node_2.end_point.__class__.__name__ == 'BranchPoint':
        q_2 = 1
    elif next_node_2.end_point.__class__.__name__ == 'Joint':
        q_2 = 2

    if q_1 != q_2:
        # interrupt this because the two streets cannot be 
        # compatible with each other if their endpoints are different
        return None

    # now that we identified the corresponding endpoints on each
    # graph, let's use this data to identify the new edges
    slot_1 = next_node_1.slot
    next_edge_1_R = next_node_1.end_point.streets[
        (slot_1 + q_1) % (3 * q_1)
    ]
    next_edge_1_L = next_node_1.end_point.streets[
        (slot_1 + 2 * q_1) % (3 * q_1)
    ]
    slot_2 = next_node_2.slot
    next_edge_2_R = next_node_2.end_point.streets[
        (slot_2 + q_2) % (3 * q_2)
    ]
    next_edge_2_L = next_node_2.end_point.streets[
        (slot_2 + 2 * q_2) % (3 * q_2)
    ]
    # add both new entries to the dictionary
    if (
        have_compatible_endpoints(next_edge_1_L, next_edge_2_L) and
        have_compatible_endpoints(next_edge_1_R, next_edge_2_R)        
    ):
        edge_dic[next_edge_1_L.label] = next_edge_2_L.label
        edge_dic[next_edge_1_R.label] = next_edge_2_R.label
    else:
        # interrupt this because the two streets cannot be 
        # compatible with each other if their endpoints are different
        return None

    # move on to the next step, next generation that is.
    next_end_point_1_L = other_node(next_edge_1_L, next_node_1.end_point)
    next_end_point_2_L = other_node(next_edge_2_L, next_node_2.end_point)

    next_end_point_1_R = other_node(next_edge_1_R, next_node_1.end_point)
    next_end_point_2_R = other_node(next_edge_2_R, next_node_2.end_point)

    # print '1 the dictionary is : {}'.format(edge_dic)

    new_entries_L = grow_dict(
        edge_dic, next_end_point_1_L, next_end_point_2_L
    )
    if new_entries_L == None:
        # there is an incompatibility, stop here
        return None
    else:
        # merge with previous dictionary
        new_dic = dict(new_entries_L.items() + edge_dic.items())
        edge_dic = new_dic
        # print '2 the dictionary is : {}'.format(edge_dic)

    new_entries_R = grow_dict(
        edge_dic, next_end_point_1_R, next_end_point_2_R
    )
    if new_entries_R == None:
        # there is an incompatibility, stop here
        return None
    else:
        # merge with previous dictionary
        new_dic = dict(new_entries_R.items() + edge_dic.items())
        edge_dic = new_dic
        # print '3 the dictionary is : {}'.format(edge_dic)

    return edge_dic

def match_graphs_by_edges(g1, g2, all_perms=None):
    """
    Does what the name says.
    Try to match graph g1 to g2 by taking the first edge of g1 
    with fixed orientation and trying to match onto every edge of g2 
    with both orientations, then follow through nodes of each 
    graph in parallel, stopping if a contradiction is found.

    An option to return ALL permutations which work is available.
    """
    if all_perms is None:
        all_perms=False
    elif all_perms is True:
        ok_perms = []

    for start_2 in g2.streets.keys():
        for orientation_2 in [+1, -1]:
            # ### DEBUG BEGIN
            # print 'trying {} with orientation {}'.format(start_2, orientation_2)
            # ### DEBUG END
            edge_dic = match_graphs_from_starting_point(
                g1, g2, start_2, orientation_2
            )
            if edge_dic is not None:
                # ### DEBUG
                # print 'found permutations!'
                # print edge_dic
                # ###
                if all_perms is False:
                    return edge_dic
                elif all_perms is True:
                    ok_perms.append(edge_dic)

    return ok_perms


def are_equivalent(graph_1, graph_2, faces=False, edges=False):
    """
    *** DEPRECATED ***
    Determine if two graphs are equivalent by checking permutations 
    of edges and/or faces.
    The former can be done in two ways:
    a- trying all edge permutations (usually much slower, should avoid)
    b- trying to match streets according to the graph structure (following
        nodes of the graph)
    method 'a' is now deprecated, but still available. 
    method 'b' is always checked.
    """
    f_check = True
    e_check = True

    edge_dictionary = match_graphs_by_edges(graph_1, graph_2)
    if edge_dictionary is not None:
        # (double) check if the dictionary of edges provides the correct
        # branhc point and joint structure
        check = validate_edge_perm(graph_1, graph_2, edge_dictionary)
        if check is not None:
            print '\nOK found a dictionary for edges'
            pass
        else:
            return None
    else:
        return None

    if faces is False and edges is False:
        raise Exception('No checks made.')

    elif faces is True and edges is False:
        f_check = are_equivalent_by_face_perm(graph_1, graph_2)
        return f_check

    elif faces is False and edges is True:
        e_check = are_equivalent_by_edge_perm(graph_1, graph_2)
        return e_check

    elif faces is True and edges is True:
        f_check = are_equivalent_by_face_perm(graph_1, graph_2)
        e_check = are_equivalent_by_edge_perm(graph_1, graph_2)
        if f_check is None or e_check is None:
            return None
        else:
            return f_check + e_check


def are_equivalent_as_graphs(graph_1, graph_2):
    """
    Determine if two graphs are equivalent by matches edges of one 
    onto the other, following the graph structure on each side, node-by-node.
    Returns ALL permutations which make two graphs identical.
    """
    edge_dictionary = match_graphs_by_edges(graph_1, graph_2, all_perms=True)
    if edge_dictionary is not None:
        ok_perms = []
        # (double) check if the dictionary of edges provides the correct
        # branch point and joint structure.
        # Select the valid ones
        for perm in edge_dictionary:
            check = validate_edge_perm(graph_1, graph_2, perm)
            if check is not None:
                # print 'OK found a dictionary for edges'
                ok_perms.append(perm)
                pass
            else:
                # print 'Discard permutation.'
                pass
        if len(ok_perms) > 0:
            print (
                'found {} inequivalent permutations, presumably the graph '
                'has an underlying Z_{} symmetry'.format(
                    len(ok_perms), len(ok_perms)
                )
            )
            return ok_perms
        else:
            return None
    else:
        return None

def determine_perm_cycles(before, after):
    """
    Determine the cycles which occur in a permutation.
    Assume that 'before' and 'after' are two lists containing
    the same elements in different orders.
    """
    l_max = len(before)
    
    # create a replacement dictionary
    rep = {before[i] : after[i] for i in range(len(before))}

    cycles = []
    for e in before:
        counter = 0
        x = e
        for i in range(l_max):
            x = rep[x]
            counter += 1
            if x == e:
                cycles.append(counter)
                break

    return cycles

def have_compatible_endpoints(street_1, street_2):
    """
    Check what the name says.
    """
    ends_1 = [pt.end_point.__class__.__name__ for pt in street_1.endpoints]
    ends_2 = [pt.end_point.__class__.__name__ for pt in street_2.endpoints]

    if ends_1 == ends_2 or ends_1[::-1] == ends_2:
        return True
    else:
        return False

def plot_BPS_graph(graph):
    """
    """
    shifts = [[x, y] for x in [-1, 0, 1] for y in [-1, 0, 1]]

    # plot the branch points
    bp_x_s = [v[0] for v in graph.bp_coordinates.values()]
    bp_y_s = [v[1] for v in graph.bp_coordinates.values()]
    for s in shifts:
        plt.plot(
            [v + s[0] for v in bp_x_s],
            [v + s[1] for v in bp_y_s],
            'rx', mew=3, ms=5,
        )

    # plot the joints
    j_x_s = [v[0] for v in graph.j_coordinates.values()]
    j_y_s = [v[1] for v in graph.j_coordinates.values()]
    for s in shifts:
        plt.plot(
            [v + s[0] for v in j_x_s],
            [v + s[1] for v in j_y_s],
            'bo'
        )

    # plot edges
    for s in graph.streets.values():
        [[x_0, y_0], [x_1, y_1]] = graph.e_coordinates[s.label]

        for d in shifts:
            if d == [0, 0]:
                marker = 'r-'
            else:
                marker = 'b-'
            plt.plot(
                [x_0 + d[0], x_1 + d[0]],
                [y_0 + d[1], y_1 + d[1]],
                marker
            )
            plt.text((x_0 + x_1) /2 + d[0], (y_0 + y_1) / 2 + d[1], s.label)

    # # plot face labels
    # # DEBUG: temporarily disabled, need to take into account shifts too.
    # #
    # for f_k, f_v in graph.faces.iteritems():
    #     f_nodes = [node_coordinates(node, graph) for node in f_v.node_sequence]
    #     x_avg = sum([x[0] for x in f_nodes]) / len(f_nodes)
    #     y_avg = sum([x[1] for x in f_nodes]) / len(f_nodes)

    #     for d in shifts:
    #         plt.text(
    #             x_avg + d[0], y_avg + d[1], f_k
    #         )

    plt.axis([-0.1,1.1,-0.1,1.1])


def closest_on_covering_space(xy_0, xy_1, in_out=None):
    """
    the option 'out' forces the closest match 
    to be chosen outside the fundamental region
    """
    [x_0, y_0] = xy_0
    [x_1, y_1] = xy_1

    if in_out is 'out':
        deltas = [-1, 1]
    elif in_out is 'in':
        deltas = [0]
    elif in_out is None:
        deltas = [-1, 0, 1]

    d_x = abs(x_1 - x_0)
    s_x = 0
    for shift in deltas:
        if abs(x_1 + shift - x_0) <= d_x:
            s_x = shift
            d_x = abs(x_1 + shift - x_0)

    d_y = abs(y_1 - y_0)
    s_y = 0
    for shift in deltas:
        if abs(y_1 + shift - y_0) <= d_y:
            s_y = shift
            d_y = abs(y_1 + shift - y_0)

    return [x_1 + s_x, y_1 + s_y]


def find_self_permutations(graph, fundamental_region):
    """
    Find nontrivial permutations that leave a graph invariant.
    'fundamental_region' is given as a pair of lists of edges representing 
    [one, tau]
    """
    one, tau = fundamental_region
    self_perms = are_equivalent_as_graphs(graph, graph)
    trivial_perm = {s : s for s in graph.streets.keys()}
    nontrivial_self_perms = []
    for p in self_perms:
        if p == trivial_perm:
            pass
        else:
            new_one = [p[edge] for edge in one]
            new_tau = [p[edge] for edge in tau]
            new_modular_parameter = compute_modular_parameter(
                new_one, new_tau, graph
            )
            nontrivial_self_perms.append([p, new_modular_parameter])
    
    return nontrivial_self_perms


def find_sequence_completion(
    graph, partial_sequence, n_moves, one, tau, save_files=None,
    avoid_last_n_moves=None, avoid_last_n_mutations=None, 
    drop_if_trivial_perm=None,
    ):
    """
    Given a BPS graph and a sequence of moves, 
    find all possible completions which bring 
    the BPS graph back to itself, starting from 
    the given partial sequence.
    """
    if drop_if_trivial_perm is None:
        drop_if_trivial_perm = False

    w_partial, mutation_seq = apply_sequence(
        graph, partial_sequence, 
        save_plot=None, 
        include_mutation_sequence=True
    )

    max_n_moves = n_moves
    seq = find_invariant_sequences(
        w_partial, max_n_moves, level=0, ref_graph=graph, 
        edge_cycle_min_length=0,
        min_n_cooties=0,
        fundamental_region=[one, tau],
        drop_if_trivial_perm=drop_if_trivial_perm,
        avoid_last_n_moves=avoid_last_n_moves,
        avoid_last_n_mutations=avoid_last_n_mutations,
    )

    modular_parameter = compute_modular_parameter(one, tau, graph)
    possible_modular_parameters = []

    time = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    if len(seq) > 0:
        if save_files is True or type(save_files) is str:
            if save_files is True:
                mydir = os.path.join(
                    os.getcwd(), 'mcg_moves', 
                    (
                        time + '_SEQ_AUTOCOMP'
                    )
                )
            elif type(save_files) is str:
                mydir = os.path.join(
                    os.getcwd(), 'mcg_moves', 
                    (
                        time + '_' + save_files
                    )
                )

            os.makedirs(mydir)

            # save info about sequences of moves in a text file
            text_file = open(mydir + '/sequence_data.txt', 'w')
            text_file.write('\t\tSequence data\n\n')
            text_file.write(
                'Modular parameter of the original torus: {}'
                '\none : {} = {}\ntau : {} = {}\n\n'.format(
                    modular_parameter, 
                    one, sum_up_edges(one, graph), 
                    tau, sum_up_edges(tau, graph)
                )
            )
            text_file.write('Quiver:\n{}\n\n'.format(graph.seed.quiver))
            
            text_file.write('Found {} sequences which complete {}.\n\n'.format(
                len(seq), partial_sequence
            ))
            text_file.write('Will now study their properties.\n\n')
        
        else:
            mydir=None

        
        for i_s, s in enumerate(seq):
            if mydir is not None:
                mydir_i = os.path.join(mydir, 'sequence_{}'.format(i_s))
                os.makedirs(mydir_i)
            else:
                mydir_i = None

            tot_seq = partial_sequence + s.moves_sequence
            w_new, tot_mutation_seq = apply_sequence(
                graph, tot_seq, 
                save_plot=mydir_i, 
                include_mutation_sequence=True
            )
            text_file.write(
                '\n\n----------------------------------------------\n\n'
            )
            
            text_file.write(
                'Sequence # {} is:\n{}\n\nIt corresponds to the mutations\n{}\n'
                .format(i_s, tot_seq, tot_mutation_seq)
            )

            perms = are_equivalent_as_graphs(graph, w_new)
            if perms is not None:
                text_file.write(
                    '\nThe original graph is matched by applying '
                    'the following permutations\n\n'
                )

                for j, p in enumerate(perms):
                    text_file.write('\nPermutation #{}:\n{}\n\n'.format(j, p))

                    new_one = [p[edge] for edge in one]
                    new_tau = [p[edge] for edge in tau]
                    new_modular_parameter = compute_modular_parameter(
                        new_one, new_tau, w_new
                    )
                    if (
                        c_num_is_in_list(
                            new_modular_parameter, 
                            possible_modular_parameters,
                            0.0001
                        ) is False
                    ):
                        possible_modular_parameters.append(
                            new_modular_parameter
                        )
                    text_file.write(
                        '\nNew modular parameter: {}\n'
                        .format(new_modular_parameter)
                    )
                    text_file.write(
                        'The new one : {} = {}\nThe new tau : {} = {}\n'.format(
                            sum_up_edges(new_one, w_new), new_one,
                            sum_up_edges(new_tau, w_new), new_tau
                        )
                    )
                    text_file.write(
                        'The permutation of homology classes \n{}\n'.format(
                            quiver_node_permutation(graph, w_new, p) 
                        )
                    )
            else:
                text_file.write('found no permutations')

        text_file.write(
            '\n\nPossible new modular parameters:\n{}'.format(
                possible_modular_parameters
            )
        )

    return seq


def c_num_is_in_list(c_num, l, eps):
    """
    determine whether a complex number is in a list or not, 
    given a tolerance 'eps'
    """
    ans = False
    for x in l:
        if abs(x - c_num) < eps:
            ans = True

    return ans


def symmetry_degree(graph):
    """
    return n for a graph with Z_n symmetry
    """
    edge_dictionary = match_graphs_by_edges(graph, graph, all_perms=True)
    if edge_dictionary is not None:
        ok_perms = []
        # (double) check if the dictionary of edges provides the correct
        # branch point and joint structure.
        # Select the valid ones
        for perm in edge_dictionary:
            check = validate_edge_perm(graph, graph, perm)
            if check is not None:
                # print 'OK found a dictionary for edges'
                ok_perms.append(perm)
                pass
            else:
                # print 'Discard permutation.'
                pass
    return len(ok_perms)


def find_graphs_with_symmetry(
    prev_graph, depth, level=None, ref_graph=None, sequence=None,
    mutation_sequence=None,
    edge_cycle_min_length=None, min_n_cooties=None,
    fundamental_region=None, drop_if_trivial_perm=None,
    avoid_last_n_moves=None, avoid_last_n_mutations=None,
):
    """
    Find a sequence of flips or cootie moves, up to length 
    'depth' which takes a BPS graph into a (usually new) BPS back 
    with a higher symmetry (Z3).
    The criterion for the automorphism is that the two graphs must 
    coincide up to a permutation of the streets.
    We try to match streets on both graphs by following the node structure.
    Only retain sequences in which streets are permuted with cycles 
    of a given minimal length.

    We avoid sequences whose corresponding mutations include two identical
    consecutive elements.

    'level' is an auxiliary variable, it is used to keep track of recursion.
    'ref_graph' is the original graph, a sequence must eventually reproduce it.
    'sequence' is also an auxiliary variable, it keeps track of how many moves
    have been done, in fact level=len(sequence)
    'mutation_sequence' contains the same data as 'sequence', except that
    it translates it into actual quiver mutations, that is: it doesn't
    record cooties, and it translates flips on streets into mutations on 
    quvier nodes.
    'edge_cycle_min_length' determines whether to retain a valid sequence, 
    based on how the edges of the graph are permuted: if it's a cyclic 
    permutation of length greater than a certain minimal value
    'min_n_cooties' determines whether to retain a valid sequence based on the
    number of cootie moves it contains
    'fundamental_region' is a list [one, tau] which contains information about 
    the fundamental region of the torus where the reference BPS graph sits.
    For example, in the [2,1] torus we would have
    one = ['p_1', 'p_8', 'p_9', 'p_2']
    tau = ['p_2', 'p_5', 'p_4', 'p_1']
    'drop_if_trivial_perm' drops a sequence of moves if there exists a trivial 
    permutation which relates its outcome to the original graph.

    """
    # ### DEBUG TEMP
    # temp_perms = are_equivalent_as_graphs(prev_graph, prev_graph)
    # if temp_perms is not None:
    #     if len(temp_perms) > 2:
    #         print '\n\n\nTHIS GRAPHS HAS A Z{} SYMMETRY\n{}\n\n\n'.format(
    #             len(temp_perms), sequence
    #             )
    # #### END DEBUG
    if edge_cycle_min_length is None:
        edge_cycle_min_length = 0

    if min_n_cooties is None:
        min_n_cooties = 0

    if level is None:
        level = 0

    if drop_if_trivial_perm is None:
        drop_if_trivial_perm = False

    if level == depth:
        # print 'max depth reached'
        return []

    if sequence is None:
        sequence = []

    if mutation_sequence is None:
        mutation_sequence = []
    
    if ref_graph is None:
        ref_graph = prev_graph

    # the graph (as an object) may change during successive iterations
    # hence record now the important data
    mutable_faces = prev_graph.mutable_faces
    mutable_edges = prev_graph.mutable_edges

    good_sequences = []

    for f in mutable_faces:
        # avoid repeating the last move in the sequence
        # that would be a trivial result (it's involutive)
        if len(sequence) > 0:
            if f == sequence[-1]:
                # print 'cootie on {} is the same as last move'.format(f)
                continue

        # avoid repeating any move in the last k ones, 
        # where k is given as an argument of this function
        if avoid_last_n_moves is not None:            
            prev_n = min(len(sequence), avoid_last_n_moves)
            if f in sequence[-prev_n:]:
                continue

        new_graph = cootie_face(prev_graph, f)
        ### DEBUG
        # if sequence == ['p_6', 'p_2'] and f == 'f_1':
        #     print '\n\n\nlevel {}: cootie on face {}'.format(level, f)
        #     print '\n\n\tNEW GRAPH\n\n'
        #     new_graph.print_info()
        ### DEBUG END
        new_sequence = [m for m in sequence] + [f]
        # Cooties are not recorded as mutations
        new_mutation_sequence = [m for m in mutation_sequence]

        # Check if there is a higher (at least Z3) symmetry
        sym_degree = symmetry_degree(new_graph)

        if sym_degree > 2:
            print (
                'Sequence {} produces a BPS graph with Z_{} symmetry'
                .format(
                    new_sequence, sym_degree
                )
            )
            good_sequences.append(new_sequence)
                
        
        # If there was no higher symmetry, keep going deeper.
        
        # print 'Will try going deeper with: {}'.format(new_sequence)
        deeper_sequences = find_graphs_with_symmetry(
            new_graph,
            depth,
            level=level+1,
            ref_graph=ref_graph,
            sequence=new_sequence,
            mutation_sequence=new_mutation_sequence,
            edge_cycle_min_length=edge_cycle_min_length,
            min_n_cooties=min_n_cooties,
            avoid_last_n_moves=avoid_last_n_moves,
            avoid_last_n_mutations=avoid_last_n_mutations,
        )
        good_sequences += deeper_sequences

    for e in mutable_edges:
        # Note: this includes both flips on I-webs
        # and 'H-flips' on the middle edges of H-webs.

        # avoid repeating the last move in the sequence
        # that would be a trivial result (it's involutive)
        # print 'studying mutation on {}'.format(e)
        if len(sequence) > 0:
            if e == sequence[-1]:
                # print 'flipping {} is the same as last move'.format(e)
                continue

        # avoid repeating any move in the last k ones, 
        # where k is given as an argument of this function
        if avoid_last_n_moves is not None:
            prev_n = min(len(sequence), avoid_last_n_moves)
            if e in sequence[-prev_n:]:
                continue

        new_graph = flip_edge(prev_graph, e)
        ### DEBUG
        # if sequence == ['p_6'] and e == 'p_2':
        #     print '\n\n\nlevel {}: flip on face {}'.format(level, e)
        #     print '\n\n\tNEW GRAPH\n\n'
        #     new_graph.print_info()
        ### DEBUG END
        new_sequence = [m for m in sequence] + [e]

        # if the edge is an I-web (as opposed to the middle of an H-web)
        # record this into the mutation sequence as well.
        p = prev_graph.streets[e]
        ep0, ep1 = p.endpoints
        if (
            ep0.end_point.__class__.__name__ == 'BranchPoint' and 
            ep1.end_point.__class__.__name__ == 'BranchPoint'
        ):
            gamma = (
                [key for key, value 
                in prev_graph.homology_dictionary().iteritems() 
                if e in value][0]
            )
            # if gamma is the same as the last entry of the mutation sequence, 
            # stop exploring these sequence and its descendants.
            if len(mutation_sequence) > 0:
                if gamma == mutation_sequence[-1]:
                    continue
                # additionally check if we are repeating one of the last n 
                # quiver mutations
                if avoid_last_n_mutations is not None:
                    prev_n = min(
                        len(mutation_sequence), avoid_last_n_mutations
                    )
                    if gamma in mutation_sequence[-prev_n:]:
                        continue
                    else:
                        pass
                else:
                    pass

            new_mutation_sequence = [m for m in mutation_sequence] + [gamma]
        # otherwise, if it's a H-web, just keep the mutation sequence as it was
        elif (
            ep0.end_point.__class__.__name__ == 'Joint' and 
            ep1.end_point.__class__.__name__ == 'Joint'
        ):
            new_mutation_sequence = [m for m in mutation_sequence] 

        # Check if there is a higher (at least Z3) symmetry
        sym_degree = symmetry_degree(new_graph)

        if sym_degree > 2:
            print (
                'Sequence {} produces a BPS graph with Z_{} symmetry'
                .format(
                    new_sequence, sym_degree
                )
            )
            good_sequences.append(new_sequence)
                
        
        # If there was no higher symmetry, keep going deeper.
        
        # print 'Will try going deeper with: {}'.format(new_sequence)
        deeper_sequences = find_graphs_with_symmetry(
            new_graph,
            depth,
            level=level+1,
            ref_graph=ref_graph,
            sequence=new_sequence,
            mutation_sequence=new_mutation_sequence,
            edge_cycle_min_length=edge_cycle_min_length,
            min_n_cooties=min_n_cooties,
            avoid_last_n_moves=avoid_last_n_moves,
            avoid_last_n_mutations=avoid_last_n_mutations,
        )
        good_sequences += deeper_sequences

    

    return good_sequences






# # --------- torus graph with two joints and two branch points ----------

# streets = ['p_1', 'p_2', 'p_3', 'p_4', 'p_5', 'p_6']

# branch_points = {
#   'b_1': ['p_2', 'p_5', 'p_6'],
#   'b_2': ['p_3', 'p_4', 'p_5'],}

# joints = {
#     'j_1' : ['p_1', None, 'p_2', None, 'p_3', None],
#     'j_2' : ['p_1', None, 'p_4', None, 'p_6', None]
# }

# homology_classes = {
#   'gamma_1' : ['p_1', 'p_2', 'p_3', 'p_6', 'p_4'],
#   'gamma_2' : ['p_5']}


# # --------------- N-punctured sphere -- 4 punctures -----------------

# streets = ['p_1', 'p_2', 'p_3', 'p_4', 'p_5', 'p_6']

# branch_points = {
#   'b_1': ['p_1', 'p_2', 'p_3'],
#   'b_2': ['p_2', 'p_4', 'p_5'],
#   'b_3': ['p_1', 'p_6', 'p_4'],
#   'b_4': ['p_3', 'p_5', 'p_6'],}

# joints = {}

# homology_classes = {
#   'gamma_1' : ['p_1'],
#   'gamma_2' : ['p_2'],
#   'gamma_3' : ['p_3'],
#   'gamma_4' : ['p_4'],
#   'gamma_5' : ['p_5'],
#   'gamma_6' : ['p_6'],}


# # --------------- T3 -----------------


# streets = ['p_1','p_2','p_3','p_4','p_5','p_6','p_7','p_8','p_9','p_10','p_11','p_12', 'p_13', 'p_14', 'p_15']

# branch_points = {
#   'b_1' : ['p_1', 'p_7', 'p_8'],
#   'b_2' : ['p_2', 'p_10', 'p_12'],
#   'b_3' : ['p_3', 'p_11', 'p_9'],
#   'b_4' : ['p_4','p_8','p_7'],
#   'b_5' : ['p_5', 'p_9', 'p_11'],
#   'b_6' : ['p_6', 'p_12', 'p_10']}

# joints = {
#   'j_1': ['p_1', 'p_15', 'p_2', 'p_13', 'p_3', 'p_14'],
#   'j_2': ['p_4', 'p_14', 'p_5', 'p_13', 'p_6', 'p_15'],}

# homology_classes = {
#   'gamma_1' : ['p_1', 'p_2', 'p_3'],
#   'gamma_2' : ['p_4', 'p_5', 'p_6'],
#   'gamma_3' : ['p_7'],
#   'gamma_4' : ['p_8'],
#   'gamma_5' : ['p_9'],
#   'gamma_6' : ['p_10'],
#   'gamma_7' : ['p_11'],
#   'gamma_8' : ['p_12'],
#   'gamma_9' : ['p_3', 'p_5', 'p_15'],
#   'gamma_10' : ['p_2', 'p_6', 'p_14'],
#   'gamma_11' : ['p_1', 'p_4', 'p_13'],}


# # --------------- [2,1]-punctured torus (FROM REDUCTION) -----------------

# streets = ['p_1', 'p_2', 'p_3', 'p_4', 'p_5', 'p_6', 'p_7', 'p_8', 'p_9']

# branch_points = {
#   'b_1': ['p_1', 'p_2', 'p_3'],
#   'b_2': ['p_1', 'p_4', 'p_8'],
#   'b_3': ['p_2', 'p_5', 'p_9'],
#   'b_4': ['p_3', 'p_6', 'p_7'],}

# joints = {
#     'j_1': ['p_4', None, 'p_5', None, 'p_6', None],
#     'j_2': ['p_7', None, 'p_8', None, 'p_9', None],
# }

# homology_classes = {
#   'gamma_1' : ['p_2'],
#   'gamma_2' : ['p_4', 'p_5', 'p_6'],
#   'gamma_3' : ['p_1'],
#   'gamma_4' : ['p_7', 'p_8', 'p_9'],
#   'gamma_5' : ['p_3'],
# }

# bp_coordinates = {
#   'b_1': [0.0, 0.0],
#   'b_2': [0.28, 0.0],
#   'b_3': [0.0, 0.31],
#   'b_4': [0.62, 0.57],
# }

# j_coordinates = {
#     'j_1': [0.32, 0.63],
#     'j_2': [0.62, 0.28],
# }

# e_coordinates = {
#     'p_1': [bp_coordinates['b_1'], bp_coordinates['b_2']],
#     'p_2': [bp_coordinates['b_1'], bp_coordinates['b_3']], 
#     'p_3': [bp_coordinates['b_4'], [bp_coordinates['b_1'][0] + 1, bp_coordinates['b_1'][1] + 1]], 
#     'p_4': [j_coordinates['j_1'], [bp_coordinates['b_2'][0], bp_coordinates['b_2'][1] + 1] ], 
#     'p_5': [j_coordinates['j_1'], bp_coordinates['b_3']],
#     'p_6': [j_coordinates['j_1'], bp_coordinates['b_4']],
#     'p_7': [j_coordinates['j_2'], bp_coordinates['b_4']],
#     'p_8': [j_coordinates['j_2'], bp_coordinates['b_2']],
#     'p_9': [j_coordinates['j_2'], [bp_coordinates['b_3'][0] + 1, bp_coordinates['b_3'][1]]],
# }

# one = ['p_1', 'p_8', 'p_9', 'p_2']
# tau = ['p_2', 'p_5', 'p_4', 'p_1']




# # --------------- [2,1]-punctured torus SELF-GLUING -----------------

# streets = ['p_1', 'p_2', 'p_3', 'p_4', 'p_5', 'p_6', 'p_7', 'p_8', 'p_9']

# branch_points = {
#   'b_1': ['p_1', 'p_2', 'p_3'],
#   'b_2': ['p_5', 'p_6', 'p_7'],
#   'b_3': ['p_2', 'p_4', 'p_8'],
#   'b_4': ['p_1', 'p_6', 'p_9'],}

# joints = {
#     'j_1': ['p_7', None, 'p_9', None, 'p_8', None],
#     'j_2': ['p_3', None, 'p_5', None, 'p_4', None],
# }

# homology_classes = {
#   'gamma_1' : ['p_1'],
#   'gamma_2' : ['p_2'],
#   'gamma_3' : ['p_3', 'p_4', 'p_5'],
#   'gamma_4' : ['p_7', 'p_8', 'p_9'],
#   'gamma_5' : ['p_6'],
# }

# bp_coordinates = {
#   'b_1': [0.711, 0.618],
#   'b_2': [0.512, 0.417],
#   'b_3': [0.613, 0.716],
#   'b_4': [0.414, 0.515],
# }

# j_coordinates = {
#     'j_1': [0.512, 0.715],
#     'j_2': [0.713, 0.514],
# }

# e_coordinates = {
#     'p_1': [bp_coordinates['b_1'], [bp_coordinates['b_4'][0] + 1, bp_coordinates['b_4'][1]]],
#     'p_2': [bp_coordinates['b_1'], bp_coordinates['b_3']], 
#     'p_3': [bp_coordinates['b_1'], j_coordinates['j_2']], 
#     'p_4': [j_coordinates['j_2'], [bp_coordinates['b_3'][0], bp_coordinates['b_3'][1] - 1] ], 
#     'p_5': [bp_coordinates['b_2'], j_coordinates['j_2']],
#     'p_6': [bp_coordinates['b_2'], bp_coordinates['b_4']],
#     'p_7': [j_coordinates['j_1'], [bp_coordinates['b_2'][0], bp_coordinates['b_2'][1] + 1]],
#     'p_8': [j_coordinates['j_1'], bp_coordinates['b_3']],
#     'p_9': [j_coordinates['j_1'], bp_coordinates['b_4']],
# }

# one = ['p_1', 'p_6', 'p_5', 'p_3']
# tau = ['p_7', 'p_6', 'p_9']




# # --------------- [3,1]-punctured torus FROM GUESSED REDUCTION (WORKS) -----------------

# streets = ['p_1', 'p_2', 'p_3', 'p_4', 'p_5', 'p_6', 'p_7', 'p_8', 'p_9', 'p_10', 'p_11', 'p_12', 'p_13', 'p_14', 'p_15']

# branch_points = {
#   'b_1': ['p_1', 'p_4', 'p_15'],
#   'b_2': ['p_4', 'p_11', 'p_5'],
#   'b_3': ['p_8', 'p_9', 'p_10'],
#   'b_4': ['p_8', 'p_7', 'p_3'],
#   'b_5': ['p_13', 'p_6', 'p_14'],
#   'b_6': ['p_13', 'p_12', 'p_2'],}

# joints = {
#     'j_1': ['p_1', None, 'p_3', None, 'p_2', None],
#     'j_2': ['p_10', None, 'p_11', None, 'p_12', None],
#     'j_3': ['p_9', None, 'p_15', None, 'p_14', None],
#     'j_4': ['p_5', None, 'p_6', None, 'p_7', None],
# }

# homology_classes = {
#   'gamma_1' : ['p_1', 'p_2', 'p_3'],
#   'gamma_2' : ['p_4'],
#   'gamma_3' : ['p_9', 'p_14', 'p_15'],
#   'gamma_4' : ['p_10', 'p_11', 'p_12'],
#   'gamma_5' : ['p_5', 'p_6', 'p_7'],
#   'gamma_6' : ['p_13'],
#   'gamma_7' : ['p_8'],
# }

# bp_coordinates = {
#   'b_1': [0.27, 0.23],
#   'b_2': [0.51, 0.45],
#   'b_3': [0.49, 0.0],
#   'b_4': [0.78, 0.0],
#   'b_5': [0.0, 0.61],
#   'b_6': [0.0, 0.7],
# }

# j_coordinates = {
#     'j_1': [0.0, 0.0],
#     'j_2': [0.24, 0.0],
#     'j_3': [0.0, 0.26],
#     'j_4': [0.73, 0.76],
# }

# e_coordinates = {
#     'p_1': [j_coordinates['j_1'], bp_coordinates['b_1']],
#     'p_2': [bp_coordinates['b_6'], [j_coordinates['j_1'][0], j_coordinates['j_1'][1] + 1]], 
#     'p_3': [bp_coordinates['b_4'], [j_coordinates['j_1'][0] + 1, j_coordinates['j_1'][1]]], 
#     'p_4': [bp_coordinates['b_1'], bp_coordinates['b_2']], 
#     'p_5': [j_coordinates['j_4'], bp_coordinates['b_2']],
#     'p_6': [j_coordinates['j_4'], [bp_coordinates['b_5'][0] + 1, bp_coordinates['b_5'][1]]],
#     'p_7': [j_coordinates['j_4'], [bp_coordinates['b_4'][0], bp_coordinates['b_4'][1] + 1]],
#     'p_8': [bp_coordinates['b_3'], bp_coordinates['b_4']],
#     'p_9': [bp_coordinates['b_3'], [j_coordinates['j_3'][0] + 1, j_coordinates['j_3'][1]]],
#     'p_10': [bp_coordinates['b_3'], j_coordinates['j_2']],
#     'p_11': [bp_coordinates['b_2'], j_coordinates['j_2']],
#     'p_12': [bp_coordinates['b_6'], [j_coordinates['j_2'][0], j_coordinates['j_2'][1] + 1]],
#     'p_13': [bp_coordinates['b_5'], bp_coordinates['b_6']],
#     'p_14': [bp_coordinates['b_5'], j_coordinates['j_3']],
#     'p_15': [bp_coordinates['b_1'], j_coordinates['j_3']],
# }

# one = ['p_1', 'p_4', 'p_11', 'p_10', 'p_8', 'p_3']
# tau = ['p_1', 'p_15', 'p_14', 'p_13', 'p_2']







# # --------------- [3,1]-punctured torus SELF-GLUING -----------------

# streets = ['p_1', 'p_2', 'p_3', 'p_4', 'p_5', 'p_6', 'p_7', 'p_8', 'p_9', 'p_10', 'p_11', 'p_12', 'p_13', 'p_14', 'p_15']

# branch_points = {
#   'b_1': ['p_1', 'p_2', 'p_3'],
#   'b_2': ['p_8', 'p_14', 'p_15'],
#   'b_3': ['p_5', 'p_13', 'p_12'],
#   'b_4': ['p_2', 'p_4', 'p_10'],
#   'b_5': ['p_7', 'p_8', 'p_9'],
#   'b_6': ['p_1', 'p_5', 'p_6'],}

# joints = {
#     'j_1': ['p_9', None, 'p_10', None, 'p_11', None],
#     'j_2': ['p_13', None, 'p_6', None, 'p_7', None],
#     'j_3': ['p_3', None, 'p_15', None, 'p_4', None],
#     'j_4': ['p_12', None, 'p_11', None, 'p_14', None],
# }

# homology_classes = {
#   'gamma_1' : ['p_1'],
#   'gamma_2' : ['p_2'],
#   'gamma_3' : ['p_3', 'p_4', 'p_15'],
#   'gamma_4' : ['p_9', 'p_10', 'p_11', 'p_12', 'p_14'],
#   'gamma_5' : ['p_13', 'p_7', 'p_6'],
#   'gamma_6' : ['p_5'],
#   'gamma_7' : ['p_8'],
# }

# bp_coordinates = {
#   'b_1': [0.87, 0.64],
#   'b_2': [0.62, 0.43],
#   'b_3': [0.35, 0.17],
#   'b_4': [0.66, 0.85],
#   'b_5': [0.43, 0.60],
#   'b_6': [0.18, 0.33],
# }

# j_coordinates = {
#     'j_1': [0.51, 0.78],
#     'j_2': [0.23, 0.52],
#     'j_3': [0.82, 0.47],
#     'j_4': [0.53, 0.21],
# }


# e_coordinates = {
#     'p_1': [bp_coordinates['b_1'], [bp_coordinates['b_6'][0] + 1, bp_coordinates['b_6'][1]]],
#     'p_2': [bp_coordinates['b_1'], bp_coordinates['b_4']],
#     'p_3': [bp_coordinates['b_1'], j_coordinates['j_3']],
#     'p_4': [bp_coordinates['b_4'], [j_coordinates['j_3'][0], j_coordinates['j_3'][1] + 1]],
#     'p_5': [bp_coordinates['b_3'], bp_coordinates['b_6']],
#     'p_6': [bp_coordinates['b_6'], j_coordinates['j_2']],
#     'p_7': [bp_coordinates['b_5'], j_coordinates['j_2']],
#     'p_8': [bp_coordinates['b_5'], bp_coordinates['b_2']],
#     'p_9': [bp_coordinates['b_5'], j_coordinates['j_1']],
#     'p_10': [bp_coordinates['b_4'], j_coordinates['j_1']],
#     'p_11': [j_coordinates['j_4'], [j_coordinates['j_1'][0], j_coordinates['j_1'][1] - 1]],
#     'p_12': [bp_coordinates['b_3'], j_coordinates['j_4']],
#     'p_13': [bp_coordinates['b_3'], [j_coordinates['j_2'][0], j_coordinates['j_2'][1] - 1]],
#     'p_14': [bp_coordinates['b_2'], j_coordinates['j_4']],
#     'p_15': [bp_coordinates['b_2'], j_coordinates['j_3']],
# }

# one = ['p_5', 'p_12', 'p_14', 'p_15', 'p_3', 'p_1']
# tau = ['p_11', 'p_14', 'p_8', 'p_9']






# # --------------- [3,1]-punctured torus SELF-GLUING (MODIFIED, NOT WORKING)-----------------

# streets2 = ['p_1', 'p_2', 'p_3', 'p_4', 'p_5', 'p_6', 'p_7', 'p_8', 'p_9', 'p_10', 'p_11', 'p_12', 'p_13', 'p_14', 'p_15']

# branch_points2 = {
#   'b_1': ['p_1', 'p_2', 'p_3'],
#   'b_2': ['p_8', 'p_14', 'p_15'],
#   'b_3': ['p_5', 'p_13', 'p_12'],
#   'b_4': ['p_3', 'p_4', 'p_6'],
#   'b_5': ['p_8', 'p_7', 'p_9'],
#   'b_6': ['p_10', 'p_11', 'p_13'],}

# joints2 = {
#     'j_1': ['p_2', None, 'p_10', None, 'p_9', None],
#     'j_2': ['p_12', None, 'p_7', None, 'p_6', None],
#     'j_3': ['p_14', None, 'p_5', None, 'p_4', None],
#     'j_4': ['p_15', None, 'p_11', None, 'p_1', None],
# }

# homology_classes2 = {
#   'gamma_1' : ['p_3'],
#   'gamma_2' : ['p_8'],
#   'gamma_3' : ['p_13'],
#   'gamma_4' : ['p_2', 'p_9', 'p_10'],
#   'gamma_5' : ['p_1', 'p_15', 'p_11'],
#   'gamma_6' : ['p_4', 'p_5', 'p_14'],
#   'gamma_7' : ['p_6', 'p_7', 'p_12'],
# }

# bp_coordinates2 = {
#   'b_1': [0.86, 0.66],
#   'b_2': [0.61, 0.41],
#   'b_3': [0.36, 0.16],
#   'b_4': [0.66, 0.86],
#   'b_5': [0.41, 0.61],
#   'b_6': [0.16, 0.36],
# }

# j_coordinates2 = {
#     'j_1': [0.5, 0.8],
#     'j_2': [0.2, 0.5],
#     'j_3': [0.8, 0.5],
#     'j_4': [0.5, 0.2],
# }


# e_coordinates2 = None
# # e_coordinates = {
# #     'p_1': [bp_coordinates['b_1'], j_coordinates['j_4']],
# #     'p_2': [bp_coordinates['b_1'], j_coordinates['j_1']],
# #     'p_3': [bp_coordinates['b_1'], bp_coordinates['b_4']],
# #     'p_4': [bp_coordinates['b_4'], [j_coordinates['j_3'][0], j_coordinates['j_3'][1] + 1]],
# #     'p_5': [bp_coordinates['b_3'], j_coordinates['j_3']],
# #     'p_6': [bp_coordinates['b_6'], j_coordinates['j_2']],
# #     'p_7': [bp_coordinates['b_5'], j_coordinates['j_2']],
# #     'p_8': [bp_coordinates['b_5'], bp_coordinates['b_2']],
# #     'p_9': [bp_coordinates['b_5'], j_coordinates['j_1']],
# #     'p_10': [bp_coordinates['b_4'], j_coordinates['j_1']],
# #     'p_11': [j_coordinates['j_4'], [j_coordinates['j_1'][0], j_coordinates['j_1'][1] - 1]],
# #     'p_12': [bp_coordinates['b_3'], j_coordinates['j_4']],
# #     'p_13': [bp_coordinates['b_3'], [j_coordinates['j_2'][0], j_coordinates['j_2'][1] - 1]],
# #     'p_14': [bp_coordinates['b_2'], j_coordinates['j_4']],
# #     'p_15': [bp_coordinates['b_2'], j_coordinates['j_3']],
# # }

# one2 = ['p_5', 'p_12', 'p_14', 'p_15', 'p_3', 'p_1']
# tau2 = ['p_11', 'p_14', 'p_8', 'p_9']









# --------------- [4,1]-punctured torus SELF-GLUING -----------------

streets = ['p_' + str(i + 1) for i in range(21)]

branch_points = {
  'b_1': ['p_1', 'p_2', 'p_3'],
  'b_2': ['p_5', 'p_6', 'p_7'],
  'b_3': ['p_8', 'p_15', 'p_16'],
  'b_4': ['p_17', 'p_18', 'p_19'],
  'b_5': ['p_2', 'p_4', 'p_11'],
  'b_6': ['p_6', 'p_10', 'p_12'],
  'b_7': ['p_15', 'p_14', 'p_21'],
  'b_8': ['p_18', 'p_20', 'p_1'],}

joints = {
    'j_1': ['p_9', None, 'p_10', None, 'p_11', None],
    'j_2': ['p_13', None, 'p_14', None, 'p_12', None],
    'j_3': ['p_19', None, 'p_20', None, 'p_21', None],
    'j_4': ['p_3', None, 'p_5', None, 'p_4', None],
    'j_5': ['p_7', None, 'p_8', None, 'p_9', None],
    'j_6': ['p_13', None, 'p_16', None, 'p_17', None],
}

homology_classes = {
  'gamma_1' : ['p_1'],
  'gamma_2' : ['p_2'],
  'gamma_3' : ['p_3', 'p_4', 'p_5'],
  'gamma_4' : ['p_6'],
  'gamma_5' : ['p_7', 'p_8', 'p_9', 'p_10', 'p_11'],
  'gamma_6' : ['p_12', 'p_13', 'p_14', 'p_16', 'p_17'],
  'gamma_7' : ['p_15'],
  'gamma_8' : ['p_19', 'p_20', 'p_21'],
  'gamma_9' : ['p_18'],
}

bp_coordinates = {
  'b_1': [0.81, 0.61],
  'b_2': [0.71, 0.51],
  'b_3': [0.61, 0.41],
  'b_4': [0.51, 0.31],
  'b_5': [0.61, 0.81],
  'b_6': [0.51, 0.71],
  'b_7': [0.41, 0.61],
  'b_8': [0.31, 0.51],
}

j_coordinates = {
    'j_1': [0.5, 0.8],
    'j_2': [0.4, 0.7],
    'j_3': [0.3, 0.6],
    'j_4': [0.8, 0.5],
    'j_5': [0.7, 0.4],
    'j_6': [0.6, 0.3],
}


e_coordinates = {
    'p_1': [bp_coordinates['b_1'], [bp_coordinates['b_8'][0] + 1, bp_coordinates['b_8'][1]]],
    'p_2': [bp_coordinates['b_1'], bp_coordinates['b_5']],
    'p_3': [bp_coordinates['b_1'], j_coordinates['j_4']],
    'p_4': [bp_coordinates['b_5'], [j_coordinates['j_4'][0], j_coordinates['j_4'][1] + 1]],
    'p_5': [bp_coordinates['b_2'], j_coordinates['j_4']],
    'p_6': [bp_coordinates['b_6'], bp_coordinates['b_2']],
    'p_7': [bp_coordinates['b_2'], j_coordinates['j_5']],
    'p_8': [bp_coordinates['b_3'], j_coordinates['j_5']],
    'p_9': [j_coordinates['j_5'], [j_coordinates['j_1'][0], j_coordinates['j_1'][1] - 1]],
    'p_10': [bp_coordinates['b_6'], j_coordinates['j_1']],
    'p_11': [j_coordinates['j_1'], bp_coordinates['b_5']],
    'p_12': [bp_coordinates['b_6'], j_coordinates['j_2']],
    'p_13': [j_coordinates['j_2'], [j_coordinates['j_6'][0], j_coordinates['j_6'][1] + 1]],
    'p_14': [bp_coordinates['b_7'], j_coordinates['j_2']],
    'p_15': [bp_coordinates['b_3'], bp_coordinates['b_7']],
    'p_16': [bp_coordinates['b_3'], j_coordinates['j_6']],
    'p_17': [bp_coordinates['b_4'], j_coordinates['j_6']],
    'p_18': [bp_coordinates['b_4'], bp_coordinates['b_8']],
    'p_19': [bp_coordinates['b_4'], [j_coordinates['j_3'][0], j_coordinates['j_3'][1] - 1]],
    'p_20': [bp_coordinates['b_8'], j_coordinates['j_3']],
    'p_21': [bp_coordinates['b_7'], j_coordinates['j_3']],
}

one = ['p_1', 'p_20', 'p_21', 'p_14', 'p_12', 'p_10', 'p_11', 'p_2']
tau = ['p_19', 'p_18', 'p_20']










w = BPSgraph(
    branch_points=branch_points, 
    streets=streets, 
    joints=joints, 
    homology_classes=homology_classes,
    faces=None,
    bp_coordinates=bp_coordinates,
    j_coordinates=j_coordinates,
    e_coordinates=e_coordinates,
)




####################################################
####################################################

####################################################
####################################################


# analyze all sequences which give back the BPS graph
max_n_moves = 5
last_moves_to_avoid = 1
last_mutations_to_avoid = 1
SAVE_PLOTS = False
seq = find_invariant_sequences(
    w, max_n_moves, level=0, ref_graph=w, 
    edge_cycle_min_length=0,
    min_n_cooties=0,
    fundamental_region=[one, tau],
    drop_if_trivial_perm=False,
    avoid_last_n_moves=last_moves_to_avoid,
    avoid_last_n_mutations=last_mutations_to_avoid,
)

# print 'Found {} sequences.'.format(len(seq))

# old_modular_parameter = compute_modular_parameter(one, tau, w)
# old_one = sum_up_edges(one, w)[0] + 1j * sum_up_edges(one, w)[1]
# old_tau = sum_up_edges(tau, w)[0] + 1j * sum_up_edges(tau, w)[1]
# # #
# # S_old_modular_parameter = -1 / old_modular_parameter
# # S_old_one = old_tau 
# # S_old_tau = -1 * old_one 
# # S_inv_old_one = -1 * old_tau 
# # S_inv_old_tau = old_one
# # #
# # L_old_modular_parameter = old_modular_parameter / (1 + old_modular_parameter)
# # L_old_one = old_one + old_tau
# # L_old_tau = old_tau
# # L_inv_old_one = old_one - old_tau
# # L_inv_old_tau = old_tau
# # #
# # R_old_modular_parameter = 1 + old_modular_parameter
# # R_old_one = old_one 
# # R_old_tau = old_tau + old_one
# # R_inv_old_one = old_one 
# # R_inv_old_tau = old_tau - old_one

# # prepare a directory
# mydir = os.path.join(
#     os.getcwd(), 'mcg_moves', 
#     (datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S') + '_{}_moves'.format(max_n_moves))
# )
# os.makedirs(mydir)

# # save info about sequences of moves in a text file
# text_file = open(mydir + '/sequence_data.txt', 'w')
# text_file.write('\t\tSequence data\n\n')
# text_file.write(
#     'Avoiding the last {} moves and last {} mutations.\n\n'.format(
#         last_moves_to_avoid, last_mutations_to_avoid
#     )
# )
# text_file.write(
#     'Modular parameter of the original torus: {}'
#     '\none : {} = {}\ntau : {} = {}\n\n'.format(
#         old_modular_parameter, 
#         one, sum_up_edges(one, w), 
#         tau, sum_up_edges(tau, w)
#     )
# )
# text_file.write('Quiver:\n{}\n\n'.format(w.seed.quiver))

# # find permutations which leave the graph invariant
# text_file.write(
#     '\n----------------------------------\n'
#     'Nontrivial permutations giving isomorphic graphs\n'
# )
# self_permutations = find_self_permutations(w, [one, tau])
# for p, mod_param in self_permutations:
#     hc_perm_dic = {}
#     for p_k, p_v in p.iteritems():
#         source = (
#             [h_k for h_k, h_v in w.basis_homology_classes.iteritems()
#             if p_k in [s.label for s in h_v.streets]][0]
#         )
#         target = (
#             [h_k for h_k, h_v in w.basis_homology_classes.iteritems()
#             if p_v in [s.label for s in h_v.streets]][0]
#         )
#         hc_perm_dic[source] = target

#     text_file.write('\t{}\n'.format(p))
#     text_file.write('\tmodular parameter : {}\n'.format(mod_param))
#     text_file.write('\thomology class permutation : {}\n'.format(hc_perm_dic))
    
# # S_move_candidates = []
# # L_move_candidates = []
# # R_move_candidates = []
# # S_inv_move_candidates = []
# # L_inv_move_candidates = []
# # R_inv_move_candidates = []

# # Keep track of how many sequences give the same quiver mutation sequence
# # only give one sequence of moves that corresponds to a given mutation
# # sequence.
# seen_by_mutations = {}
# # keep track of all possible modular parameters associated to a sequence 
# # of mutations
# mutation_seq_modular_parameters = {}
# # keep track of all modular parameters which are seen
# all_modular_parameters = []

# # keep track of modular transformations matrices (better than modular parameters)
# all_modular_transformations = []
# # keep track of all possible mutations that give the same modular transformation
# modular_transf_mutations = {}

# for i, s in enumerate(seq):
#     if str(s.quiver_mutations) not in seen_by_mutations.keys():
#         # initiate the count for graph moves corresponding to 
#         # this mutation sequence.
#         # use a list of mutation turned into a string for dictionary keys...
#         seen_by_mutations[str(s.quiver_mutations)] = [i]

#         text_file.write(
#             '\n\n-------------------------'
#             '\nSequence #{} : {}\nMutations : {}'.format(
#                 i, s.moves_sequence, s.quiver_mutations
#             )
#         )
#         text_file.write(
#             '\n\t\tPermutation data\n(Note: this is in the form '
#                 '(.., old_edge : new_edge, ..)'
#                 '\nAll permutations are included\n'
#         )
#         if SAVE_PLOTS is True:
#             sub_dir = os.path.join(mydir, 'sequence_'+str(i))
#             os.makedirs(sub_dir)
#             new_graph = apply_sequence(w, s.moves_sequence, save_plot=sub_dir)
#         else:
#             new_graph = apply_sequence(w, s.moves_sequence)
#         # Each sequence of moves comes with a set of permutations
#         # that relate the final graph to the original one.
#         # Must consider each of them
#         for j in range(len(s.edge_permutations)):
#             # the edge permutation
#             p = s.edge_permutations[j] 
#             ### DEBUG: disabling this for now
#             # # the homology class permutation
#             # hc_p = 'homology permutation not implemented'
#             hc_p = s.homology_permutations[j] 
#             ### DEBUG END

#             # use the infor about the funadmental region to compute the 
#             # new modular parameter, by tracking where the edges 
#             # ended up.
#             # Before the moves, the two parameteters [1, \tau] were 
#             # specified by a sequence of edges each.
#             # Now tracking where those edges ended up we are
#             # able to compute the new [1', \tau'].
#             # First of all, we translate the old sequence of edges
#             # into the new sequence, by applying the permutation 
#             # dictionary
#             new_one = [p[edge] for edge in one]
#             new_tau = [p[edge] for edge in tau]
#             new_modular_parameter = compute_modular_parameter(
#                 new_one, new_tau, new_graph
#             )
#             new_one_value = (
#                 sum_up_edges(new_one, new_graph, int_round=True)[0] + 
#                 1j * sum_up_edges(new_one, new_graph, int_round=True)[1]
#             )
#             new_tau_value = (
#                 sum_up_edges(new_tau, new_graph, int_round=True)[0] +
#                 1j * sum_up_edges(new_tau, new_graph, int_round=True)[1]
#             )
#             # record the new modular parameter produced by this 
#             # sequence of moves
#             mutation_seq_modular_parameters[str(s.quiver_mutations)] = (
#                 [new_modular_parameter]
#             )
#             # also record the modular parameter in overall list
#             if c_num_is_in_list(
#                 new_modular_parameter, 
#                 all_modular_parameters,
#                 EPSILON
#             ) is False:
#                 all_modular_parameters.append(new_modular_parameter)

#             # also record the modular transformation in overall list,
#             # and record which mutations give the same modular transformatilon
#             modular_tmn = [
#                     [
#                         sum_up_edges(new_tau, new_graph, int_round=True)[1], 
#                         sum_up_edges(new_tau, new_graph, int_round=True)[0]
#                     ], 
#                     [
#                         sum_up_edges(new_one, new_graph, int_round=True)[1], 
#                         sum_up_edges(new_one, new_graph, int_round=True)[0]
#                     ]
#                 ]

#             if str(modular_tmn) not in all_modular_transformations:
#                 all_modular_transformations.append(str(modular_tmn))
#                 modular_transf_mutations[str(modular_tmn)] = (
#                     [str(s.quiver_mutations)]
#                 )
#             else: 
#                 modular_transf_mutations[str(modular_tmn)].append(
#                     str(s.quiver_mutations)
#                 )


#             # # check if this is a candidate for n S, L or R move:
#             # #
#             # # if abs(new_modular_parameter - S_old_modular_parameter) < EPSILON:
#             # #     S_move_candidates.append([i, j])
#             # # if abs(new_modular_parameter - L_old_modular_parameter) < EPSILON:
#             # #     L_move_candidates.append([i, j])
#             # # if abs(new_modular_parameter - R_old_modular_parameter) < EPSILON:
#             # #     R_move_candidates.append([i, j])
#             # if (
#             #     abs(new_one_value - S_old_one) < EPSILON and 
#             #     abs(new_tau_value - S_old_tau) < EPSILON
#             # ):
#             #     S_move_candidates.append([i, j])
#             # if (
#             #     abs(new_one_value - S_inv_old_one) < EPSILON and 
#             #     abs(new_tau_value - S_inv_old_tau) < EPSILON
#             # ):
#             #     S_inv_move_candidates.append([i, j])
#             # if (
#             #     abs(new_one_value - L_old_one) < EPSILON and 
#             #     abs(new_tau_value - L_old_tau) < EPSILON
#             # ):
#             #     L_move_candidates.append([i, j])
#             # if (
#             #     abs(new_one_value - L_inv_old_one) < EPSILON and 
#             #     abs(new_tau_value - L_inv_old_tau) < EPSILON
#             # ):
#             #     L_inv_move_candidates.append([i, j])
#             # if (
#             #     abs(new_one_value - R_old_one) < EPSILON and 
#             #     abs(new_tau_value - R_old_tau) < EPSILON
#             # ):
#             #     R_move_candidates.append([i, j])
#             # if (
#             #     abs(new_one_value - R_inv_old_one) < EPSILON and 
#             #     abs(new_tau_value - R_inv_old_tau) < EPSILON
#             # ):
#             #     R_inv_move_candidates.append([i, j])

#             text_file.write(
#                 '\n\tpermutation #{}'
#                 '\n\t---------------'
#                 '\n\ton edges:\n\t{}'
#                 '\n\ton homology basis:\n\t{}'
#                 '\n\tmodular parameter : {}'.format(
#                     j, p, hc_p, new_modular_parameter
#                 )
#             )
#             text_file.write(
#                 '\n\tThe new one : {} = {}\n\tThe new tau : {} = {}'.format(
#                     sum_up_edges(new_one, new_graph, int_round=True), new_one, 
#                     sum_up_edges(new_tau, new_graph, int_round=True), new_tau
#                 )
#             )
#             text_file.write(
#                 '\n\tModular transformation: \n{}\n'.format(modular_tmn)
#             )
            
#     else:
#         # add the sequence to the list of moves sequences that give
#         # a certain mutation sequence
#         seen_by_mutations[str(s.quiver_mutations)].append(i)
#         # record which modular parameter is produced by this sequence of moves
#         new_graph = apply_sequence(w, s.moves_sequence)
        
#         for j in range(len(s.edge_permutations)):
#             # the edge permutation
#             p = s.edge_permutations[j] 
#             ### DEBUG: disabling this for now
#             # # the homology class permutation
#             # hc_p = 'homology permutation not implemented'
#             hc_p = s.homology_permutations[j] 
#             ### DEBUG END

#             # use the infor about the funadmental region to compute the 
#             # new modular parameter, by tracking where the edges 
#             # ended up.
#             # Before the moves, the two parameteters [1, \tau] were 
#             # specified by a sequence of edges each.
#             # Now tracking where those edges ended up we are
#             # able to compute the new [1', \tau'].
#             # First of all, we translate the old sequence of edges
#             # into the new sequence, by applying the permutation 
#             # dictionary
#             new_one = [p[edge] for edge in one]
#             new_tau = [p[edge] for edge in tau]
#             new_modular_parameter = compute_modular_parameter(
#                 new_one, new_tau, new_graph
#             )
#             # record the new modular parameter produced by this 
#             # sequqnece of moves
#             if (
#                 new_modular_parameter not in 
#                 mutation_seq_modular_parameters[str(s.quiver_mutations)]
#             ):
#                 mutation_seq_modular_parameters[str(s.quiver_mutations)].append(
#                     new_modular_parameter
#                 )

# # Write the tally of how many moves seem to give the same sequence of 
# # mutations
# text_file.write('\n\n\n\tTally of mutations\n\n')
# for mut_seq in seen_by_mutations.keys():
#     text_file.write(
#         'This mutation sequence : {}\nwas seen {} times. '
#         'In moves sequences : {}\n'
#         .format(
#             mut_seq, len(seen_by_mutations[mut_seq]), 
#             seen_by_mutations[mut_seq]
#         )
#     )
#     text_file.write(
#         'Associated modular parameters: {}\n\n'.format(
#             [round(mp.real, 2) + 1j * round(mp.imag, 2) 
#             for mp in mutation_seq_modular_parameters[mut_seq]]
#         )
#     )


# # Write the tally of which mutations seem to give the same modular transformation
# # avoiding the identity transformation
# text_file.write(
#     '\n\n\tTally of modular transformations (skipping the identity)'
# )
# for mod_transf in modular_transf_mutations.keys():
#     if (
#         (eval(mod_transf) != [[1, 0], [0, 1]]) and 
#         (eval(mod_transf) != [[-1, 0], [0, -1]])
#     ):
#         text_file.write(
#             '\n\nThis modular transformation : \n{}\n'
#             'Corresponds to the following mutations\n'
#             .format(numpy.array(eval(mod_transf)))
#         )
#         for mut in sorted(modular_transf_mutations[mod_transf], key=len):
#             text_file.write('{}\n'.format(mut))


# # # write in the text file the candidates for S, L, R transformations,
# # # and their inverses
# # text_file.write('\n\n-----------------------\n\t Candidates for S-move\n')
# # for [i, j] in S_move_candidates:
# #     text_file.write('\nmove #{}, permutation #{}'.format(i, j))

# # text_file.write('\n\n-----------------------\n\t Candidates for S^-1-move\n')
# # for [i, j] in S_inv_move_candidates:
# #     text_file.write('\nmove #{}, permutation #{}'.format(i, j))

# # text_file.write('\n\n-----------------------\n\t Candidates for L-move\n')
# # for [i, j] in L_move_candidates:
# #     text_file.write('\nmove #{}, permutation #{}'.format(i, j))

# # text_file.write('\n\n-----------------------\n\t Candidates for L^-1-move\n')
# # for [i, j] in L_inv_move_candidates:
# #     text_file.write('\nmove #{}, permutation #{}'.format(i, j))

# # text_file.write('\n\n-----------------------\n\t Candidates for R-move\n')
# # for [i, j] in R_move_candidates:
# #     text_file.write('\nmove #{}, permutation #{}'.format(i, j))

# # text_file.write('\n\n-----------------------\n\t Candidates for R^-1-move\n')
# # for [i, j] in R_inv_move_candidates:
# #     text_file.write('\nmove #{}, permutation #{}'.format(i, j))


# text_file.write('\n\n-----------------------\n\t All modular transformations\n')
# for mt in all_modular_transformations:
#     text_file.write('\n{}\n'.format(numpy.array(eval(mt))))



# ####################################################
# ####################################################

# ####################################################
# ####################################################



# ### STUDY A PARTICULAR SEQUENCE
# seq_0 = ['p_13', 'p_8', 'f_3', 'p_13', 'p_2', 'f_2', 'p_14', 'f_2', 'p_2', 'p_14', 'p_2', 'p_13']
# seq_1 = ['p_4', 'p_8', 'f_2', 'p_3', 'f_4']
# seq_2 = ['p_13', 'f_4']
# seq_3 = ['f_0', 'p_9', 'p_3']
# seq_4 = ['p_6', 'p_2', 'f_1']
# seq_5 = ['p_2', 'p_5', 'p_8', 'f_1', 'f_3', 'p_8']
# seq_6 = seq_5 + ['p_1', 'p_4', 'p_1', 'p_4', 'p_1']

# ### THIS IS A GOOD SEQUENCE for the [3,1] graph derived by reduction of the puncture (not by self-gluing)
# ### it should correspond to the SL2Z element (1, 0 // -1, 1) = L^{-1}
# seq_7 = [
# 'p_13', 'p_8', 'f_4', 'p_2', 'p_14', 
# 'p_13', 'f_3', 'p_9', 'p_5', 'p_14', 
# 'p_13', 'p_10', 'f_0',
# 'p_14', 'p_3', 'f_4', 'p_4'
# ]

# #### THIS GIVES THE S-transformation! For the [3,1] graph obtainedd by reduction
# seq_8 = ['p_4', 'f_2', 'p_10', 'p_8', 'p_3', 'p_10', 'f_4', 'p_11', 'p_5', 'p_3',
# 'p_8', 'p_14', 'f_0',
# 'p_1', 'f_2', 'p_13']

# #### This takes the [3,1] selfglued graph and transforms it into the [3,1] graph obtained by reduction of the puncture
# seq_9 = ['p_11', 'f_0', 'p_1', 'p_5', 'p_3', 'f_4', 'p_7']

# #### The L-move for the [3,1] selfglued graph
# seq_10 = ['p_11', 'f_0', 'f_2', 'p_11', 'p_6', 'p_3']

# #### The (S L^{-1}) -move for the [3,1] selfglued graph
# ###seq_11=['p_1', 'p_5', 'p_2', 'f_4', 'p_15', 'p_7', 'p_11', 'f_0', 'f_2', 'p_8', 'p_1']
# seq_11 = ['p_11', 'p_1', 'p_5', 'p_2', 'f_4', 'p_15', 'f_0', 'p_7', 'f_2', 'p_8', 'p_1']

# #### L-move sequence for [4,1] selfglued graph
# seq_12 = ['p_9', 'p_13'] + ['f_4', 'f_1', 'f_2'] + ['p_9', 'p_13'] + ['p_20', 'p_3']

# #### The (S L^{-1}) -move for the [4,1] selfglued graph NOT WORKING(?)
# # seq_13 = (
# #     ['p_1', 'p_18', 'p_2', 'f_6', 'p_21', 'p_5', 'p_9', 'p_13'] +
# #     ['f_4', 'f_1', 'f_2'] + ['p_1'] #+ ['f_1'] + ['p_20', 'p_3'] + ['f_6']
# #     #+ ['p_7', 'p_14'] + ['p_15', 'p_6'] + ['p_1']
# # )

# #### An absolutely crazy sequence for the [4,1] selfglued graph
# # seq_13 = (
# #     ['p_1', 'p_18', 'p_2', 'f_6', 'p_9', 'p_13', 'f_1', 'p_5', 'p_21'] +
# #     ['f_4', 'f_2', 'f_1', 'p_6', 'p_15', 'p_1', 'p_3', 'p_20', 'f_5'] + 
# #     ['f_0', 'p_14', 'p_7', 'p_8', 'p_12', 'f_3', 'p_6', 'p_15', 'p_1'] +
# #     [ 'f_1', 'p_18', 'p_2', 'f_3', 'p_14', 'p_7']#, 'p_1']
# # )

# #### trying to find another SL(2,Z) move for the [4,1] selfglued graph, inspired by SL^-1 move for [3,1]
# seq_14 = ['p_13', 'p_9', 'p_1', 'p_18', 'p_2', 'f_6', 'p_5', 'f_2', 'p_21', 'f_4', 'p_6', 'p_15', 'p_1']
# seq_15 = ['p_13', 'p_9', 'p_1', 'p_18', 'p_2', 'f_6', 'p_5', 'p_21']


# seq_16 = ['p_13', 'f_4', 'p_9', 'f_4', 'f_1', 'f_2', 'f_1', 'p_1', 'p_18', 'f_2', 'p_1', 'p_13', 'p_9', 'p_18']#, 'p_1']



# # mydir = os.path.join(
# #     os.getcwd(), 'mcg_moves', 
# #     (datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S') + '_DEBUG')
# # )
# # # os.makedirs(mydir)

# w1, mutation_seq = apply_sequence(
#     w, seq_17, 
#     save_plot=None,
#     # save_plot=mydir, 
#     include_mutation_sequence=True
# )

# # # # print w1.compute_quiver()

# # print w1.mutable_faces
# # print w1.mutable_edges

# # perms = are_equivalent_as_graphs(w, w1)
# # # print len(perms)
# # print perms

# # # print '\n\n'
# # # w1.faces['f_1'].print_info()
# # # print '\n\n'
# # # w1.faces['f_2'].print_info()
# # # print '\n\n'
# # # w1.faces['f_4'].print_info()


# ##

# find_sequence_completion(
#     w, seq_16, 1, one, tau, save_files=True, 
#     drop_if_trivial_perm=False,
#     avoid_last_n_moves=None, avoid_last_n_mutations=1,
# )







##############


########### Build a partial sequence based on some criterion
# import random
# def build_partial_sequence(graph, length):
#     """
#     Building a partial sequence based on certain criteria
#     """
#     new_graph = graph
#     seq = []
#     for i in range(length):
#         forbidden_moves = seq[-4:]
#         moves = new_graph.mutable_faces + new_graph.mutable_edges
#         allowed_moves = [m for m in moves if m not in forbidden_moves]
#         new_move = random.choice(allowed_moves)
#         seq.append(new_move)
#         new_graph = apply_sequence(
#             new_graph, [new_move], 
#             save_plot=None,
#             # save_plot=mydir, 
#             include_mutation_sequence=False
#         )
#     return seq

# for i in range(30):
#     partial_seq = build_partial_sequence(w, 18)
#     print 'The partial sequence is: {}'.format(
#             partial_seq
#     )

#     partial_seq_comps = find_sequence_completion(
#         w, partial_seq, 5, one, tau, save_files=True,
#         drop_if_trivial_perm=False,
#         avoid_last_n_moves=None, avoid_last_n_mutations=None,
#     )

########### Build a partial sequence based on some criterion,
########### and incrementally look for completions of the sequence
# import random
# def venture(graph, in_length, fin_length, completion_depth):
#     """
#     Building a partial sequence based on certain criteria 
#     and checking at each step (between in_length, fin_length)
#     whether there is a completion (with up to 'completion_depth'
#     number of moves).
#     """
#     new_graph = graph
#     seq = []
#     for i in range(fin_length):
#         forbidden_moves = seq[-4:]
#         moves = new_graph.mutable_faces + new_graph.mutable_edges
#         allowed_moves = [m for m in moves if m not in forbidden_moves]
#         if len(allowed_moves) == 0:
#             allowed_moves = moves

#         new_move = random.choice(allowed_moves)
#         seq.append(new_move)
#         new_graph = apply_sequence(
#             new_graph, [new_move], 
#             save_plot=None,
#             # save_plot=mydir, 
#             include_mutation_sequence=False
#         )

#         if i >= in_length:
#             partial_seq_comps = find_sequence_completion(
#                 w, seq, completion_depth, one, tau, 
#                 save_files='venture_length_{}'.format(i),
#                 drop_if_trivial_perm=True,
#                 avoid_last_n_moves=None, avoid_last_n_mutations=None,
#             )


# min_seq_length = 15
# max_seq_length = 25
# completion_depth = 5

# n_tries = 30
# for i in range(n_tries):
#     venture(w, min_seq_length, max_seq_length, completion_depth)



##############################
##############################
#   Find a sequence leading to a Z3 symmetric graph

# max_n_moves = 10
# seq = find_graphs_with_symmetry(w, max_n_moves)


# # prepare a directory
# mydir = os.path.join(
#     os.getcwd(), 'mcg_moves', 
#     (datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S') + 'symm_seq__{}_moves'.format(max_n_moves))
# )
# os.makedirs(mydir)

# # save info about sequences of moves in a text file
# text_file = open(mydir + '/symmetric_sequence_data.txt', 'w')
# text_file.write('\t\tSymmetric sequence data\n\n')
# text_file.write(
#     'Number of streets of the original graph: {}'
#     .format(len(streets))
# )
# text_file.write('\n\nSymmetric sequences:\n')
# for s in seq:
#     text_file.write('\n{}'.format(s))



##############################

###########

# ###### OLD-DEPRECATED Try to find a completion of the sequence #####

# max_n_moves = 11
# last_moves_to_avoid = 5
# seq = find_invariant_sequences(
#     w1, max_n_moves, level=0, ref_graph=w, 
#     edge_cycle_min_length=0,
#     min_n_cooties=0,
#     fundamental_region=[one, tau],
#     drop_if_trivial_perm=True,
#     avoid_last_n_moves=last_moves_to_avoid,
#     avoid_last_n_mutations=last_mutations_to_avoid,
# )
# print 'Found {} sequences.'.format(len(seq))
# for s in seq:
#     s.print_info()

######## Analyze the sequence: give quiver permutation and SL2Z transf.

# print 'Sequence\n{}\ncorresponds to the mutations\n{}\n'.format(
#     seq_12, mutation_seq
# )

# perms = are_equivalent_as_graphs(w, w1)
# if perms is not None:
#     print '\nFound the following permutations'

#     for p in perms:
#         print p

#     for j, p in enumerate(perms):
#         new_one = [p[edge] for edge in one]
#         new_tau = [p[edge] for edge in tau]
#         new_modular_parameter = compute_modular_parameter(
#             new_one, new_tau, w1
#         )
#         print (
#             '\nModular parameter of this sequence '
#             'with permutation #{}: {}\n'
#             .format(j, new_modular_parameter)
#         )
#         print (
#             'The new one : {} = {}\nThe new tau : {} = {}'.format(
#                 new_one, sum_up_edges(new_one, w1), 
#                 new_tau, sum_up_edges(new_tau, w1)
#             )
#         )
#         print (
#             'The permutation of homology classes \n{}'.format(
#                 quiver_node_permutation(w, w1, p) 
#             )
#         )
# else:
#     print 'found no permutations'

# #######


# for p in perms:
#     print p
#     print determine_perm_cycles(p.keys(), p.values())

# print quiver_node_permutation(w, w1 ,
#     perms[0]
#     )
# print quiver_node_permutation(w, w1, 
#     perms[1]
#     )
### DEBUG END


# self_permutations = find_self_permutations(w, [one, tau])
# print self_permutations
# for p, mod_param in self_permutations:
#     print 'The permutation of homology classes for \n{}'.format(p)
#     hc_perm_dic = {}
#     for p_k, p_v in p.iteritems():
#         source = (
#             [h_k for h_k, h_v in w.basis_homology_classes.iteritems()
#             if p_k in [s.label for s in h_v.streets]][0]
#         )
#         target = (
#             [h_k for h_k, h_v in w.basis_homology_classes.iteritems()
#             if p_v in [s.label for s in h_v.streets]][0]
#         )
#         hc_perm_dic[source] = target
#     print hc_perm_dic
    
    




# ####################################################
# ###     Just find permutations which leave a graph invariant 
# w2 = apply_sequence(w, ['p_13', 'f_4', 'p_4', 'p_5', 'p_4', 'p_5', 'p_4', 'f_4', 'p_13'])
# perms = are_equivalent_as_graphs(w, w2)
# print '\nFound {} permutations'.format(len(perms))
# for p in perms:
#     print p
#     print determine_perm_cycles(p.keys(), p.values())

# # streets = w.streets.keys()
# # # s_0 = streets[0]
# # s_0 = 'p_14'
# # epts = [ep.end_point.label for ep in w.streets[s_0].endpoints]

# # print 'edges of graph w'
# # print streets
# # print 'endpoints of {} are {}'.format(s_0, epts)

# # print 
# # perm = match_graphs_from_starting_point(
# #     w, w, s_0, -1
# # )
# # print perm


# # --------------- [3,1]-punctured torus SELF-GLUING 3 twist cable (NOT WORKING)-----------------

# streets_3 = ['p_1', 'p_2', 'p_3', 'p_4', 'p_5', 'p_6', 'p_7', 'p_8', 'p_9', 'p_10', 'p_11', 'p_12', 'p_13', 'p_14', 'p_15']

# branch_points_3 = {
#   'b_1': ['p_1', 'p_2', 'p_3'],
#   'b_2': ['p_8', 'p_14', 'p_15'],
#   'b_3': ['p_5', 'p_13', 'p_12'],
#   'b_4': ['p_5', 'p_1', 'p_10'],
#   'b_5': ['p_7', 'p_8', 'p_9'],
#   'b_6': ['p_2', 'p_6', 'p_13'],}

# joints_3 = {
#     'j_1': ['p_9', None, 'p_10', None, 'p_4', None],
#     'j_2': ['p_11', None, 'p_6', None, 'p_7', None],
#     'j_3': ['p_3', None, 'p_15', None, 'p_4', None],
#     'j_4': ['p_12', None, 'p_11', None, 'p_14', None],
# }

# homology_classes_3 = {
#   'gamma_1' : ['p_1'],
#   'gamma_2' : ['p_2'],
#   'gamma_3' : ['p_3', 'p_4', 'p_15','p_9', 'p_10'],
#   'gamma_4' : ['p_7', 'p_6', 'p_11', 'p_12', 'p_14'],
#   'gamma_5' : ['p_13'],
#   'gamma_6' : ['p_5'],
#   'gamma_7' : ['p_8'],
# }

# bp_coordinates_3 = {
#   'b_1': [0.27, 0.23],
#   'b_2': [0.51, 0.45],
#   'b_3': [0.49, 0.0],
#   'b_4': [0.78, 0.0],
#   'b_5': [0.0, 0.61],
#   'b_6': [0.0, 0.7],
# }

# j_coordinates_3 = {
#     'j_1': [0.0, 0.0],
#     'j_2': [0.24, 0.0],
#     'j_3': [0.0, 0.26],
#     'j_4': [0.73, 0.76],
# }

# ### NOTE: unreliable tau etc
# e_coordinates_3 = None

# one_3 = ['p_1', 'p_10', 'p_4' 'p_3']
# tau_3 = ['p_1', 'p_10', 'p_4' 'p_3']




# # --------------- [3,1]-punctured torus FROM GUESSED REDUCTION (not working) -----------------

# streets = ['p_1', 'p_2', 'p_3', 'p_4', 'p_5', 'p_6', 'p_7', 'p_8', 'p_9', 'p_10', 'p_11', 'p_12', 'p_13', 'p_14', 'p_15']

# branch_points = {
#   'b_1': ['p_1', 'p_4', 'p_15'],
#   'b_2': ['p_4', 'p_11', 'p_5'],
#   'b_3': ['p_8', 'p_9', 'p_10'],
#   'b_4': ['p_8', 'p_7', 'p_3'],
#   'b_5': ['p_13', 'p_6', 'p_14'],
#   'b_6': ['p_13', 'p_12', 'p_2'],}

# joints = {
#     'j_1': ['p_1', None, 'p_3', None, 'p_2', None],
#     'j_2': ['p_10', None, 'p_11', None, 'p_12', None],
#     'j_3': ['p_9', None, 'p_15', None, 'p_14', None],
#     'j_4': ['p_5', None, 'p_6', None, 'p_7', None],
# }

# homology_classes = {
#   'gamma_1' : ['p_1', 'p_2', 'p_3'],
#   'gamma_2' : ['p_4'],
#   'gamma_3' : ['p_9', 'p_14', 'p_15'],
#   'gamma_4' : ['p_10', 'p_11', 'p_12'],
#   'gamma_5' : ['p_5', 'p_6', 'p_7'],
#   'gamma_6' : ['p_13'],
#   'gamma_7' : ['p_8'],
# }

# bp_coordinates = {
#   'b_1': [0.27, 0.23],
#   'b_2': [0.51, 0.45],
#   'b_3': [0.49, 0.0],
#   'b_4': [0.78, 0.0],
#   'b_5': [0.0, 0.61],
#   'b_6': [0.0, 0.7],
# }

# j_coordinates = {
#     'j_1': [0.0, 0.0],
#     'j_2': [0.24, 0.0],
#     'j_3': [0.0, 0.26],
#     'j_4': [0.73, 0.76],
# }

# e_coordinates = {
#     'p_1': [j_coordinates['j_1'], bp_coordinates['b_1']],
#     'p_2': [bp_coordinates['b_6'], [j_coordinates['j_1'][0], j_coordinates['j_1'][1] + 1]], 
#     'p_3': [bp_coordinates['b_4'], [j_coordinates['j_1'][0] + 1, j_coordinates['j_1'][1]]], 
#     'p_4': [bp_coordinates['b_1'], bp_coordinates['b_2']], 
#     'p_5': [j_coordinates['j_4'], bp_coordinates['b_2']],
#     'p_6': [j_coordinates['j_4'], [bp_coordinates['b_5'][0] + 1, bp_coordinates['b_5'][1]]],
#     'p_7': [j_coordinates['j_4'], [bp_coordinates['b_4'][0], bp_coordinates['b_4'][1] + 1]],
#     'p_8': [bp_coordinates['b_3'], bp_coordinates['b_4']],
#     'p_9': [bp_coordinates['b_3'], [j_coordinates['j_3'][0] + 1, j_coordinates['j_3'][1]]],
#     'p_10': [bp_coordinates['b_3'], j_coordinates['j_2']],
#     'p_11': [bp_coordinates['b_2'], j_coordinates['j_2']],
#     'p_12': [bp_coordinates['b_6'], [j_coordinates['j_2'][0], j_coordinates['j_2'][1] + 1]],
#     'p_13': [bp_coordinates['b_5'], bp_coordinates['b_6']],
#     'p_14': [bp_coordinates['b_5'], j_coordinates['j_3']],
#     'p_15': [bp_coordinates['b_1'], j_coordinates['j_3']],
# }

# one = ['p_1', 'p_4', 'p_11', 'p_10', 'p_8', 'p_3']
# tau = ['p_1', 'p_15', 'p_14', 'p_13', 'p_2']




# # --------------- [3,1]-punctured torus SELF-GLUING (NOT WORKING) -----------------

# streets_1 = ['p_1', 'p_2', 'p_3', 'p_4', 'p_5', 'p_6', 'p_7', 'p_8', 'p_9', 'p_10', 'p_11', 'p_12', 'p_13', 'p_14', 'p_15']

# branch_points_1 = {
#   'b_1': ['p_1', 'p_2', 'p_3'],
#   'b_2': ['p_8', 'p_14', 'p_15'],
#   'b_3': ['p_5', 'p_11', 'p_12'],
#   'b_4': ['p_2', 'p_13', 'p_10'],
#   'b_5': ['p_7', 'p_8', 'p_9'],
#   'b_6': ['p_4', 'p_5', 'p_6'],}

# joints_1 = {
#     'j_1': ['p_9', None, 'p_10', None, 'p_11', None],
#     'j_2': ['p_1', None, 'p_6', None, 'p_7', None],
#     'j_3': ['p_3', None, 'p_15', None, 'p_4', None],
#     'j_4': ['p_12', None, 'p_13', None, 'p_14', None],
# }

# homology_classes_1 = {
#   'gamma_1' : ['p_1', 'p_6', 'p_7'],
#   'gamma_2' : ['p_2'],
#   'gamma_3' : ['p_3', 'p_4', 'p_15'],
#   'gamma_4' : ['p_12', 'p_13', 'p_14'],
#   'gamma_5' : ['p_9', 'p_10', 'p_11'],
#   'gamma_6' : ['p_5'],
#   'gamma_7' : ['p_8'],
# }

# bp_coordinates_1 = {
#   'b_1': [0.27, 0.23],
#   'b_2': [0.51, 0.45],
#   'b_3': [0.49, 0.0],
#   'b_4': [0.78, 0.0],
#   'b_5': [0.0, 0.61],
#   'b_6': [0.0, 0.7],
# }

# j_coordinates_1 = {
#     'j_1': [0.0, 0.0],
#     'j_2': [0.24, 0.0],
#     'j_3': [0.0, 0.26],
#     'j_4': [0.73, 0.76],
# }

# ### NOTE: unreliable tau etc
# e_coordinates_1 = None

# one_1 = ['p_1', 'p_6', 'p_4', 'p_3']
# tau_1 = ['p_10', 'p_11', 'p_12', 'p_13']