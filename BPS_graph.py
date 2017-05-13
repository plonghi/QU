import sys, itertools, re, numpy, os, datetime, matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mcsn import MCSN, Street, Joint, BranchPoint
from decimal import Decimal

from config import MCSNConfig

# an auxiliary parameter for checking whether, modulo [1, 1]
# two vectors in the plane are the same
EPSILON = 0.0001


class BPSgraph(MCSN):
    """
    A class for BPS graphs.
    The attribute 'mutable_edges' collects (labels of) edges
    which can be flipped.
    The attribute 'mutable_faces' collects (labels of) faces
    on which a cootie move can be performed.

    'e_coordinates' is a dictionary related to edges; 
    it contains coordinates of edges

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
        bp_coordinates={},
        j_coordinates={},
        e_coordinates={},
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
        self.e_coordinates = e_coordinates
        self.faces = {}
        self.mutable_faces = []
        self.mutable_edges = []
        self.check_edge_node_coordinates()
        self.determine_faces()
        for f in self.faces.values():
            f.determine_neighbors(self)
        self.tessellation_data = {
            f.label : f.neighbors for f in self.faces.values()
        }
        self.determine_mutable_elements()
    

    def determine_faces(self):
        # an array of all initial points
        # for building face paths (these are instances of 
        # StreetEndPoint)
        initial_points = []
        for s in self.streets.values():
            initial_points += s.endpoints

        # record all face paths while they are build to avoid
        # building same face twice from different initial points
        all_face_paths = []

        # build all face paths
        i_f = 0 # a counter
        for ep in initial_points:
            if ep in all_face_paths:
                continue
            else:
                f_label = 'f_'+str(i_f)
                new_face = Face(ep, label=f_label)
                self.faces[f_label] = new_face
                all_face_paths += new_face.full_sequence
                i_f += 1

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
            if can_be_flipped(s_v) is True:
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
                raise Exception(
                    'Coordinates of {} seems not to match ' 
                    'those of its endpoints. The edge has {}'
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
        return BPSgraph(
            branch_points=new_branch_points, 
            streets=new_streets, 
            joints=new_joints, 
            homology_classes=new_homology_classes,
            bp_coordinates=new_bp_coordinates,
            j_coordinates=new_j_coordinates,
            e_coordinates=new_e_coordinates,
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
    d = numpy.linalg.norm(numpy.array(xy_A_mod) - numpy.array(xy_B_mod))
    d_mod_1 = min(
        [abs(shift_number_into_interval_01(d)), 
        abs(shift_number_into_interval_01(-d))]
    )
    if d_mod_1 < eps:
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
    at the edge whose label is given as an argument
    """
    if edge not in graph.streets.keys():
        raise Exception('This edge is not part of the graph.')
    else:
        # the street
        p = graph.streets[edge]

    ep0, ep1 = p.endpoints

    # check that we can perform a flip on this edge: must be an I-web
    if can_be_flipped(p) is False:
        raise Exception('Cannot flip on edge {}.'.format(edge))

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

    b_0 = ep0.end_point.label
    slot_0 = ep0.slot
    old_node_xy_0 = graph.bp_coordinates[b_0]
    b_1 = ep1.end_point.label
    slot_1 = ep1.slot
    old_node_xy_1 = graph.bp_coordinates[b_1]

    del new_branch_points[b_0]
    del new_branch_points[b_1]

    # now I label the streets meeting at the two endpoints as follows:
    # at endpoint ep0 I have in ccw order: {p, s00, s01}
    # at endpoint ep1 I have in ccw order: {p, s10, s11}
    # then I must make new branch points where I have, in ccw order
    # {p, s11, s00} at b_0 and {p, s01, s10} at b_1
    # let's get the labels fo these streets then:
    s00 = ep0.end_point.streets[(slot_0 + 1) % 3].label
    s01 = ep0.end_point.streets[(slot_0 + 2) % 3].label
    s10 = ep1.end_point.streets[(slot_1 + 1) % 3].label
    s11 = ep1.end_point.streets[(slot_1 + 2) % 3].label

    # now add the new keys to the dictionary of branch points 
    new_branch_points[b_0] = [p.label, s11, s00]
    new_branch_points[b_1] = [p.label, s01, s10]

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
    # correspond to the coordinate of b_0 (resp 1), before flipping.
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

    # Also update positions of the branch points, now using
    # the separation induced by the edge coordinates, though!
    # TODO: get rid of this step, by automatically computing positions
    # of branch points and joints from the edge coordinates.
    # TODO: DECIDE whether keep all nodes within the unit square, 
    # or let them go around. Probably it's ok to have either choice.
    new_node_xy_0 = new_edge_xy_0
    new_node_xy_1 = new_edge_xy_1

    new_bp_coordinates = graph.bp_coordinates.copy()
    new_bp_coordinates[b_0] = new_node_xy_0
    new_bp_coordinates[b_1] = new_node_xy_1

    new_j_coordinates = graph.j_coordinates.copy()

    # finally return the new BPS graph
    return BPSgraph(
        branch_points=new_branch_points, 
        streets=new_streets, 
        joints=new_joints, 
        homology_classes=new_homology_classes,
        bp_coordinates=new_bp_coordinates,
        j_coordinates=new_j_coordinates,
        e_coordinates=new_e_coordinates,
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

    print (
        '\nWARNING: '
        'Homology classes not handled correctly in cootie.\n'
    )

    return BPSgraph(
        branch_points=new_branch_points, 
        streets=new_streets, 
        joints=new_joints, 
        homology_classes=new_homology_classes,
        bp_coordinates=new_bp_coordinates,
        j_coordinates=new_j_coordinates,
        e_coordinates=new_e_coordinates,
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

    # start by checking equality of branch points
    bp_1_perm = replace(s_1, p, bp_1)
    if equivalent_as_dictionaries(bp_1_perm, bp_2) is True:
        # then also check equality of joints
        j_1_perm = replace(s_1, p, j_1)
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
    """
    if len(old_vars) != len(new_vars):
        raise Exception('Replacement impossible')

    # First, let's turn the dictionary into a string
    dic_str = str(dic)

    # then, let's perform the subs, first create a replacement dictionary
    rep = {old_vars[i] : new_vars[i] for i in range(len(old_vars))}

    # taken from 
    # http://stackoverflow.com/questions/6116978/python-replace-multiple-strings
    rep = dict((re.escape(k), v) for k, v in rep.iteritems())
    pattern = re.compile("|".join(rep.keys()))
    new_dic_str = pattern.sub(
        lambda m: rep[re.escape(m.group(0))], dic_str
    )

    # finally, turn the string back into a dictionary
    return eval(new_dic_str)


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


def find_invariant_sequences(
    graph, depth, level=None, ref_graph=None, sequence=None,
    edge_cycle_min_length=None, min_n_cooties=None,
):
    """
    Find all sequences of flips or cootie moves, up to length 
    'depth' which take a BPS graph back to itself, up to an automorphism.
    The criterion for the automorphism is that the two graphs must 
    coincide up to a permutation of the streets.
    We try to match streets on both graphs by following the node structure.
    Only retain sequences in which streets are permuted with cycles 
    of a given minimal length.

    'level' is an auxiliary variable, it is used to keep track of recursion.
    'ref_graph' is the original graph, a sequence must eventually reproduce it.
    'sequence' is also an auxiliary variable, it keeps track of how many moves
    have been done, in fact level=len(sequence)
    'edge_cycle_min_length' determines whether to retain a valid sequence, 
    based on how the edges of the graph are permuted: if it's a cyclic 
    permutation of length greater than a certain minimal value
    'min_n_cooties' determines whether to retain a valid sequence based on the
    number of cootie moves it contains

    """
    if edge_cycle_min_length is None:
        edge_cycle_min_length = 0

    if min_n_cooties is None:
        min_n_cooties = 0

    if level is None:
        level = 0

    if level == depth:
        # print 'max depth reached'
        return []

    if sequence is None:
        sequence = []

    # print '\nStudying sequence {}'.format(sequence)
    
    if ref_graph is None:
        ref_graph = graph

    self_similar_graphs = []

    for f in graph.mutable_faces:
        # avoid repeating the last move in the sequence
        # that would be a trivial result (it's involutive)
        if len(sequence) > 0:
            if f == sequence[-1]:
                # print 'cootie on {} is the same as last move'.format(f)
                continue

        g_new = cootie_face(graph, f)
        new_sequence = sequence + [f]
        if have_same_face_types(ref_graph, g_new) is True:
            perms = are_equivalent_as_graphs(ref_graph, g_new)
            if perms is not None:
                # get ALL valid permutation dictionaries of edges
                
                # Now among all the permutations which would 
                # make two graphs coincide, we take the minimal one
                # only. Because we want to discard two graphs which 
                # can be related to each other by 'small permutations'
                min_perm = []
                n_edges = len(graph.streets.keys())
                min_l = n_edges
                # find the minimal length of all the maximal 
                # cycles within each permutation
                for p in perms:
                    # recall each permutation is a dictionary
                    p_k = p.keys()
                    p_v = [p[key] for key in p_k]
                    cycles = determine_perm_cycles(p_k, p_v)
                    if max(cycles) < min_l:
                        min_l = max(cycles)
                        min_perm = p
                # now let's see if there is a permutation with 
                # all cycles below the threshold cycle_min_length
                # in that case, we don't keep this graph.
                if min_l < edge_cycle_min_length:
                        continue

                else:
                    print (
                        'This is a good sequence: {}'.format(new_sequence)
                    )
                    self_similar_graphs.append([new_sequence, p])
        else:
            # print 'Will try going deeper with: {}'.format(new_sequence)
            deeper_sequences = find_invariant_sequences(
                g_new,
                depth,
                level=level+1,
                ref_graph=ref_graph,
                sequence=new_sequence,
                edge_cycle_min_length=edge_cycle_min_length,
                min_n_cooties=min_n_cooties,
            )
            self_similar_graphs += deeper_sequences

    for e in graph.mutable_edges:
        # avoid repeating the last move in the sequence
        # that would be a trivial result (it's involutive)
        # print 'studying mutation on {}'.format(e)
        if len(sequence) > 0:
            if e == sequence[-1]:
                # print 'flipping {} is the same as last move'.format(e)
                continue

        g_new = flip_edge(graph, e)
        new_sequence = sequence + [e]
        
        if have_same_face_types(ref_graph, g_new) is True:
            perms = are_equivalent_as_graphs(ref_graph, g_new)
            if perms is not None:
                # get ALL valid permutation dictionaries of edges
                
                # Now among all the permutations which would 
                # make two graphs coincide, we take the minimal one
                # only. Because we want to discard two graphs which 
                # can be related to each other by 'small permutations'
                min_perm = []
                n_edges = len(graph.streets.keys())
                min_l = n_edges
                # find the minimal length of all the maximal 
                # cycles within each permutation
                for p in perms:
                    # recall each permutation is a dictionary
                    p_k = p.keys()
                    p_v = [p[key] for key in p_k]
                    cycles = determine_perm_cycles(p_k, p_v)
                    if max(cycles) < min_l:
                        min_l = max(cycles)
                        min_perm = p
                # now let's see if there is a permutation with 
                # all cycles below the threshold cycle_min_length
                # in that case, we don't keep this graph.
                if min_l < edge_cycle_min_length:
                        continue

                else:
                    
                    print (
                        'This is a good sequence: {}'.format(new_sequence)
                    )
                    self_similar_graphs.append([new_sequence, p])
        else:
            # print 'Will try going deeper with: {}'.format(new_sequence)
            deeper_sequences = find_invariant_sequences(
                g_new,
                depth,
                level=level+1,
                ref_graph=ref_graph,
                sequence=new_sequence,
                edge_cycle_min_length=edge_cycle_min_length,
                min_n_cooties=min_n_cooties,
            )
            self_similar_graphs += deeper_sequences

    # determine which sequences to retain, based on several factors
    faces = graph.faces.keys()
    # edges = graph.streets.keys()

    # first drop empty ones
    selected_ssg = [ssg for ssg in self_similar_graphs if len(ssg) > 0]
    # require a certain number of cootie moves
    selected_ssg = [
        ssg for ssg in selected_ssg if 
        (   
            move_count(ssg[0], faces) >= min_n_cooties
        )
    ]

    return selected_ssg

def apply_sequence(graph, seq, save_plot=None):
    """
    Perform a sequence of moves.
    Gives an option to save the plots of all steps, 
    takes the directory path as argument for this.
    """
    if save_plot is not None:
        plot = plot_BPS_graph(graph)
        plt.savefig(os.path.join(save_plot, '0.png'))
        plt.clf()  # Clear the figure for the next loop

    faces = [k for k in graph.faces.keys()]
    edges = [k for k in graph.streets.keys()]

    # Note: this takes precisely the old graph, it doesn't create a new one
    # for now. This is OK. New graphs will be created at each flip and cootie.
    new_graph = graph 
    for i, el in enumerate(seq):
        if el in faces:
            new_graph = cootie_face(new_graph, el)
        elif el in edges:
            new_graph = flip_edge(new_graph, el)
        else:
            print 'got this command: {}'.format(el)
            raise Exception('Unknown step.')
        if save_plot is not None:
            plot = plot_BPS_graph(new_graph)
            plt.savefig(os.path.join(save_plot, str(i + 1)+'.png'))
            plt.clf()  # Clear the figure for the next loop

    return new_graph


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
    on each graph.
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
    # add new entries to the dictionary
    if have_compatible_endpoints(
        next_edge_1_L, next_edge_2_L
    ):
        edge_dic[next_edge_1_L.label] = next_edge_2_L.label
    else:
        # interrupt this because the two streets cannot be 
        # compatible with each other if their endpoints are different
        return None

    if have_compatible_endpoints(
        next_edge_1_R, next_edge_2_R
    ):
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

def match_graphs_by_edges(g1, g2, all_perms=False):
    """
    Does what the name says.
    Try to match graph g1 to g2 by taking the first edge of g1 
    with fixed orientation and trying to match onto every edge of g2 
    with both orientations, then follow through nodes of each 
    graph in parallel, stopping if a contradiction is found.

    An option to return ALL permutations whch work is available.
    """
    if all_perms is True:
        ok_perms = []

    for start_2 in g2.streets.keys():
        for orientation_2 in [+1, -1]:
            edge_dic = match_graphs_from_starting_point(
                g1, g2, start_2, orientation_2
            )
            if edge_dic is not None:
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
    Returns all permutations which make two graphs identical.
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
                print 'OK found a dictionary for edges'
                ok_perms.append(perm)
                pass
            else:
                print 'Discard permutation.'
                pass
        if len(ok_perms) > 0:
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


# --------------- [2,1]-punctured torus -----------------

streets = ['p_1', 'p_2', 'p_3', 'p_4', 'p_5', 'p_6', 'p_7', 'p_8', 'p_9']

branch_points = {
  'b_1': ['p_1', 'p_2', 'p_3'],
  'b_2': ['p_1', 'p_4', 'p_8'],
  'b_3': ['p_2', 'p_5', 'p_9'],
  'b_4': ['p_3', 'p_6', 'p_7'],}

joints = {
    'j_1': ['p_4', None, 'p_5', None, 'p_6', None],
    'j_2': ['p_7', None, 'p_8', None, 'p_9', None],
}

homology_classes = {
  'gamma_1' : ['p_1'],
  'gamma_2' : ['p_2'],
  'gamma_3' : ['p_3'],
  'gamma_4' : ['p_4', 'p_5', 'p_6'],
  'gamma_5' : ['p_7', 'p_8', 'p_9'],
}

bp_coordinates = {
  'b_1': [0.0, 0.0],
  'b_2': [0.28, 0.0],
  'b_3': [0.0, 0.31],
  'b_4': [0.62, 0.57],
}

j_coordinates = {
    'j_1': [0.32, 0.63],
    'j_2': [0.62, 0.28],
}

e_coordinates = {
    'p_1': [bp_coordinates['b_1'], bp_coordinates['b_2']],
    'p_2': [bp_coordinates['b_1'], bp_coordinates['b_3']], 
    'p_3': [bp_coordinates['b_4'], [bp_coordinates['b_1'][0] + 1, bp_coordinates['b_1'][1] + 1]], 
    'p_4': [j_coordinates['j_1'], [bp_coordinates['b_2'][0], bp_coordinates['b_2'][1] + 1] ], 
    'p_5': [j_coordinates['j_1'], bp_coordinates['b_3']],
    'p_6': [j_coordinates['j_1'], bp_coordinates['b_4']],
    'p_7': [j_coordinates['j_2'], bp_coordinates['b_4']],
    'p_8': [j_coordinates['j_2'], bp_coordinates['b_2']],
    'p_9': [j_coordinates['j_2'], [bp_coordinates['b_3'][0] + 1, bp_coordinates['b_3'][1]]],
}




# # --------------- [3,1]-punctured torus -----------------

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



w = BPSgraph(
    branch_points=branch_points, 
    streets=streets, 
    joints=joints, 
    homology_classes=homology_classes,
    bp_coordinates=bp_coordinates,
    j_coordinates=j_coordinates,
    e_coordinates=e_coordinates,
)

# analyze all sequences which give back the BPS graph
max_n_moves = 5
seq = find_invariant_sequences(
    w, max_n_moves, level=0, ref_graph=w, 
    edge_cycle_min_length=2, 
    min_n_cooties=1,
)
# print [s[0] for s in seq]
print 'Found {} sequences.'.format(len(seq))

# prepare a directory
mydir = os.path.join(
    os.getcwd(), 'mcg_moves', 
    (datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S') + '_{}_moves'.format(max_n_moves))
)
os.makedirs(mydir)

# save info about sequences
text_file = open(mydir + '/sequence_data.txt', 'w')
text_file.write('\t\tSequence data\n\n')
for i, s in enumerate(seq):
    text_file.write('\nsequence {} : {}'.format(i, s[0]))

# save plots of sequences
for i, s in enumerate(seq):
    sub_dir = os.path.join(mydir, 'sequence_'+str(i))
    os.makedirs(sub_dir)
    apply_sequence(w, s[0], save_plot=sub_dir)

