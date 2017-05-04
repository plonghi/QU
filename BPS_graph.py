import sys, itertools, re

from mcsn import MCSN, Street, Joint, BranchPoint
from solitons import (
    copy_of_soliton, join_dashes_at_branch_point, SolitonPath, Dash, 
    growing_pairs_from_dashes, growing_clusters
)
from growth_rules import (
    # bp_type_1_growth, j_type_3_growth, 
    grow_soliton_once, grow_soliton,
)
from soliton_data import SolitonData, NetworkSolitonContent

from intersections import get_dash_nodes, compute_self_intersections

from config import MCSNConfig



class BPSgraph(MCSN):
    """
    A class for BPS graphs.
    The attribute 'mutable_edges' collects (labels of) edges
    which can be flipped.
    The attribute 'mutable_faces' collects (labels of) faces
    on which a cootie move can be performed.
    """
    def __init__(
        self, 
        label='no_label', 
        branch_points={}, 
        streets=[], 
        joints={}, 
        homology_classes={}
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
        self.faces = {}
        self.mutable_faces = []
        self.mutable_edges = []
        self.determine_faces()
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
                new_face = Face(ep)
                self.faces['f_'+str(i_f)] = new_face
                all_face_paths += new_face.full_sequence
                i_f += 1

    def print_face_info(self):
        print '\nThere are {} faces, they are:'.format(len(self.faces))
        for i, f in enumerate(self.faces.values()):
            print '\nFace {}.'.format(i)
            f.print_info()

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


class Face():
    """
    A class for faces of BPS graphs
    """
    def __init__(self, endpoint):
        self.full_sequence = []
        self.node_sequence = []
        self.street_sequence = []
        [
            self.full_sequence, 
            self.node_sequence, 
            self.street_sequence
        ] = build_face(endpoint)

        # the face type is a list of letters like ['b', 'b', 'j',...]
        # which records the types of nodes at the feca boundary
        # (b is for branch point, j for joint)
        self.face_type = []
        for el in self.full_sequence:
            if el.__class__.__name__ == 'StreetEndpoint':
                if el.end_point.__class__.__name__ == 'BranchPoint':
                    self.face_type += ['b']
                elif el.end_point.__class__.__name__ == 'Joint':
                    self.face_type += ['j']
                else:
                    pass

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

    new_streets = graph.streets.keys()
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
    b_1 = ep1.end_point.label
    slot_1 = ep1.slot

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

    return BPSgraph(
        branch_points=new_branch_points, 
        streets=new_streets, 
        joints=new_joints, 
        homology_classes=new_homology_classes
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

    new_streets = graph.streets.keys()
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

    print (
        '\n\t **** WARNING **** \n'
        'Homology classes not handled correctly in cootie for the moment.\n'
    )

    return BPSgraph(
        branch_points=new_branch_points, 
        streets=new_streets, 
        joints=new_joints, 
        homology_classes=new_homology_classes
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

def are_equivalent(graph_1, graph_2):
    """
    Determine whether two graphs are equivalent, in the 
    sense that there exists a permutation of the edges
    which makes them identical.
    If they are inequivalent, returns None.
    If they are equivalent, returns the permutation.
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

    # Now apply a permutation to the set of streets s_1
    # and recast branch points bp_1 and joints j_1
    # in the new labels
    # If they both coincide with bp_2 and j_2 respectively,
    # the two graphs are equivalent!

    # create all permutations of the list of streets
    all_perms = map(list, list(itertools.permutations(
        s_1
    )))

    for p in all_perms:
        # start by checking equality of branch points
        bp_1_perm = replace(s_1, p, bp_1)
        if equivalent_as_dictionaries(bp_1_perm, bp_2) is True:
            # then also check equality of joints
            j_1_perm = replace(s_1, p, j_1)
            if equivalent_as_dictionaries(j_1_perm, j_2) is True:
                return [s_1, p]

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


def find_invariant_sequences(
    graph, depth, level=0, ref_graph=None, sequence=None,
):
    """
    Find all sequences of flips or cootie moves, up to length 
    'depth' which take a BPS graph back to itself, up to an automorphism.
    The criterion for the automorphism is that the two graphs must 
    coincide up to a permutation of the streets.
    Since the permutation of streets can take much time, 
    to speed up the comparison we first make weaker checks.
    For instance we check that two graphs have the same types of faces.
    """

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
            perm = are_equivalent(ref_graph, g_new)
            if perm is not None:
                print 'This is a good sequence: {}'.format(new_sequence)
                self_similar_graphs.append([new_sequence, perm])
        else:
            # print 'Will try going deeper with: {}'.format(new_sequence)
            deeper_sequences = find_invariant_sequences(
                g_new,
                depth,
                level=level+1,
                ref_graph=ref_graph,
                sequence=new_sequence
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
            # print 'HAVE SAME FACES with sequence {}'.format(new_sequence)
            perm = are_equivalent(ref_graph, g_new)
            if perm is not None:
                print 'This is a good sequence: {}'.format(new_sequence)
                self_similar_graphs.append([new_sequence, perm])
        else:
            # print 'Will try going deeper with: {}'.format(new_sequence)
            deeper_sequences = find_invariant_sequences(
                g_new,
                depth,
                level=level+1,
                ref_graph=ref_graph,
                sequence=new_sequence
            )
            self_similar_graphs += deeper_sequences

    return [ssg for ssg in self_similar_graphs if len(ssg) > 0]


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



w = BPSgraph(
    branch_points=branch_points, 
    streets=streets, 
    joints=joints, 
    homology_classes=homology_classes
)

# # now reshuffle the streets of the network

# old_vars = streets
# perms = map(list, list(itertools.permutations(
#         old_vars
#     )))
# new_vars = perms[10]
# print 'old vars = {}'.format(old_vars)
# print 'new vars = {}'.format(new_vars)
# bp1 = replace(old_vars, new_vars, branch_points)
# j1 = replace(old_vars, new_vars, joints)

# w1 = BPSgraph(
#     branch_points=bp1, 
#     streets=streets, 
#     joints=j1, 
#     homology_classes=homology_classes
# )

# print 'the two graphs have same faces: {}'.format(have_same_face_types(w, w1))
# print 'the two graphs are equivalent: {}'.format(are_equivalent(w, w1))


# w.print_face_info()

# w.print_mutable_info()

# w1 = flip_edge(w, 'p_6')
# w2 = flip_edge(w, 'p_4')
# w3 = flip_edge(w, 'p_5')
# print have_same_face_types(w, w3)

# w2 = cootie_face(w, 'f_1')

# w1.print_face_info()

# print have_same_face_types(w, w)
# print have_same_face_types(w, w1)
# print have_same_face_types(w, w2)
# print have_same_face_types(w1, w2)

seq = find_invariant_sequences(w, 5, level=0, ref_graph=w,)
print seq

# w.print_face_info()

# w_1 = flip_edge(w, 'p_2')
# w_2 = flip_edge(w_1, 'p_1')
# w_3 = cootie_face(w_2, 'f_1')
# w_4 = flip_edge(w_3, 'p_4')
# w_5 = flip_edge(w_4, 'p_2')
# print have_same_face_types(w, w_5)
# print are_equivalent(w, w_5)

# w_1 = flip_edge(w, 'p_1')
# w_2 = flip_edge(w_1, 'p_5')

# print have_same_face_types(w, w_2)
# print are_equivalent(w, w_2)
# # w_2.print_face_info()

# # now reshuffle the streets of the network



# graph_1 = w
# graph_2 = w_2

# s_1 = graph_1.streets.keys()
# bp_1 = {
#     b.label : [s.label for s in b.streets] 
#     for b in graph_1.branch_points.values()
# }
# j_1 = {
#     j.label : [get_label(s) for s in j.streets] 
#     for j in graph_1.joints.values()
# }

# s_2 = graph_2.streets.keys()
# bp_2 = {
#     b.label : [s.label for s in b.streets] 
#     for b in graph_2.branch_points.values()
# }
# j_2 = {
#     j.label : [get_label(s) for s in j.streets] 
#     for j in graph_2.joints.values()
# }

# old_vars = ['p_1', 'p_2', 'p_3', 'p_4', 'p_5', 'p_6']
# new_vars = ['p_1', 'p_4', 'p_2', 'p_6', 'p_5', 'p_3']
# print 'old vars = {}'.format(old_vars)
# print 'new vars = {}'.format(new_vars)
# new_bp = replace(old_vars, new_vars, bp_1)

# print new_bp
# print bp_2
# print equivalent_as_dictionaries(new_bp, bp_2)


# print 'the two graphs have same faces: {}'.format(have_same_face_types(w, w1))
# print 'the two graphs are equivalent: {}'.format(are_equivalent(w, w1))



# print '\n\n-------------------------------------------------------'
# print '\nSoliton Data'
# Q1 = SolitonData(label='Q_1', network=w , street=s3, resolution='british')
# Q1.initialize()
# # # Q1.print_info(full_path=True)
# Q1.grow(n_steps=4)
# Q1.print_info(full_path=True, soliton_paths=True)


