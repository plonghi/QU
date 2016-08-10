from sympy import symbols

TYPE_3_JOINTS = (
    # 3-way
    [
        ['Street', 'NoneType', 'Street', 'NoneType', 'Street', 'NoneType'],
        ['NoneType', 'Street', 'NoneType', 'Street', 'NoneType', 'Street'],
    ]
)

TYPE_4_A_JOINTS = (
    # 4-way - peace sign
    [
        ['NoneType', 'Street', 'NoneType', 'Street', 'Street', 'Street'],
        ['Street', 'NoneType', 'Street', 'Street', 'Street', 'NoneType'],
        ['NoneType', 'Street', 'Street', 'Street', 'NoneType', 'Street'],
        ['Street', 'Street', 'Street', 'NoneType', 'Street', 'NoneType'],
        ['Street', 'Street', 'NoneType', 'Street', 'NoneType', 'Street'],
        ['Street', 'NoneType', 'Street', 'NoneType', 'Street', 'Street'],
    ]
)

TYPE_4_B_JOINTS = (
    # 4-way - X sign
    [
        ['Street', 'Street', 'NoneType', 'Street', 'Street', 'NoneType'],
        ['Street', 'NoneType', 'Street', 'Street', 'NoneType', 'Street'],
        ['NoneType', 'Street', 'Street', 'NoneType', 'Street', 'Street'],
    ]
)

TYPE_5_JOINTS = (
    # 5-way
    [
        ['NoneType', 'Street', 'Street', 'Street', 'Street', 'Street'],
        ['Street', 'NoneType', 'Street', 'Street', 'Street', 'Street'],
        ['Street', 'Street', 'NoneType', 'Street', 'Street', 'Street'],
        ['Street', 'Street', 'Street', 'NoneType', 'Street', 'Street'],
        ['Street', 'Street', 'Street', 'Street', 'NoneType', 'Street'],
        ['Street', 'Street', 'Street', 'Street', 'Street', 'NoneType'],
    ]
)

TYPE_6_JOINTS = (
    # 6-way
    [
        ['Street', 'Street', 'Street', 'Street', 'Street', 'Street'],
    ]
)

ALLOWED_JOINTS = (
    TYPE_3_JOINTS +
    TYPE_4_A_JOINTS +
    TYPE_4_B_JOINTS +
    TYPE_5_JOINTS +
    TYPE_6_JOINTS
)

TYPE_1_BRANCH_POINTS = (
    [
        ['Street', 'NoneType', 'NoneType'],
        ['NoneType', 'Street', 'NoneType'],
        ['NoneType', 'NoneType', 'Street']
    ]
)

TYPE_2_BRANCH_POINTS = (
    [
        ['Street', 'Street', 'NoneType'],
        ['NoneType', 'Street', 'Street'],
        ['Street', 'NoneType', 'Street']
    ]
)

TYPE_3_BRANCH_POINTS = (
    [
        ['Street', 'Street', 'Street'],
    ]
)

ALLOWED_BRANCH_POINTS = (
    TYPE_1_BRANCH_POINTS + TYPE_2_BRANCH_POINTS + TYPE_3_BRANCH_POINTS
)


class BranchPoint:
    """
    The Branch Point class.

    All streets ending on a branch point must be given 
    when the object is instanced.
    The ordering is important, it is understood to be 
    counter-clockwise, as in 1204.4824 figure 47 on page 74.
    """

    def __init__(self, label=None, streets=None):
        self.label = label
        if streets is None:
            self.streets = [None for i in range(3)]
        elif len(streets) < 3:
            self.streets = streets + [None for i in range(3 - len(streets))]
        elif len(streets) == 3:
            self.streets = streets
        else:
            raise Exception(
                'Cannot handle branch points with more than three streets.'
            )
        self.type = self.determine_type()
        self.available_slots = self.determine_available_slots()

    def street_position(self, street):
        """
        Note, there may be cases in which a street ends twice on a 
        branch point or joint.
        """
        positions = []
        for i, s in enumerate(self.streets):
            if street == s:
                positions.append(i)
        return positions
        
        # try:
        #   ind = self.streets.index(street)
        # except ValueError:
        #   ind = None
        # return ind

    def determine_type(self):
        bp_type = [st.__class__.__name__ for st in self.streets]
        if bp_type in TYPE_1_BRANCH_POINTS:
            return 'type_1_branch_point'
        elif bp_type in TYPE_2_BRANCH_POINTS:
            return 'type_2_branch_point'
        elif bp_type in TYPE_3_BRANCH_POINTS:
            return 'type_3_branch_point'
        else:
            return None

    def print_type(self):
        print '{} is a {}'.format(self.label, self.type)

    def determine_available_slots(self):
        bp_type = [st.__class__.__name__ for st in self.streets]
        return [i for i, x in enumerate(bp_type) if x == 'Street']

    def print_info(self):
        self.print_type()
        slot_labels = []
        for s in self.streets:
            if s is None:
                slot_labels.append('empty')
            else:
                slot_labels.append(s.label)
        print 'Streets attached (c.c.w.): {}\n'.format(slot_labels)


class Street:
    def __init__(self, label=None):
        self.label = label
        self.endpoints = []
        pass

    def add_endpoint(self, end_point):
        if len(self.endpoints) == 2:
            raise Exception(
                'Street {} cannot have more than two endpoints.'
                .format(self.label)
            )
        else:
            self.endpoints.append(end_point)

    def initial_point(self):
        return self.endpoints[0]
        
    def final_point(self):
        return self.endpoints[1]

    def print_info(self):
        print 'Street {} is oriented {} ---> {}\n'.format(
            self.label, 
            self.initial_point().end_point.label, 
            self.final_point().end_point.label
        )


class StreetEndpoint:
    def __init__(self, street=None, end_point=None, slot=None):
        self.street = street
        self.end_point = end_point
        self.slot = slot


class Joint:
    """
    All streets ending on a joint must be
    given when the object is instanced.
    The ordering is important, it is understood to be 
    counter-clockwise, as in 1204.4824 figure 46 on page 73.
    Very important: All 6 streets must be specified! If a street 
    is one-way, the corresponding entry will be the keyword "None".
    """
    def __init__(self, label=None, streets=None):
        self.label = label
        if streets is None:
            self.streets = [None for i in range(6)]
        elif len(streets) == 6:
            self.streets = streets
        else:
            raise Exception(
                'Must specify 6 streets for each joint.'
            )
        self.type = self.determine_type()
        self.available_slots = self.determine_available_slots()

    def street_position(self, street):
        """
        Note, there may be cases in which a street ends twice on a 
        branch point or joint.
        """
        positions = []
        for i, s in enumerate(self.streets):
            if street == s:
                positions.append(i)

        return positions

        # try:
        #   ind = self.streets.index(street)
        # except ValueError:
        #   ind = None
        # return ind

    def determine_type(self):
        j_type = [st.__class__.__name__ for st in self.streets]
        if j_type in TYPE_3_JOINTS:
            return 'type_3_joint'
        elif j_type in TYPE_4_A_JOINTS:
            return 'type_4_A_joint'
        elif j_type in TYPE_4_B_JOINTS:
            return 'type_4_B_joint'
        elif j_type in TYPE_5_JOINTS:
            return 'type_5_joint'
        elif j_type in TYPE_6_JOINTS:
            return 'type_6_joint'
        else:
            return None

    def print_type(self):
        print '{} is a {}'.format(self.label, self.type)

    def determine_available_slots(self):
        j_type = [st.__class__.__name__ for st in self.streets]
        return [i for i, x in enumerate(j_type) if x == 'Street']

    def print_info(self):
        self.print_type()
        slot_labels = []
        for s in self.streets:
            if s is None:
                slot_labels.append('empty')
            else:
                slot_labels.append(s.label)
        print 'Streets attached (c.c.w.): {}\n'.format(slot_labels)


# Important: must develop this class, 
# to extract the charges to assign to each
# lift of solitons.
class StringUnit:
    def __init__(self, label=None):
        self.label = label
        pass


class MCSN:
    """
    The network is constructed by specifying 
    the streets attached to each joint and branch point.
    Then, the method 'attach_streets()' is called.

    Arguments: 
    
    streets = ['p_1', 'p_2', ...] 
        (a list of labels for network streets)
    
    branch_points = {
        'b_1': ['p_1', 'p_2', 'p_3'],
        'b_2': ['p_1', 'p_3'],
        'b_3': ['p_2'],
        ...
    }
        (a dictionary of lists of labels, streets attached to a branch point 
        are ordered counter-clockwise)
    
    joints = {
        'j_1': ['p_1', 'p_2', None, 'p_3', None, None],
        ...
    }
        (a dictionary of lists of labels, streets attached to a joint 
        are ordered counter-clockwise)
    
    homology_classes = {
        'gamma_1' : ['p_1'],
        'gamma_2' : ['p_2', 'p_3'],
        ...
    }
        (a dictionary of lists of labels, each homology class
        is characterized as the sum of lifts of a set of streets.
        The ordering doesn't matter)

    """
    def __init__(
        self, 
        label='no_label', 
        branch_points={}, 
        streets=[], 
        joints={}, 
        homology_classes={}
    ):
        self.label = label
        
        self.streets = (
            {s_label: Street(label=s_label) for s_label in streets}
        )
        
        self.branch_points = {}
        for bp_lbl in branch_points.keys():
            # list of streets at the branch point
            streets_list = []
            for str_lbl in branch_points[bp_lbl]:
                if str_lbl is None:
                    streets_list.append(None)
                else:
                    streets_list.append(self.streets[str_lbl])

            # add the branch point
            self.branch_points[bp_lbl] = BranchPoint(
                label=bp_lbl, streets=streets_list
            )

        self.joints = {}
        for j_lbl in joints.keys():
            # list of streets at the joint
            streets_list = []
            for str_lbl in joints[j_lbl]:
                if str_lbl is None:
                    streets_list.append(None)
                else:
                    streets_list.append(self.streets[str_lbl])

            # add the joint
            self.joints[j_lbl] = Joint(
                label=j_lbl, streets=streets_list
            )

        self.basis_homology_classes = (
            {hc_lbl: HomologyClass(
                atomic_labels=basis_atomic_labels(
                    i, homology_classes.keys()
                ), 
                streets=[
                    self.streets[str_lbl] 
                    for str_lbl in homology_classes[hc_lbl]
                ]
            ) for i, hc_lbl in enumerate(homology_classes.keys())}
        )
        
        # create the trivial homology class
        self.trivial_homology_class = HomologyClass(
            atomic_labels={
                lbl: 0 for lbl in self.basis_homology_classes.keys()
            },
            streets=[],
        )

        self.attach_streets()
        self.check_network()

    def attach_streets(self):
        for b_pt in self.branch_points.values():
            for i, street in enumerate(b_pt.streets):
                if street is not None:
                    street.add_endpoint(
                        StreetEndpoint(
                            street=street, end_point=b_pt, slot=i,
                        )
                    )
        
        for j_pt in self.joints.values():
            for i, street in enumerate(j_pt.streets):
                if street is not None:
                    street.add_endpoint(
                        StreetEndpoint(
                            street=street, end_point=j_pt, slot=i,
                        )
                    )

    def check_network(self):
        # Check that every street has exactly two endpoints
        for street in self.streets.values():
            if len(street.endpoints) == 2:
                pass
            else:
                raise Exception(
                    'Street {} does not have two endpoints, but {}'
                    .format(
                        street.label, [pt.end_point for pt in street.endpoints]
                    )
                )
        print 'All streets have two well-defined endpoints.'

        # Check that every joint is of an allowed type
        for j_pt in self.joints.values():
            j_type = j_pt.type
            is_allowed = False
            for a_joint in ALLOWED_JOINTS:
                if j_type is not None:
                    is_allowed = True
                    break
                else:
                    pass

            if is_allowed is False:
                raise Exception(
                    'Joint {} is of type {}'.format(j_pt.label, j_type)
                )
        print 'All joints are of a well-defined type.'

        # Check that homology classes include all streets of the network
        # and that no street is repeated more than once
        if len(self.basis_homology_classes) > 0:
            # a list with all streets in all homology classes
            all_hom_streets = []
            for hc in self.basis_homology_classes.values():
                all_hom_streets += hc.streets

            # check that every street of the network is contained in at 
            # least one homology class
            for s in self.streets.values():
                if s not in all_hom_streets:
                    raise Exception(
                        'Street {} is not contained in any homology class'
                        .format(s.label)
                    )

            # check that no street is contained in more than one homology class
            if not len(all_hom_streets) == len(self.streets.values()):
                print '\nWARNING: Two or more homology classes share a street.'
                # raise Exception(
                #     'Two or more homology classes share a street.'
                # )

            print 'All homology classes are well-defined.'
        else:
            print 'No homology classes have been defined.'

        # Print info about homology classes and the formal variables
        print '\nDictionary of formal variables and homology basis:'
        for hc_k in sorted(self.basis_homology_classes.keys()):
            hc_v = self.basis_homology_classes[hc_k]
            print '{} ---> {}'.format(hc_v.label, hc_v.symbol)
        
    def add_street(self):
        pass

    def add_joint(self):
        pass

    def add_branch_point(self):
        pass

    def print_info(self):
        print '\nSpectral Network Data\n-------------------------'

        print '\nSTREETS:\n'
        for s in self.streets.values():
            s.print_info()

        print '\nBRANCH POINTS:\n'
        for bp in self.branch_points.values():
            bp.print_info()

        print '\nJOINTS:\n'
        for j in self.joints.values():
            j.print_info()

        print '\nHOMOLOGY CLASSES:\n'
        for hc in self.basis_homology_classes.values():
            hc.print_info()


class HomologyClass:
    def __init__(self, atomic_labels={}, streets=[],):
        """
        Atomic labels is a dictionary based on a choice of basis for the
        homology lattice
        {'gamma_1': 1, 'gamma_2': 0, 'gamma_3': 2}
            (for gamma_1 + 2 * gamma_3)
        """

        self.atomic_labels = atomic_labels
        self.label = homology_label(self.atomic_labels)

        self.streets = sorted(streets)
        # TODO: get rid of this one below
        self.symbol_long = symbols('X_(' + self.label + ')')

        # create a list of variables 
        # [X1, X2, X3, ...]
        # defined as sympy symbols
        # they will correspond to the atomic 
        # labels for homology classes in the same ordering
        variables = [symbols('X'+str(i+1)) for i in range(len(atomic_labels))]
        sorted_keys = sorted(atomic_labels.keys())
        X_expression = 1
        for i, k_i in enumerate(sorted_keys):
            X_expression *= variables[i] ** atomic_labels[k_i]
        self.symbol = X_expression

    def __add__(self, other):
        new_atomic_labels = {
            a_l_k : self.atomic_labels[a_l_k] + other.atomic_labels[a_l_k]
            for a_l_k in self.atomic_labels.keys()
        }
        new_streets = self.streets + other.streets
        return HomologyClass(
            atomic_labels=new_atomic_labels, streets=new_streets
        )

    def print_info(self):
        print 'Homology class {} uses symbol {}.'.format(
            self.label, self.symbol
        )
        print 'Streets: {}\n'.format(
            [s.label for s in self.streets]
        )


def basis_atomic_labels(i, hc_lbls):
    """
    Returns the atomic labels for the i-th basis element of the 
    homology lattice e.g.
    {'gamma_1': 0, 'gamma_2': 1, 'gamma_3': 0, 'gamma_4': 0, }
    for the 2nd element (i.e. i=1)
    """
    value_list = [0 for j in range(len(hc_lbls))]
    value_list[i]=1
    return {lbl: value_list[j] for j, lbl in enumerate(hc_lbls)}


def homology_label(atomic_labels):
    """
    Atomic labels is a dictionary based on a choice of basis for the
    homology lattice
    {'gamma_1': 1, 'gamma_2': 0, 'gamma_3': 2}
        (for gamma_1 + 2 * gamma_3)
    """

    # initialize the label
    nonzero_labels = [
        i for i, a_l_k in enumerate(atomic_labels.keys()) 
        if atomic_labels[a_l_k]!=0
    ]
    if len(nonzero_labels) == 0:
        # the trivial homology class
        return ''
    else:
        first_nonzero_label_i = nonzero_labels[0]

        first_label = atomic_labels.keys()[first_nonzero_label_i]
        first_coeff = atomic_labels[first_label]
        
        if first_coeff == 1:
            label = first_label
        else:
            label = str(first_coeff) + '_*_' + first_label


        # add the remaining homology basis elements
        for a_l_k in atomic_labels.keys()[first_nonzero_label_i+1:]:
            if atomic_labels[a_l_k] == 0:
                pass
            elif atomic_labels[a_l_k] == 1:
                label += '_+_' + a_l_k
            elif atomic_labels[a_l_k] < 0:
                label += str(atomic_labels[a_l_k]) + '_*_' + a_l_k
            else:
                label += '_+_' + str(atomic_labels[a_l_k]) + '_*_' + a_l_k

        return label


# w = MCSN()
# # w.check_network()
# w.attach_streets()
# w.check_network()
# s1 = w.joints['j_1'].streets[0]
# print w.joints['j_1'].street_position(w.streets['p_2'])
# print w.branch_points['b_1'].street_position(w.streets['p_2'])
# print s1.initial_point().end_point
# print w.joints['j_1'].street_position(w.streets['p_1'])