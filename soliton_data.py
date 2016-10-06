"""
Module to handle soliton generating functions
"""

from solitons import SolitonPath, ClosedSoliton
from growth_rules import grow_soliton
from sympy import symbols

# Default maximal growth level for solitons.
N_MAX_GROWTH = 10


class SolitonData:
    """
    A container for soliton data supported on a street.
    For a two-way street of type ij/ji it will collect 
    all the solitons of each type.
    It will also provide a functionality for growing solitons,
    either by a number of steps, or until all of them are 'complete'.

    The soliton data of a street will moreover depend on a choice of
    resolution for the critical network, the resolution may be either 
    'american' or 'british'.
    """
    def __init__(self, label=None, network=None, street=None, resolution=None):
        self.label = label
        self.network = network
        self.street = street
        self.resolution = resolution

        if resolution is None:
            raise ValueError('Must specify the resolution of the network.')
        elif resolution=='american' or resolution=='british':
            pass
        else:
            raise ValueError('Unknown resolution: {}'.format(resolution))
        
        # given the orientation of the street, let that be the 
        # orientation of the underlying one-way ij-street. 
        # Then co-oriented solitons will be of ij type, while
        # anti-oriented solitons will be of ji type.
        self.co_oriented_solitons = []
        self.anti_oriented_solitons = []
        self.closed_solitons = []
        self.closed_solitons_reversed = []
        self.Q = 0
        self.Q_y = 0
        self.Q_y_reversed = 0

    def initialize(self):
        # The initial co-oriented soliton
        ij_sol = SolitonPath(
            label=self.label + '_ij_sol',resolution=self.resolution,
        )
        ij_source_pt = self.street.initial_point().end_point
        ij_source_slot = self.street.initial_point().slot
        ij_sol.create(
            street=self.street, 
            source_pt=ij_source_pt, 
            slot=ij_source_slot,
        )

        # The initial anti-oriented soliton
        ji_sol = SolitonPath(
            label=self.label + '_ji_sol', resolution=self.resolution,
        )
        ji_source_pt = self.street.final_point().end_point
        ji_source_slot = self.street.final_point().slot
        ji_sol.create(
            street=self.street, 
            source_pt=ji_source_pt, 
            slot=ji_source_slot,
        )

        self.co_oriented_solitons = [ij_sol]
        self.anti_oriented_solitons = [ji_sol]

    def grow(self, n_steps=N_MAX_GROWTH):
        # grow co-oriented solitons
        new_co_or_sols = []
        for sol in self.co_oriented_solitons:
            new_sols = grow_soliton(
                sol, n_steps=n_steps, resolution=self.resolution
            )
            new_co_or_sols += new_sols

        self.co_oriented_solitons = new_co_or_sols

        # grow anti-oriented solitons
        new_anti_or_sols = []
        for sol in self.anti_oriented_solitons:
            new_sols = grow_soliton(
                sol, n_steps=n_steps, resolution=self.resolution
            )
            new_anti_or_sols += new_sols

        self.anti_oriented_solitons = new_anti_or_sols

        # compute closed solitons
        self.compute_closed_solitons()

    def print_info(self, full_path=False, soliton_paths=False, writhes=False):
        print (
            '\n\nData of street {}'
            '\nIn {} resolution\n----------------------------'
            .format(self.street.label, self.resolution)
        )

        if soliton_paths is True:
            print (
                '\nCo-oriented solitons on {}:\n----------------------------'
                .format(self.street.label)
            )

            for i, s in enumerate(self.co_oriented_solitons):
                print '{}.'.format(i + 1)
                s.print_info(full_path=full_path)
            print (
                '\nAnti-oriented solitons on {}:\n------------------------------'
                .format(self.street.label)
            )
            for i, s in enumerate(self.anti_oriented_solitons):
                print '{}.'.format(i + 1)
                s.print_info(full_path=full_path)

            # print (
            #     '\nClosed solitons on {}:\n------------------------------'
            #     .format(self.street.label)
            # )

            # for i, s in enumerate(self.closed_solitons):
            #     print '{}.'.format(i + 1)
            #     s.print_info()

            # print (
            #     '\nClosed solitons on {} with opposite choice of basepoint:'
            #     '\n------------------------------'
            #     .format(self.street.label)
            # )

            # for i, s in enumerate(self.closed_solitons_reversed):
            #     print '{}.'.format(i + 1)
            #     s.print_info()

        if soliton_paths is True or writhes is True:
            print (
                '\nClosed solitons on {}:\n------------------------------'
                .format(self.street.label)
            )
            for i, s in enumerate(self.closed_solitons):
                print '{}.'.format(i + 1)
                s.print_info()

            print (
                '\nClosed solitons on {} with opposite choice of basepoint:'
                '\n------------------------------'
                .format(self.street.label)
            )
            for i, s in enumerate(self.closed_solitons_reversed):
                print '{}.'.format(i + 1)
                s.print_info()

        # print '\n4d Soliton generating function (without spin)'
        # print self.Q
        print '\n4d Soliton generating function (with spin)'
        print self.Q_y
        print '\n4d Soliton generating function with opposite basepoint (with spin)'
        print self.Q_y_reversed

    def compute_closed_solitons(self):
        self.closed_solitons = []
        self.closed_solitons_reversed = []

        complete_co_oriented_sols = (
            [s for s in self.co_oriented_solitons if s.is_complete is True]
        )
        complete_anti_oriented_sols = (
            [s for s in self.anti_oriented_solitons if s.is_complete is True]
        )

        for sol_a in complete_co_oriented_sols:
            for sol_b in complete_anti_oriented_sols:
                self.closed_solitons.append(
                    ClosedSoliton(
                        label=sol_a.label + '_+_' + sol_b.label, 
                        soliton_a=sol_a, 
                        soliton_b=sol_b, 
                        # soliton_a=sol_b, 
                        # soliton_b=sol_a, 
                        network=self.network,
                        resolution=self.resolution
                    )
                )
                self.closed_solitons_reversed.append(
                    ClosedSoliton(
                        label=sol_a.label + '_+_' + sol_b.label, 
                        # soliton_a=sol_a, 
                        # soliton_b=sol_b, 
                        soliton_a=sol_b, 
                        soliton_b=sol_a, 
                        network=self.network,
                        resolution=self.resolution
                    )
                )
                
        
        self.Q = 1
        self.Q_y = 1
        self.Q_y_reversed = 1
        for c_sol in self.closed_solitons:
            self.Q += c_sol.homology_class.symbol
            self.Q_y += (
                (symbols('y') ** c_sol.writhe) * c_sol.homology_class.symbol
            )
        for c_sol in self.closed_solitons_reversed:
            self.Q_y_reversed += (
                (symbols('y') ** c_sol.writhe) * c_sol.homology_class.symbol
            )

    def compute_open_solitons(self):
        """
        Compute the generating function of complete open solitons
        of both types on a 2-way streets, including their writhe.
        """
        # TODO
        pass


class NetworkSolitonContent:
    """
    Container class for a Network's soliton content
    """
    def __init__(self, network=None, iterations=None):
        self.american_data = None
        self.british_data = None
        self.network = network
        self.iterations = iterations

        if (self.network is not None) and (self.iterations is not None):
            self.compute()

    def compute(self):
        self.american_data = (
            [SolitonData(
                label='Q_('+street_label+')', 
                network=self.network, 
                street=street, 
                resolution='american'
            ) for street_label, street in self.network.streets.iteritems()]
        )

        self.british_data = (
            [SolitonData(
                label='Q_('+street_label+')', 
                network=self.network, 
                street=street, 
                resolution='british'
            ) for street_label, street in self.network.streets.iteritems()]
        )

        for q in self.american_data:
            q.initialize()
            q.grow(n_steps=self.iterations)

        for q in self.british_data:
            q.initialize()
            q.grow(n_steps=self.iterations)

    def print_info(self, full_path=False, soliton_paths=False, writhes=False):
        for q in self.american_data:
            q.print_info(
                full_path=full_path, 
                soliton_paths=soliton_paths, 
                writhes=writhes
            )

        for q in self.british_data:
            q.print_info(
                full_path=full_path, 
                soliton_paths=soliton_paths, 
                writhes=writhes
            )







