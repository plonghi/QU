"""
Module to handle soliton generating functions
"""

from solitons import SolitonPath 
from growth_rules import grow_soliton

# Default maximal growth level for solitons.
N_MAX_GROWTH = 10

class SolitonData:
	"""
	A container for soliton data supported on a street.
	For a two-way street of type ij/ji it will collect 
	all the solitons of each type.
	It will also provide a functionality for growing solitons,
	either by a number of steps, or until all of them are 'complete'.
	"""
	def __init__(self, label=None, network=None, street=None):
		self.label = label
		self.network = network
		self.street = street
		
		# given the orientation of the street, let that be the 
		# orientation of the underlying one-way ij-street. 
		# Then co-oriented solitons will be of ij type, while
		# anti-oriented solitons will be of ji type.
		self.co_oriented_solitons = []
		self.anti_oriented_solitons = []

	def initialize(self):
		# The initial co-oriented soliton
		ij_sol = SolitonPath(label=self.label+'ij_sol')
		ij_source_pt = self.street.initial_point().end_point
		ij_source_slot = self.street.initial_point().slot
		ij_sol.create(
			street=self.street, 
			source_pt=ij_source_pt, 
			slot=ij_source_slot
		)

		# The initial anti-oriented soliton
		ji_sol = SolitonPath(label=self.label+'ji_sol')
		ji_source_pt = self.street.final_point().end_point
		ji_source_slot = self.street.final_point().slot
		ji_sol.create(
			street=self.street, 
			source_pt=ji_source_pt, 
			slot=ji_source_slot
		)

		self.co_oriented_solitons = [ij_sol]
		self.anti_oriented_solitons = [ji_sol]

	def grow(self, n_steps=N_MAX_GROWTH):
		# grow co-oriented solitons
		new_co_or_sols = []
		for sol in self.co_oriented_solitons:
			new_sols = grow_soliton(sol, n_steps=n_steps)
			new_co_or_sols += new_sols

		self.co_oriented_solitons = new_co_or_sols

		# grow anti-oriented solitons
		new_anti_or_sols = []
		for sol in self.anti_oriented_solitons:
			new_sols = grow_soliton(sol, n_steps=n_steps)
			new_anti_or_sols += new_sols

		self.anti_oriented_solitons = new_anti_or_sols



