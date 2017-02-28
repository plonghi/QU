# reference:
# https://networkx.readthedocs.io/en/stable/reference/drawing.html#module-networkx.drawing.nx_pylab
# bit NOTE: cannot rely on automatic positioning, because the ordering of edges at branch points actually matters!
# so instead you could employ an interative applet that prompts you to choose positions of nodes, and lets you drag them, and then positions can be read off from there.

from mcsn import MCSN, Street, Joint, BranchPoint
import networkx as nx
import matplotlib.pyplot as plt


def draw_network(w):
	G = nx.Graph()
	# edge labels
	edge_labels = w.streets.keys()
	# node labels	
	temp = w.branch_points.copy()
	temp.update(w.joints)
	node_labels = dict((v,k) for k,v in temp.iteritems())

	node_list = w.branch_points.values() + w.joints.values()
	G.add_nodes_from(node_list)
	edge_list = [[ep.end_point for ep in s.endpoints] for s in w.streets.values()]
	G.add_edges_from(edge_list)

	print node_labels
	print edge_labels
	print node_list
	print edge_list

	nx.draw_circular(G)
	nx.draw_networkx_labels(
		G, labels=node_labels, pos=nx.circular_layout(G),
		font_size=16
	)
	# nx.draw_networkx_edge_labels(
	# 	G, labels=edge_labels, pos=nx.circular_layout(G)
	# )
	plt.show()

	return G

### SU(2) N_f=4

# ----- Create a spectral network ------
# 	  	   SU(2) N_f=4 network


streets = ['p_'+str(i+1) for i in range(6)]
branch_points = {
	'b_1': ['p_1', 'p_2', 'p_3'],
	'b_2': ['p_2', 'p_4', 'p_5'],
	'b_3': ['p_1', 'p_6', 'p_4'],
	'b_4': ['p_3', 'p_5', 'p_6'],
}

joints = {}

homology_classes = {
	'gamma_1' : ['p_1'],
	'gamma_2' : ['p_2'],
	'gamma_3' : ['p_3'],
	'gamma_4' : ['p_4'],
	'gamma_5' : ['p_5'],
	'gamma_6' : ['p_6'],
}

w = MCSN(
	branch_points=branch_points, 
	streets=streets, 
	joints=joints, 
	homology_classes=homology_classes
)

draw_network(w)

# G = nx.Graph()
# nodes_list = w.branch_points.values() + w.joints.values()
# G.add_nodes_from(nodes_list)
# edge_list = [[ep.end_point for ep in s.endpoints] for s in w.streets.values()]
# G.add_edges_from(edge_list)

# nx.draw_random(G)
# plt.show()
# print nodes_list
# print edge_list
