import sympy, numpy

class Seed():
    """
    A class for cluster seed.

    'quiver' is a numpy matrix 
    'nodes' is a list of labels for nodes, like ['gamma_1', ...]
    """
    def __init__(
        self, quiver=None, nodes=None, y_variables=None,
    ):
        self.quiver = quiver
        self.nodes = nodes
        self.dim = len(nodes)
        # create the y-variables related to the quiver, 
        # by default they are just [y_1, ...] when the seed is created 
        # otherwise they are specified from a previous seed
        if y_variables is None:
            self.y_variables = sympy.symbols('y_1:' + str(self.dim + 1))
        else:
            self.y_variables = y_variables

    def copy(self):
        return Seed(
            quiver=self.quiver, 
            nodes=self.nodes, 
            y_variables=self.y_variables
        )

    def mutate(self, node):
        # which node of the quiver we mutate on
        k = self.nodes.index(node)
        new_quiver = self.quiver.copy()
        new_y_variables = [y_v for y_v in self.y_variables]

        # mutate the quiver
        for r in range(self.dim):
            for c in range(self.dim):
                if r == k or c == k:
                    new_quiver[r, c] = -1 * self.quiver[r, c]
                else:
                    new_quiver[r, c] = self.quiver[r, c] + (
                        max(self.quiver[r, k], 0) * max(self.quiver[k, c], 0)
                        - max(self.quiver[c, k], 0) * max(self.quiver[k, r], 0)
                    )
        
        # mutate the y-variables
        for i in range(self.dim):
            if i != k:
                new_y_variables[i] = sympy.simplify(
                    self.y_variables[i] * (
                        self.y_variables[k] ** max(self.quiver[i, k], 0)
                    ) * (
                        1 + (self.y_variables[k] ** (-1))
                    ) ** self.quiver[i, k]
                )
            elif i == k:
                new_y_variables[i] = self.y_variables[i] ** (-1)

        self.quiver = new_quiver
        self.y_variables = new_y_variables

    def permute(self, perm_dic):
        """
        the permutation dictionary is of the form 
        {..., 'gamma_i' : 'gamma_j', }
        meaning that we send \gamma_i to \gamma_j
        """
        new_quiver = self.quiver.copy()
        # convert the permutation dictionary in the basis of the quiver
        # for ex {..., gamma_i : gamma_j ...}
        # translates into
        # {... i : j ...}
        node_dic = {
            self.nodes.index(key) : self.nodes.index(value) 
            for key, value in perm_dic.iteritems()
        }
        for r in range(self.dim):
            for c in range(self.dim):
                new_quiver[node_dic[r], node_dic[c]] = self.quiver[r, c]

        # build the inverse dictionary
        node_dic_inv = {val : key for key, val in node_dic.iteritems()}
        new_y_variables = [self.y_variables[node_dic_inv[i]] for i in range(self.dim)]
        
        self.quiver = new_quiver
        self.y_variables = new_y_variables

    def print_info(self):
        print '\nBasis cycles:\n{}'.format(self.nodes)
        print 'Quiver pairing matrix:\n{}'.format(self.quiver)
        print 'Cluster variables:\n{}'.format(self.y_variables)
  
### The matrix for the [3,1] puncture
# Q = numpy.matrix([
#     [0, 1, -1, -1, -1, 1, 1],
#     [-1, 0, 1, 1, -1, 0, 0],
#     [1, -1, 0, 1, -1, 1, -1],
#     [1, -1, -1, 0, 1, -1, 1],
#     [1, 1, 1, -1, 0, -1, -1],
#     [-1, 0, -1, 1, 1, 0, 0],
#     [-1, 0, 1, -1, 1, 0, 0]
#     ])

# G = ['gamma_'+str(i + 1) for i in range(7)]


# S = Seed(quiver=Q, nodes=G)
# print S.nodes
# print S.y_variables
# print S.quiver

# # S.mutate('gamma_1')
# S.permute({
#     'gamma_1' : 'gamma_3', 
#     'gamma_3' : 'gamma_1', 
#     'gamma_2' : 'gamma_2',
#     'gamma_4' : 'gamma_4',
#     'gamma_5' : 'gamma_5',
#     'gamma_6' : 'gamma_6',
#     'gamma_7' : 'gamma_7',
#     })
# print S.y_variables
# print S.quiver
# print S.nodes



# ### The matrix for the [2,1] ppuncture
# Q = numpy.matrix([
# [0, -1, 1, 1, -1],
# [1, 0, 1, -3, 1],
# [-1, -1, 0, 1, 1],
# [-1, 3, -1, 0, -1],
# [1, -1, -1, 1, 0]
# ])

# G = ['gamma_'+str(i + 1) for i in range(5)]

# S_inv_sequence = ['gamma_1', 'gamma_3', 'gamma_4', 'gamma_1']
# S_sequence = S_inv_sequence[::-1]
# S_inv_perm = {
#     'gamma_1' : 'gamma_3',
#     'gamma_2' : 'gamma_5',
#     'gamma_3' : 'gamma_2',
#     'gamma_4' : 'gamma_1',
#     'gamma_5' : 'gamma_4',
# }
# S_perm = {val : key for key, val in S_inv_perm.iteritems()}

# a = Seed(quiver=Q, nodes=G)
# a.print_info()

# def S_inv(seed):
#     new_seed = seed.copy()
#     for g in S_inv_sequence:
#         new_seed.mutate(g)
#     new_seed.permute(S_inv_perm)
#     return new_seed

# def S(seed):
#     new_seed = seed.copy()
#     for g in S_sequence:
#         new_seed.mutate(g)
#     new_seed.permute(S_perm)
#     return new_seed

# b = S(S_inv(a))
# b.print_info()

