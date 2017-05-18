import sympy, numpy

class Seed():
    """
    A class for cluster seed.

    'quiver' is a numpy matrix 
    'nodes' is a list of labels for nodes, like ['gamma_1', ...]
    """
    def __init__(
        self, quiver=None, nodes=None,
    ):
        self.quiver = quiver
        self.nodes = nodes
        self.dim = len(nodes)
        # create the y-variables related to the quiver, 
        # when the class is created
        self.y_variables = sympy.symbols('y_1:' + str(self.dim + 1))

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
                new_y_variables[i] = self.y_variables[i] * (
                    self.y_variables[k] ** max(self.quiver[i, k], 0)
                ) * (
                    1 + (self.y_variables[k] ** (-1))
                ) ** self.quiver[i, k]
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

        new_y_variables = [self.y_variables[node_dic[i]] for i in range(self.dim)]
        
        self.quiver = new_quiver
        self.y_variables = new_y_variables
    
Q = numpy.matrix([
    [0, 1, -1, -1, -1, 1, 1],
    [-1, 0, 1, 1, -1, 0, 0],
    [1, -1, 0, 1, -1, 1, -1],
    [1, -1, -1, 0, 1, -1, 1],
    [1, 1, 1, -1, 0, -1, -1],
    [-1, 0, -1, 1, 1, 0, 0],
    [-1, 0, 1, -1, 1, 0, 0]
    ])

G = ['gamma_'+str(i + 1) for i in range(7)]


S = Seed(quiver=Q, nodes=G)
print S.nodes
print S.y_variables
print S.quiver

# S.mutate('gamma_1')
S.permute({
    'gamma_1' : 'gamma_3', 
    'gamma_3' : 'gamma_1', 
    'gamma_2' : 'gamma_2',
    'gamma_4' : 'gamma_4',
    'gamma_5' : 'gamma_5',
    'gamma_6' : 'gamma_6',
    'gamma_7' : 'gamma_7',
    })
print S.y_variables
print S.quiver
print S.nodes