"""Implementation of the algorithm for isometric reconciliation of two
unrooted trees.
"""

from collections import namedtuple
from tree import Tree, Point

# Auxiliary named tuple Result and its documentation
Result = namedtuple('Result',
                    'x1 x2 alphamin alphamax h q1 q2 u1 u2 betamin betamax')
Result.__doc__ += "Structure for storing the result of reconciliation for one edge of S."
Result.x1.__doc__ = "Left endpoint of an edge in S"
Result.x2.__doc__ = "Right endpoint of an edge in S"
Result.alphamin.__doc__ = "Minimum value of alpha so that S(x1,x2,alpha) is a possible root"
Result.alphamax.__doc__ = "Maximum value of alpha so that S(x1,x2,alpha) is a possible root"
Result.h.__doc__ = "How far above root of S to map root of G"
Result.q1.__doc__ = "Root in G corresponding to S(x1,x2,alphamin) written as point G(a,b,gamma1)"
Result.q2.__doc__ = "Root in G corresponding to S(x1,x2,alphamax) written as point G(a,b,gamma2)"
Result.u1.__doc__ = "Left endpoint of an edge in G"
Result.u2.__doc__ = "Right endpoint of an edge in G"
Result.betamin.__doc__ \
    = "Value of beta so that S(u1,u2,beta) is a possible rooting wrt S rooted in S(x1,x2,alphamin)"
Result.betamin.__doc__ \
    = "Value of beta so that S(u1,u2,beta) is a possible rooting wrt S rooted in S(x1,x2,alphamax); betamin <= betamax"

Potential = namedtuple('Potential', 'point h leaves')
Potential.__doc__ += " Structure for storing a potential rooting of G"
Potential.point.__doc__ = "Point q, which is a potential root in tree G"
Potential.h.__doc__ = "Height h correspodning to the potential rooting"
Potential.leaves.__doc__ \
    = "Set of representative leaves, one for each edge leaving q in the subtree"



class Reconciliation:
    """Class for storing data structures necessary for reconciliation.

    Attributes:
       g gene tree (instance of Tree)
       s species tree (instance of Tree)
       potential dictionary, with keys pairs node and its parent
          defining a rooted subtree in the species tree,
          value is the potential rooting
          of this subtree (instance of Potential).
          If no potential rooting exists, the key is not inserted.
          The potential rootings are computed only if distance check passes
       min_leaf dictionary with keys pairs node and its parent,
          defining a rooted subtree of the gene tree.
          The value is the smallest leaf in the rooted subtree of g
      solution list of structures Result, one for each edge of the species
          tree with a solution

    The reconciliation is computed directly in the constructor.
    """

    def __init__(self, g, s):
        """Constructor getting gene tree g, species tree s.

        The constructor computes the reconciliation.
        """

        # store the graphs and initialize dictionaries
        self.g = g
        self.s = s
        self.potential = {}
        self.min_leaf = {}
        self.solution = []

        # min_leaves and distance check
        self.compute_min_leaves()
        if self.check_distances():
            # potential rootings
            self.compute_potential()
            # potential roots in intermediete form
            intermediate = self.find_potential_roots()
            # final form of the solution
            self.solution = self.find_solution(intermediate)

    def compute_min_leaves(self):
        """Compute min_leaf for every rootes subtree of self.g"""
        for (node, parent) in self.s.subtree_traversal:
            children = self.s.children(node, parent)
            if len(children) == 0:
                leaf = node
            else:
                leaf = min([self.min_leaf[(child, node)]
                            for child in children])
            self.min_leaf[(node, parent)] = leaf

    def check_leaf_dist(self, leaves):
        """Check distances for every pair of leaves in the input list.

        Returns true of for every pair of leaves the distance in g
        is at east as big as their distance in s. This method
        is used within check_distances to check selected sets of leaves.
        """
        for i in range(len(leaves)):
            for j in range(i+1, len(leaves)):
                dist_s = self.s.node_dist(leaves[i], leaves[j])
                dist_g = self.g.node_dist(leaves[i], leaves[j])
                if dist_g < dist_s:
                    return False
        return True

    def check_distances(self):
        """Check selected leaf distances in the trees.

        For every pair of nodes adjacent to the same internal node,
        compare distance in S and in G between their minimal leaves.
        If for any pair is the distane in G smaller, return False.
        If all tests pass, return True.
        """
        for node in self.s.preorder:
            neighbors = self.s.neighbors[node]
            # set of minimal leaves for all neighbors if node
            leaves = [self.min_leaf[(neighbor, node)]
                      for neighbor in neighbors]
            if not self.check_leaf_dist(leaves):
                return False

        # special case for a tree consisting of one edge
        if len(self.s.neighbors) == 2:
            if not self.check_leaf_dist(list(self.s.preorder)):
                return False

        # no problem found, return True
        return True

    def test_root_on_edge(self, x1, x2):
        """Check if edge (x1,x2) from tree s contains a potential root.

        If so, return Result structure containing nodes
        x1, x2 in S, min and max alpha, h and points q1 and q2
        corresponding to the endpoints of the intervals of alpha. If
        not solution exists, returns None."""

        # retrieve potential rootings for the two subtrees
        if (x1, x2) not in self.potential or (x2, x1) not in self.potential:
            return None
        p1 = self.potential[(x1, x2)]
        p2 = self.potential[(x2, x1)]

        dist_g = self.g.point_dist(p1.point, p2.point)
        dist_s = self.s.node_dist(x1, x2)
        # compute value of h from Claim 4
        h = (dist_g - dist_s + p1.h + p2.h)/2
        # possible values of alpha to place root between p1.point and p2.point
        alphamin = max(0, p1.h - h)
        alphamax = min(dist_s, p1.h-h+dist_g)

        def restrict_alpha(beta):
            """Auxiliary function to restrict alpha to a single beta."""
            nonlocal alphamin, alphamax, h, p1
            alpha = beta - h + p1.h
            alphamin = max(alphamin, alpha)
            alphamax = min(alphamax, alpha)

        # check if path p1,point - p2.point outside of subtree of G
        # containing leaves of subtree rooted at x1
        if not self.check_condition3_pair(p1.leaves, p1.point, p2.point):
            # root of G must be in p1.point, i.e. beta=0
            restrict_alpha(beta=0)

        # check if path p1,point - p2.point outside of subtree of G
        # containing leaves of subtree rooted at x2
        if not self.check_condition3_pair(p2.leaves, p2.point, p1.point):
            # root must be in p2.point, i.e. beta=dist_g
            restrict_alpha(beta=dist_g)

        # if we end up with an empty interval - no solution
        if alphamin > alphamax:
            return None

        beta1 = h - p1.h + alphamin
        beta2 = h - p1.h + alphamax
        q1 = self.g.get_point(p1.point, p2.point, beta1)
        q2 = self.g.get_point(p1.point, p2.point, beta2)

        return Result(x1, x2, alphamin, alphamax, h, q1, q2, None, None, None, None)


    def find_potential_roots(self):
        """Returns the list of edges in S with potential roots,
        and for each edge also interval of alpha, value of h and
        endpoints of the corresponding interval in G.

        However, edge in G is not located yet.  The result is a list
        of Result tuples in preorder traversal order. Values , with
        values of u1, u2, betamin, betamax are missing, to be supplied
        later.
        """

        # traverse all edges in preorder and check each of them if it
        # contains a potential root. Store found solutions
        results = []
        for node in self.s.preorder:
            parent = self.s.parent[node]
            if parent is not None:
                res = self.test_root_on_edge(node, parent)
                if res is not None:
                    results.append(res)
        return results

    def find_solution(self, results):
        """Get preliminary results, compute the final results, after locating
        the solution in tree G.

        List results contain records of type Result, but with values
        of u1, u2, betamin and betamax missing. The return value is a
        similar list with the missing values filled in. Edges of S are
        listed in results in the preorder traversal ordering.
        """

        def find_edge(edges, point1, point2):
            """Given a set of edges, find and return the one that contains
            both point1 and point2.

            Such an edge should exist.
            """
            for edge in edges:
                if self.g.between(Point(edge[0]), point1, Point(edge[1])) \
                   and self.g.between(Point(edge[0]), point2, Point(edge[1])):
                    return edge
            assert False

        def complete_result(result, edge):
            """Compute missing values in the result record,
            assuming that points q1 and q2 are located on a given edge in G.

            Returns a new record with missing values filled in."""

            # compute values of beta for q1 and q2 in the result
            betamin = self.g.point_dist(Point(edge[0]), result.q1)
            betamax = self.g.point_dist(Point(edge[0]), result.q2)
            # if betamin>betamax, use the edge in opposite order to fix this
            if betamin > betamax:
                return complete_result(result, [edge[1], edge[0]])
            # add newly computed values to the record
            # u1 and u2 are simply endpoints of the edge
            return result._replace(u1=edge[0],
                                   u2=edge[1],
                                   betamin=betamin,
                                   betamax=betamax)

        def get_edges(point):
            """Return the list of all edges of G adjacent to a given point.

            We assume that the point is given as G(u,v,alpha)
            for u and v adjacent."""

            dist = self.g.node_dist(point.node1, point.node2)
            # in point is one of the endpoints of the edge,
            # return all edges adjacent to that endpoint
            if point.alpha == 0:
                return [(point.node1, other)
                        for other in self.g.neighbors[point.node1]]
            if point.alpha == dist:
                return [(point.node2, other)
                        for other in self.g.neighbors[point.node2]]
            # otherwise if point is inside an edge, return that edge
            return [(point.node1, point.node2)]

        # body of find_solution
        # nodemap maps nodes of S to potential rootings in G
        nodemap = {}
        # new_results is a list of new records with missing values filled in
        new_results = []
        for i, result in enumerate(results):
            # first edge is found among all edges of G
            if i == 0:
                edges = self.g.edges
            else:
                # subsequent edges are adjacent to a known point from nodemap
                # search only among edges adjacent to that point
                if result.x1 in nodemap:
                    edges = get_edges(nodemap[result.x1])
                else:
                    assert result.x2 in nodemap
                    edges = get_edges(nodemap[result.x2])

            # find the edge with solution among selected edges
            edge = find_edge(edges, result.q1, result.q2)

            # complete and sotre the Result record
            new_result = complete_result(result, edge)
            new_results.append(new_result)
            # if the solution covers a node in S, add it to nodemap
            x_dist = self.s.node_dist(new_result.x1, new_result.x2)
            if new_result.alphamin == 0:
                nodemap[new_result.x1] = Point(new_result.u1, new_result.u2,
                                               new_result.betamin)
            if new_result.alphamax == x_dist:
                nodemap[new_result.x2] = Point(new_result.u1, new_result.u2,
                                               new_result.betamax)

        return new_results


    def print_solution(self, filename):
        """Print records stored in solution attribute to filename"""

        with open(filename, "w") as f:
            for result in self.solution:
                f.write("edge in S (%s,%s) alpha [%g,%g]; edge in G (%s,%s) beta [%g,%g]; h %g\n"
                  % (result.x1, result.x2, result.alphamin, result.alphamax,
                     result.u1, result.u2, result.betamin, result.betamax,
                     result.h))
            if len(self.solution)==0:
                f.write("no solution")


    def compute_potential(self):
        """For each rooted subtree of S, compute the potential
        rooting of the corresponding part of G."""

        # in bottom-up traversal of rooted subtrees of S
        for (node, parent) in self.s.subtree_traversal:
            children = self.s.children(node, parent)
            # special case for leaves
            if len(children) == 0:
                self.potential[(node, parent)] \
                    = Potential(Point(node), 0, {node})
            else:
                # internal vertices
                assert len(children) == 2
                self.potential_connect_two(node, children[0],
                                           children[1], parent)

    def potential_connect_two(self, node, child1, child2, parent):
        """Compute a potential rooting for a rooted subtree of S rooted at
        node, given subtrees for its two children.

        The rooting is stored in the dictionary potential, if it exists.
        """

        # if any of the children do not have potential rooting, no solution
        if (child1, node) not in self.potential:
            return
        if (child2, node) not in self.potential:
            return

        # get potential rootings for subtrees
        p1 = self.potential[(child1, node)]
        p2 = self.potential[(child2, node)]
        dist_g = self.g.point_dist(p1.point, p2.point)
        dist_s = self.s.node_dist(child1, child2)

        # compute h using Condition 1 of Claim 4
        h = (dist_g - dist_s + p1.h + p2.h)/2
        # compute beta using Condition 2 of Claim 4
        beta = h + self.s.node_dist(child1, node) - p1.h
        # if beta not within path p1, p2, no solution
        if beta < 0 or beta > dist_g:
            return

        # create the rooting
        point = self.g.get_point(p1.point, p2.point, beta)
        # check Condition 3 of Lemma 4 using representative elaves
        leaves = self.check_condition3(point, p1, p2)
        if leaves is None:
            return

        # store the potential rooting
        self.potential[(node, parent)] = Potential(point, h, leaves)


    def check_condition3(self, point, p1, p2):
        """Check if potential root at point satisfies Condition 3 in Lemma 4
        with respect to potential rootings p1 and p2 of the two subtrees.

        Return the set of representative leaves for point, or None if
        Condition 3 not satisfied.
        """

        def first_element(some_set):
            """Auxiliary function for getting the first element ina  set."""
            for a in some_set:
                return a

        leaves = None
        if self.g.point_dist(p1.point, p2.point) == 0:
            # case point == p1.point == p2.point
            # take all leaves from p1
            # add leaves from p2 that do not connect via the same edge
            leaves = set(p1.leaves)
            for a2 in p2.leaves:
                add = True
                for a1 in p1.leaves:
                    if not self.g.between(Point(a1), point, Point(a2)):
                        # a1 and a2 connect to point via the same edge
                        add = False
                if add:
                    leaves.add(a2)
        elif self.g.point_dist(point, p1.point) == 0:
            # point == p1.point != p2.point
            if self.check_condition3_pair(p2.leaves, p2.point, p1.point):
                leaves = set(p1.leaves).union({first_element(p2.leaves)})
        elif self.g.point_dist(point, p2.point) == 0:
            # point == p2.point != p1.point
            if self.check_condition3_pair(p1.leaves, p1.point, p2.point):
                leaves = set(p2.leaves).union({first_element(p1.leaves)})
        else:
            # point, p1.point, p2.point distinct
            if self.check_condition3_pair(p2.leaves, p2.point, p1.point) and \
               self.check_condition3_pair(p1.leaves, p1.point, p2.point):
                leaves = {first_element(p1.leaves), first_element(p2.leaves)}
        return leaves


    def check_condition3_pair(self, leaves1, point1, point2):
        """Check Condition 3 from Claim 4 when q is not point1.

        Points point1 and point2 are rootes of two subtrees in G to be
        connected. Set leaves1 contains representative leaves for all
        edges leaving point1 within its subtree. Chech that for each
        of these leaves, the path from the leaf to point2 leads
        through point1.
        """
        for a in leaves1:
            if not self.g.between(Point(a), point1, point2):
                return False
        return True

def read_tree(filename):
    edges = []
    with open(filename) as f:
        for line in f:            
            parts = line.split()
            parts[2] = int(parts[2])
            edges.append(parts)
    return edges
    

def main():
    """A simple main function which runs reconciliations on several pairs
    of trees read from examples directory and writes the results to
    the same directory.
    """

    import glob
    for filename_g in sorted(glob.glob('examples/input*-g.txt')):
        filename_s = filename_g.replace("-g.txt", "-s.txt")
        print(filename_g, filename_s)
        g = Tree(read_tree(filename_g))
        s = Tree(read_tree(filename_s))
        r = Reconciliation(g, s)
        filename_out = filename_g.replace("-g.txt", "-out.txt")
        r.print_solution(filename_out)

if __name__ == '__main__':
    main()
