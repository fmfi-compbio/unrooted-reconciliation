"""Module implementing a tree data structure for the isometric reconciliation algorithm."""

from collections import defaultdict

class Tree:
    """Class for storing a single unrooted tree and various data
    structures needed in the isometric reconciliation algorithm.

    Nodes are represented by their string names. The tree is rooted in
    an arbitrary node for simplifying some computations.
    We use a leaf as this auxiliary root.

    Attributes:

      edges: a dictionary mapping a pair of edges to its weight (length).
            Each edge is inserted twice, once in each direction.
      neighbors: a dictionary mapping node name to a list of its neighbors
      root: the name of a node serving as an auxiliary root of the tree
      depth: a dictionary mapping each node to its depth, i.e. weighted
            distance from the root
      parent: a dictionary mapping each node to its parent. Parent of the root is None.
      preorder: the list of all nodes in an preorder traversal of the tree
      subtree_traversal: list of all rooted subtrees of the unrooted tree in
          the bottom-up order, where subtrees of each subtree tree are listed before it.
    """

    def __init__(self, edges):
        """Constructor whch gets a list of edges, each edge a triple (start, end, weight).
        All necessary structures are initialized."""

        # process edges into self.edges and self.neighbors
        self.edges = {}
        self.neighbors = defaultdict(list)
        for edge in edges:
            assert len(edge) == 3
            pair = (edge[0], edge[1])
            assert pair not in edges
            self.edges[pair] = edge[2]
            pair = (edge[1], edge[0])
            assert pair not in edges
            self.edges[pair] = edge[2]
            self.neighbors[edge[0]].append(edge[1])
            self.neighbors[edge[1]].append(edge[0])

        #select a leaf as a root
        self.root = None
        for node in self.neighbors:
            if len(self.neighbors[node]) == 1:
                self.root = node
                break
        assert self.root is not None

        # initialize depth, parent and preorder by a recursive traversal
        self.depth = {}
        self.parent = {}
        self.preorder = []
        self.init_traversal(self.root, 0, None)

        #initialize subtree_traversal by two recursive traversals
        self.subtree_traversal = []
        self.postorder_subtree_traversal(self.root, None)
        self.preorder_subtree_traversal(self.root, None)

    def postorder_subtree_traversal(self, node, parent):
        """Method used to fill in the first part of subtree_traversal."""
        children = self.children(node, parent)
        for child in children:
            self.postorder_subtree_traversal(child, node)
        if parent is not None:
            self.subtree_traversal.append((node, parent))

    def preorder_subtree_traversal(self, node, parent):
        """Method used to fill in the second part of subtree_traversal."""
        children = self.children(node, parent)
        for child in children:
            self.subtree_traversal.append((node, child))
        for child in children:
            self.preorder_subtree_traversal(child, node)

    def init_traversal(self, node, depth, parent):
        """Method to recursively fill in preorder, depth and parent data structures."""
        if node in self.depth:  # attempt to call for parent
            return
        self.preorder.append(node)
        self.depth[node] = depth
        self.parent[node] = parent
        for node2 in self.neighbors[node]:
            self.init_traversal(node2, depth+self.edges[(node, node2)], node)

    def children(self, node, parent):
        """Return the set of children of a node given its parent.
        That is, return all neighbors except the parent."""
        assert node in self.neighbors
        children = list(self.neighbors[node])
        if parent is not None:
            assert parent in children
            children.remove(parent)
        return children

    def lca(self, node1, node2):
        """Find lca of two nodes (under the auxiliary rooting).
        This is here done by a simple linear-time algorithm"""

        def predecessors(node):
            """Get the list of predecessors of a given node starting in the root."""
            pred = []
            while node is not None:
                pred.append(node)
                node = self.parent[node]
            pred.reverse()
            return pred

        pred1 = predecessors(node1)
        pred2 = predecessors(node2)
        n = min(len(pred1), len(pred2))
        common = [pred1[i] for i in range(n) if pred1[i] == pred2[i]]
        return common.pop()

    def node_dist(self, node1, node2):
        """Compute distance of two nodes.
        This is accomplished using lca and depths."""
        return self.depth[node1] + self.depth[node2] \
            - 2*self.depth[self.lca(node1, node2)]

    def point_dist(self, point1, point2):
        """Compute distance between two points, each point specified on a path
        between two nodes."""

        # synchronize the points to the same pair of nodes
        (point1, point2) = point1.synchronize(point2, self)
        assert point1.node1 == point2.node1
        assert point1.node2 == point2.node2
        # distance is now absolute difference on their distances from the first point
        return abs(point1.alpha - point2.alpha)

    def get_point(self, point1, point2, alpha):
        """Create a point in distance alpha from point1
        along the path to point2."""

        # synchronize the points to the same pair of nodes
        (point1, point2) = point1.synchronize(point2, self)
        assert point1.node1 == point2.node1
        assert point1.node2 == point2.node2
        # decide if the new point is before or after point1
        dist = abs(point1.alpha - point2.alpha)
        assert alpha >= 0 and alpha <= dist
        if point1.alpha > point2.alpha:
            alpha = -alpha
        return Point(point1.node1, point1.node2,
                     point1.alpha + alpha)

    def between(self, point1, point2, point3):
        """Is point2 on the path between point1 and point3?"""
        dist13 = self.point_dist(point1, point3)
        dist12 = self.point_dist(point1, point2)
        dist23 = self.point_dist(point2, point3)
        return  dist13 == dist12 + dist23

class Point:
    """A simple class representing a pont on the tree given by a distance
    alpha from node node1 on the path towards node2.

    Attributes:
       alpha, node1, node2
    """

    def __init__(self, node1, node2=None, alpha=0):
        """A constructor storing values of all three attributes.

        If node2 is not specified, node1 as used as node2."""
        if node2 is None:
            node2 = node1
        self.node1 = node1
        self.node2 = node2
        self.alpha = alpha

    def __str__(self):
        """A simple string representation for debugging."""
        return "(%s,%s,%g)" % (str(self.node1), str(self.node2), self.alpha)

    def reverse(self, tree):
        """Return the same point represented with node1 and node2 exchanged."""
        dist = tree.node_dist(self.node1, self.node2)
        return Point(self.node2, self.node1, dist-self.alpha)

    def reroute(self, node, tree):
        """Return the same point represented so that node1 will be equal to
        node (node2 will stay the same or original node1 will be used
        as node2)."""

        # find where the path from node1 joins path node2-node.
        dist = self.dist_to_pair(self.node1, self.node2, node, tree)
        # if node joins before the point, we need to use node and node2
        # this is done by reversing the point and calling the same function
        if dist < self.alpha:
            result = self.reverse(tree).reroute(node, tree)
        else:
            # otherwise we use node1 and node, then reverse to get node1, node
            result = Point(self.node1, node, self.alpha).reverse(tree)
        return result

    @staticmethod
    def dist_to_pair(a, b, c, tree):
        """A function to compute the distance from a to path between b and c."""
        ab = tree.node_dist(a, b)
        ac = tree.node_dist(a, c)
        bc = tree.node_dist(b, c)
        return (ab+ac-bc)/2

    def synchronize(self, other, tree):
        """A method to recompute self and other point so that they use the
        same pair of nodes. A pair of new points is returned."""

        node = other.node1
        point = self.reroute(node, tree)
        dist = self.dist_to_pair(node, point.node2, other.node2, tree)
        if point.alpha <= dist:
            result = (Point(node, other.node2, point.alpha), other)
        elif other.alpha <= dist:
            result = (point, Point(node, point.node2, other.alpha))
        else:
            point1 = point.reverse(tree)
            point2 = other.reverse(tree)
            result = (Point(point1.node1, point2.node1, point1.alpha),
                      Point(point2.node1, point1.node1,
                            point2.alpha).reverse(tree))
        return result
