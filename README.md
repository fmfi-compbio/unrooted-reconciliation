# unrooted-reconciliation
A linear-time algorithm for isometric reconciliation of two unrooted trees
Algorithm by Broňa Brejová and Rastislav Královič, implementation by Broňa Brejová

This a simple proof-of-concept implementation, intended to illustrate the details of the algorithm, not for actual use. 

The input to the algorithm are two unrooted trees (a gene tree and a species tree) with branch lengths known exactly. The algorithm finds all pairs of positions where to root the trees so that the resulting reconciliation obeys the branch lengths. 

The implementation has the following limitations:

* The two trees are given in a simple format, as a list of edges, each edge given by two strings denoting node names and a real-valued weight
* It is assumed that preprocessing was already done, i.e. that the two trees have the same set of leaves, all internal nodes of the species tree have degree three, all internal nodes of the gene tree have degree at least three, all edge lengths are non-negative and the gene tree does not contain zero-length edges.
* The trees have at least two leaves each
* LCA computation is assumed in the algorithm to work in O(1) time, but only trivial O(n) implementation is provided in the code.
* For simplicity, the implementation frequently uses dctionaries; this could be avoided by numbering nodes and edges by 0,1,... and using direct adressing instead of hashing

Examples of input and output files are in the examples folder.

All internal nodes of tree S have degree three; S may have zero-length edges.
3. All internal nodes of tree G have degree at least three and bounded from above by some constant
c; all edges of G have strictly positive lengths.
