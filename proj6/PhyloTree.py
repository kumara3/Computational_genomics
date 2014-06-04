import re

class NewickParseError(Exception):
    def __init__(self, value = "Bad Newick String"):
        self.value = value

    def __str__(self):
        return str(self.value)

class TreeIterator:
    def __init__(self, t):
        self.stack = [t]

    def __iter__(self):
        return self

    def __next__(self):
        while self.stack:
            t = self.stack.pop()
            if t.isLeaf():
                return t;
            else:
                self.stack.extend([c for c,w in t.children])
        raise StopIteration

class PhyloTree:
    """Phylogenetic tree manipulation and Newick conversion"""

    # PhyloTree constructor.  You may add to this as needed, but may not change the code
    # provided.
    def __init__(self, label = "", children = []):
        """Contructor: create a root node with the specified label and child set.
        Assumes children is a list of tuples, with each tuple containing a 
        PhyloTree object and a weight."""
        if not label and not children:
            raise ValueError("Illegal tree: Empty arguments")
        if label and not type(label) == str:
            raise ValueError("Illegal tree: non-string label argument")
        if children and not type(children) == list:
            raise ValueError("Illegal tree: non-list children argument")
        self.label = label if label else ""
        self.children = children if children else []

    # 1) We will define the length of a tree as the number of leaves in the tree.
    def __len__(self):
        """Return the number of leaves in the tree."""
        return 1 if self.isLeaf() else sum([len(t) for t,w in self.children])

    def _strHelper(self):
        if self.isLeaf():
            return self.label
        else:
            return "(" + ",".join([t._strHelper() + ":" + str(w) for t,w in self.children]) + ")" + self.label

    # 2) Print the string in Newick format.
    def __str__(self):
        """Return the Newick representation of the tree"""
        return self._strHelper() + ";"

    # Uncomment this function once __str__ is written.
    def __repr__(self):
        return "PhyloTree: " + str(self)

    # EXTRA CREDIT: Iterator will visit all leaf nodes exactly once in an order of your choice.
    def __iter__(self):
        """Iterate throught the *leaves* of the tree is an order of your choice"""
        return TreeIterator(self)

    def rootLabel(self):
        """Return the label of the root."""
        return self.label

    def child(self, i):
        """Get the ith child"""
        return self.children[i][0]

    def numChildren(self):
        """Return the number of children"""
        return len(self.children)

    def isLeaf(self):
        """Indicates whether a tree has children"""
        return self.numChildren() == 0 

    # 3) Return a deep copy of the tree.  That is, the resulting tree should contain entirely
    # new node objects (so modifying a node in one tree does not modify nodes in the other), but
    # should be  isomorphic ("identical") to the source.
    def copy(self):
        """Return a deep copy of the tree"""
        return PhyloTree(label = self.label, children = [(t.copy(), w) for t,w in self.children])

    # 4) Return a set of all leaf nodes in the tree.  (The set should contain PhyloTree objects,
    # not just labels.)
    def leafSet(self):
        """Return a set of the leaf nodes"""
        return {self} if self.isLeaf() else {x for t,w in self.children for x in t.leafSet()}
    
    # 5) Return the sum of all edges in the tree.
    def treeWeight(self):
        """Return the sum of the weight of all edges in the tree"""
        return sum([r.treeWeight() + w for r,w in self.children])

    def mapDepthHelper(self, D, d):
        if self.isLeaf():
            D[self.label] = d
        else:
            [ t.mapDepthHelper(D, d+w) for t,w in self.children ]

    # 6) Return a dictionary that maps the labels of the leaves to their dpeth in the tree.
    def mapDepth(self):
        """Return a dictionary mapping leaf labels to their depth in the tree."""
        D = {}
        self.mapDepthHelper(D,0)
        return D

    # 7) Return a dictionary mapping each leaf-label of the tree to a list of
    # the nodes on the path from the root to the leaf (in order).
    def mapPath(self):
        D = {}
        path = []
        return self._mapPath(D, path)

    def _mapPath(self, D, path):
        if self.isLeaf():
            D[self.label] = path + [self]
        else:
            for t,w in self.children:
                t._mapPath(D, path + [self])
        return D

    # 8) Return the MRCA (Most Recent Common Ancestor) of the two labels.  That is,
    # returns the node with the greatest depth that is an ancestor to both leaves.  D, if provided,
    # is the result of the mapPath call (so that does not need to be re-calculated).
    # Throw a ValueError is one of the labels is not in the tree.
    def MRCA(self, label1, label2, D=None):
        if not D:
            D = self.mapPath()
        try:
            P1 = D[label1];
            P2 = D[label2];
        except:
            raise ValueError("MRCA: Label not contained in tree")

        i = 0;
        while i < len(P1) and P1[i] == P2[i]:
            i+=1

        return P1[i-1]
            
    def _find_child_index(self, y):
        """Return the i such that y is the ith child of self; return -1 is not a child"""
        for i,T in enumerate(self.children):
            if T[0] == y:
                return i
        return -1



    # 9) Returns a copy of the tree that has been reoriented based on an outgroup.
    #    Input:
    #    * outgroup_label: The label of one leaf that will be a child of the root in the
    #      reoriented tree. 
    #    * split_fraction: The split edge will have 100*split_fraction% of its weight
    #      distriubted to the r->outgroup edge, and the remaining distrubted to the other
    #      new edge.
    #    * root_label: The label of the new root node.
    def rootByOutgroup(self, outgroup_label, split_fraction = 0.5, root_label = ""):
        T = self.copy()
        MP = T.mapPath()
        P = MP[outgroup_label]

        # Flip edges along the path
        for x,y in zip(P[:-2],P[1:-1]):    # x and y will take on consecutive nodes along the path 
            # Flip directed edge (x,y) to (y,x)
            i = x._find_child_index(y)
            assert i >= 0
            w = x.children[i][1]          # record weight of x -> y edge
            x.children = x.children[:i] + x.children[i+1:]  # Remove y has a child of x
            y.children.append((x, w))     # Add x as a child of y
            
        x, y = P[-2:]              # These need to be the children of the new root
        i = x._find_child_index(y)
        w = x.children[i][1]
        x.children = x.children[:i] + x.children[i+1:]  # Remove y as a child of x
        return PhyloTree(label = root_label, children = [(x, split_fraction*w), (y, (1-split_fraction)*w)])


