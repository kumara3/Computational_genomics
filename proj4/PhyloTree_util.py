from PhyloTree import *

labels = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

def isFullBinary(T):
    """
    Input: A PhyloTree object T.
    Return: True if T is a binary tree, False otherwise.
    """
    return T.isLeaf() or (len(T.children) == 2 and all([isFullBinary(t) for t,w in T.children]))

def distFromRoot(T):
    """Return library mapping each node to the root-to-node path length."""
    D = {}
    _distFromRoot(T, D, 0)
    return D

def isUltraMetric(T):
    """
    Input: A PhyloTree object T.
    Return: True if T is an ultrametric tree, False otherwise.
    """
    return _isUltraMetric(T) >= 0

def _distFromRoot(T, D, d):
    D[T] = d
    for t,w in T.children:
        _distFromRoot(t, D, d+w)
    

def createDistMatrix(T):
    """
    Input: A PhyloTree object T.
    Return value: A disctionary D such that:
    * The keys of D are the leaf labeld of T.
    * For each leaf label pair x,y of T, D[x][y] is the path distance from x to y.
    """
    L = T.leafSet()
    DPath = T.mapPath()
    DDist= distFromRoot(T)
    
    D = {x.label:{} for x in L}
    for x in L:
        D[x.label][x.label] = 0
        for y in L:
            if x != y:
                A = T.MRCA(x.label, y.label, DPath)
                D[x.label][y.label] = DDist[x] + DDist[y] - 2*DDist[A]
    return D


def _isUltraMetric(T):
    if T.isLeaf():
        return 0
    if T.numChildren() != 2:
        return -1
    c1 = _isUltraMetric(T.child(0))
    c2 = -1 if c1 < 0 else _isUltraMetric(T.child(1))
    if c1 < 0 or c2 < 0:
        return -1
    if abs( (c1 + T.children[0][1]) - (c2 + T.children[1][1]) ) > 0.0001:
        return -1
    return c1 + T.children[0][1]


def sameTree(T1,T2):
    """
    Input: PhyloTree objects T1 and T2
    Return value: True if T1 and T2 represent the same tree with the same lecture
                  when disregarding orderings.
    Example:
    * sameTree(T1,T2) == True when T1 is (A:5,(B:10,C:20)D:1); and T2 is ((B:10,C:20)D:1,A:5);
    * sameTree(T1,T2) == False when T1 is (A:5,(B:10,C:20)D:1); and T2 is ((C:10,B:20)D:1,A:5);
    """
    return _sort_key(T) == _sort_key(T2)

def _sort_key(T):
    if T.isLeaf():
        return T.label
    return "(" + ",".join(sorted([t._strHelper() + ":" + str(w) for t,w in T.children])) + ")" + T.label
