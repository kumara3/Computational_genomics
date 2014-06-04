import argparse
import PhyloTree
import PhyloTree_sol
import PhyloTreeAlgs
import PhyloTreeAlgs_sol
from PhyloTree_util import *
from parseNewick_sol import *

def isPhyloTreeType(t):
    try:
        if isinstance(t, PhyloTree):
            return True;
    except:
        pass

    try:
        if isinstance(t, PhyloTree.PhyloTree):
            return True;
    except:
        pass

    try:
        if isinstance(t, PhyloTree_sol.PhyloTree):
            return True;
    except:
        pass

    return False


def _reorderTree(T):
    if T.isLeaf():
        T.sort_key = T.label
    else:
        for t,w in T.children:
            _reorderTree(t)
        T.children.sort(key = lambda n: n[0].sort_key)
        T.sort_key = "|".join([t.sort_key for t,w in T.children])

def compareTree(t1, t2, ignoreWeight = False):
    """Compare trees; ignore interior labels if userInteriorLabels
    is false"""
    
    _reorderTree(t1)
    _reorderTree(t2)

    return _compareTreeHelper(t1, t2, ignoreWeight)

def _compareTreeHelper(t1, t2, ignoreWeight):
    if t1.isLeaf() and t2.isLeaf():
        return t1.label == t2.label
    if t1.isLeaf() or t2.isLeaf() or t1.numChildren() != t2.numChildren():
        return False

    return all([_compareTreeHelper(r1,r2,ignoreWeight) and (ignoreWeight or abs(w2 - w1) < 0.0001) for (r1,w1),(r2,w2) in zip(t1.children, t2.children)])

def treeMatchesDist(t,D):
    """Returns True if D is the distance matrix for t"""
    D2 = createDistMatrix(t)

    if set(D.keys()) != set(D2.keys()):
        return False

    return all([abs(D2[s1][s2] - D[s1][s2]) < 0.0001 for s1 in D.keys() for s2 in D.keys() if s2 > s1])

def isPhyloTree(t):
    """Check that t is a legitimate PhyloTree object with unique leaf labels"""
    S = set()   # Leaf nodes
    return _isPhyloTreeHelper(t, S)

def _isPhyloTreeHelper(t, S):
    if not isPhyloTreeType(t):
        return False

    if t.isLeaf():
        if t.label in S:
            return False
        S.add(t.label)
        return True
    return all([_isPhyloTreeHelper(r,S) and (isinstance(w,int) or isinstance(w,float)) for r,w in t.children])

def translateTree(t):
    """Convert any PhyloTree to a PhyloTree_sol.PhyloTree object"""
    return PhyloTree_sol.PhyloTree(label = t.label, children = [(translateTree(c),w) for c,w in t.children])

#################################################
def _test_UPGMA(newick):
    """
    Input: 
    * newick: A newick string
    * student_alg: The student's UPGMA or NJ function
    * instructor_alg: The instructors UPGAM or NJ function
    Scoring: Let T be the PhyloTree return by the student's algorithm.
    * If T is not a correctly constructed PhyloTree: 0/1
    * If T is missing leaves or has extra leaves: 0/1
    * If createDistMatrix(t) is correct: 1/1
    * If createDistMatrix(t) is incorrect but T has the same shape as the insturctor solution: 0.6/1
    """
    t = parseNewick(newick)
    T = PhyloTreeAlgs.UPGMA(createDistMatrix(t))

    if not isPhyloTree(T):
        return 0
    
    if not ({x.label for x in T.leafSet()} == {x.label for x in t.leafSet()}):
        return 0

    if treeMatchesDist(T, createDistMatrix(t)):
        return 1

    T_sol = PhyloTreeAlgs_sol.UPGMA(createDistMatrix(t))

    if compareTree(T, T_sol, True):
        return 0.6

    return T, T_sol
    return 0





def test1():
    return _test_UPGMA("(A:1,B:1);")

def test2():
    return _test_UPGMA("((A:1,B:1):2,C:3);")

def test3():
    return _test_UPGMA("(C:3,(A:1,B:1):2);")

def test4():
    return _test_UPGMA("((((A:1,B:1):1,C:2):1,D:3):1,E:4);")

def test5():
    return _test_UPGMA("(((A:1,B:1):9,C:10):5,(D:5,(E:2,F:2):3):10);")


def _test_NJ(newick, outgroup):
    """
    Input: 
    * newick: A newick string
    * student_alg: The student's UPGMA or NJ function
    * instructor_alg: The instructors UPGAM or NJ function
    Scoring: Let T be the PhyloTree return by the student's algorithm.
    * If T is not a correctly constructed PhyloTree: 0/1
    * If T is missing leaves or has extra leaves: 0/1
    * If createDistMatrix(t) is correct: 1/1
    * If createDistMatrix(t) is incorrect but T has the same shape as the insturctor solution: 0.6/1
    """
    t = parseNewick(newick)
    T = PhyloTreeAlgs.NJ(createDistMatrix(t))

    if not isPhyloTree(T):
        return 0
    
    if not ({x.label for x in T.leafSet()} == {x.label for x in t.leafSet()}):
        return 0

    if treeMatchesDist(T, createDistMatrix(t)):
        return 1

    T = translateTree(T)
    T_sol = PhyloTreeAlgs_sol.NJ(createDistMatrix(t))

    return T, T_sol

    if compareTree(T.rootByOutgroup(outgroup), T_sol.rootByOutgroup(outgroup), True):
        return 0.6

    return T, T_sol
    return 0


def test6():
    return _test_NJ("((A:1,B:1):1, C:1);", "C")

def test7():
    return _test_NJ("(((A:10,B:1):20,(C:2,D:20):22):1,X:50);", "X")

def test8():
    return _test_NJ("(X:50,((C:2,D:20):22,(A:10,B:1):20):1);", "X")

def test9():
    return _test_NJ("(((A:10,C:1):1,(B:1,D:10):1):1,E:5);", "E")

def test10():
    return _test_NJ("(((((A:10,C:1):1,(B:1,D:10):1):1,E:5):1,(((A2:10,C2:1):1,(B2:1,D2:10):1):1,E2:5):1):1,X:1);", "X")

###########################################

def test(i):
    try:
        return float(eval("test%d()" % (i)))
    except:
        return 0


def run_test(show_results = True):
    T = [test(i) for i in range(1,11)]
    if show_results:
        for i,v in enumerate(T):
            print("Test %d: %5.3f" % (i+1,v))
    score = sum(T)/float(len(T))
    print("SCORE: %5.2f" % (100*score))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Testing dnaSeq project')
    parser.add_argument('-d', '--display', action = "store_true", help = "display indvidual test results", default = False)
    args = parser.parse_args()
    run_test(show_results = args.display)


