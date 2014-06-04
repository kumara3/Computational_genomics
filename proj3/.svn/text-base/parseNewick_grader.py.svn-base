import os
import PhyloTree_sol
import parseNewick_sol
import parseNewick
import PhyloTree
import argparse
from random import *
import re

labels = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

def randomTree(maxLevels, minLeafLevels, maxChildren, weightRange, level=0, usedLabels = set()):
    """Generate a random tree"""
    numChildren = randint(1,maxChildren) if level < minLeafLevels else (randint(0,maxChildren) if level < maxLevels else 0)

    if numChildren > 0:
        childNodes = [randomTree(maxLevels, minLeafLevels, maxChildren, weightRange, level+1, usedLabels) for i in range(numChildren)]
        weights = [int(uniform(*weightRange)) if random() < 0.25 else round(uniform(*weightRange), randint(0,4)) for i in range(numChildren)]
        return PhyloTree_sol.PhyloTree(children = list(zip(childNodes, weights)))
    
    while True:
        label = "".join([choice(labels) for i in range(3)])
        if not label in usedLabels:
            break
        
    usedLabels.add(label)

    return PhyloTree_sol.PhyloTree(label=label)

def addLabels(t):
    """Give a label to every node"""
    if not t.label:
        t.label = "".join([choice("abcdefghijklmnopqrstuvwxyz") for i in range(4)])
    for r,w in t.children:
        addLabels(r)

def intEdge(t):
    """Round every edge weight to an int"""
    t.children = [(r,int(w)) for r,w in t.children]
    for r,w in t.children:
        intEdge(r)


########

def reorderTree(T):
    if T.isLeaf():
        T.sort_key = T.label
    else:
        for t,w in T.children:
            reorderTree(t)
        T.children.sort(key = lambda n: n[0].sort_key)
        T.sort_key = "|".join([t.sort_key for t,w in T.children])

def compareTree(t1, t2):
    """Compare trees; ignore interior labels if userInteriorLabels
    is false"""
    
    reorderTree(t1)
    reorderTree(t2)

    return compareTreeHelper(t1, t2)


# Error set return value:
#   0 = bad match; no credit
#   1 = interior label match error
#   2 = approx. weight error
#   3 = large weight error
#   4 = bad leaf label
def compareTreeHelper(t1, t2):
    S = set()
    if t1.isLeaf() and t2.isLeaf():
        if t1.label != t2.label:
            S.add(4)   # Leaf labels must match
    elif t1.isLeaf() or t2.isLeaf():
        S.add(0)   # Tree structure must match
    else:  # Now at two interiod nodes
        for (r1,w1),(r2,w2) in zip(t1.children,t2.children):
            S |= compareTreeHelper(r1,r2)
            if 0 in S:
                break
            if abs(float(w1) - float(w2)) > 0.001:
                if abs(int(w1) - int(w2)) <= 1:
                    S.add(2)
                else:
                    S.add(3)
        if t1.label != t2.label:
            S.add(1)
    return S
    


def reduceLabel(T):
    """Reduce all interior node labels to a single character"""
    if not T.isLeaf():
        T.label = T.label[:1]
        for t,w in T.children:
            reduceLabel(t)

def resetWeights(T):
    """Change all edge weights to 0"""
    T.children = [(t,0) for t in T.children]
    for t,w in T.children:
        resetWeights(t)


def checkStr(s):
    try:
        T1 = parseNewick_sol.parseNewick(s)
        T2 = parseNewick.parseNewick(s)

        S = compareTree(T1,T2)
        if 0 in S:
            return 0

        score = 1
        if 1 in S:
            score -= 0.15
        if 2 in S:
            score -= 0.1
        if 3 in S:
            score -= 0.25
        if 4 in S:
            score -= 0.25
        return score
    except Exception as e:
        return 0

def checkRandomTree():
    reload(parseNewick)
    T1 = randomTree(4,3,3,(1,4));

    T2 = T1.copy()
    addLabels(T2)   # ensure every node has a label

    T3 = T1.copy()
    intEdge(T3)   # ensure every edge has an integer weight

    T4 = T2.copy()               # ensure every node has a label and every edge has an integer weight
    intEdge(T4);

    S = [checkStr(str(t)) for t in [T1,T2,T3,T4]]
    return 0.8*max(S) + 0.2*( (sum(S) - max(S)) / (len(S)-1))
    
def run_test(n, show_results = True, RNGseed = None):
    T = ["A;", "(A:1);", "(A:1.2);", "(A:1)B;", "(A:1.2)B;", "(A:1,B:1)C;", "(A:1,B:1);", "(((A:1)B:1)C:1)D;", "(((A:1):1):1);", "(A:1,B:1,D:1)E;"] 
    
    # Tests 1-10: Check fixed trees
    P = [checkStr(s) for s in T]

    # Test 11: Check n random trees
    S = [checkRandomTree() for i in range(n)]
    P.append(sum(S) / len(S))

    if show_results:
        print("\n".join(["Test %d: %5.3f" % (i+1, v) for i,v in enumerate(P)]) + "\n")

    score = 0.8*sum(P[:-1])/float(len(P[:-1])) + 0.2*P[-1]
    print("SCORE: %5.2f" % (100*score))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Testing dnaSeq project')
    parser.add_argument('-n', type = int, help = "Each result will be averaged over n runs", default = 100)
    parser.add_argument('-s', type = int, help = "Seed for random number generator", default = None)
    parser.add_argument('-d', '--display', action = "store_true", help = "display indvidual test results", default = False)
    args = parser.parse_args()
    run_test(args.n, RNGseed = args.s, show_results = args.display)




