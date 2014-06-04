import argparse
import PhyloTree
import PhyloTree_sol
import parseNewick_sol
import random
import sys


labels = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

def randomTree(maxLevels, minLeafLevels, maxChildren, weightRange, level=0, usedLabels = set()):
    while (True):
        node_label = "".join([random.choice(labels) for i in range(4)])
        if not node_label in usedLabels:
            break

    usedLabels.add(node_label)

    numChildren = random.randint(1,maxChildren) if level < minLeafLevels else (random.randint(0,maxChildren) if level < maxLevels else 0)

    if numChildren > 0:
        childNodes = [randomTree(maxLevels, minLeafLevels, maxChildren, weightRange, level+1, usedLabels) for i in range(numChildren)]
        weights = [int(random.uniform(*weightRange)) if random.random() < 0.25 else round(random.uniform(*weightRange), random.randint(0,4)) for i in range(numChildren)]
        return PhyloTree_sol.PhyloTree(node_label, children = list(zip(childNodes, weights)))
    
    return PhyloTree_sol.PhyloTree(label=node_label)




def addAltLabel(r):
    """Add an alternate label to the tree, useful for testing"""
    if r.isLeaf():
        if type(r.label)==str and r.label:
            r.altLabel = r.label
            return r.altLabel
        else:
            raise ValueError("Bad leaf label")
    else:
        r.altLabel = "|".join(sorted([addAltLabel(t) for t,w in r.children]))
        return r.altLabel

def sameTrees(T1, T2):
    """Determined if the trees are the same.
    Assumes every non-root node has a unique label."""
    if T1.label != T2.label:
        return False;
    if T1.numChildren() != T2.numChildren():
        return False
    if T1.isLeaf():  
        return True

    C1 = sorted(T1.children, key = lambda x: x[0].label)
    C2 = sorted(T2.children, key = lambda x: x[0].label)
    return all([sameTrees(c1,c2) and abs(float(w1) - float(w2)) < 0.0001 for (c1,w1),(c2,w2) in zip(C1,C2)])


def translateTree(r):
    # Convert a solution tree to a student tree
    return PhyloTree.PhyloTree(label = r.label, children = [(translateTree(t),w) for t,w in r.children])

def createTrees():
    t = randomTree(5, 4, 4, [5,15])
    t2 = translateTree(t)
    return t,t2

def tree2set(r):
    return {r} | {x for t,w in r.children for x in tree2set(t)}

def karroLeafSet(tree):
    return {tree} if len(tree.children)==0 else set([x for t,w in tree.children for x in karroLeafSet(t)])


###########################
# Specific test trees
A = PhyloTree_sol.PhyloTree("A")                         # A;
B = PhyloTree_sol.PhyloTree("B")                         # B;
C = PhyloTree_sol.PhyloTree("C", [(A,1),(B,1)])          # (A:1,B:1):C;
D = PhyloTree_sol.PhyloTree("D")                         # D;
E = PhyloTree_sol.PhyloTree("E", [(D,1.5)])              # (D:1.5)E;
F = PhyloTree_sol.PhyloTree("F", [(E,2), (C,3)])         # ((D:1.5)E:2,(A:1,B:1)C:3)F;
G = PhyloTree_sol.PhyloTree("G", [(A,10),(B,20),(E,30)]) # (A:10,B:20,(D:1.5)E:30)G;
U = PhyloTree_sol.PhyloTree(children = [(A,27),(B,92)])  # (A:27,B:92);
TREELIST = [A,B,C,D,E,F,G,U]
TLIST = [translateTree(x) for x in TREELIST]

#############################
def test1(t = None):
    """Test len()"""
    t,t2 = (t, translateTree(t)) if t else createTrees()
    return len(t) == len(t2)

def test2(t = None):
    """Test str()"""
    modifier = 1

    t,t2 = (t, translateTree(t)) if t else createTrees()
    s = str(t2)
    if s[-1] != ';':
        modifier = 0.95
        s += ';'

    t3 = parseNewick_sol.parseNewick(s)
    return modifier*float(sameTrees(t,t3))

def test3(t = None):
    """Test copy"""
    t,t2 = (t, translateTree(t)) if t else createTrees()
    t3 = t2.copy()

    S2 = tree2set(t2)
    S3 = tree2set(t3)
    return sameTrees(t2,t2) and not (S2 & S3)

def test4(t = None):
    """Test leafSet"""
    t,t2 = (t, translateTree(t)) if t else createTrees()
    if t2.leafSet() == karroLeafSet(t2):
        return 1
    elif set(t2.leafSet()) == set(karroLeafSet(t2)):
        return 0.95
    return 0

def test5(t = None):
    """Test treeWeight"""
    t,t2 = (t, translateTree(t)) if t else createTrees()
    return abs(t.treeWeight() - t2.treeWeight()) < 0.0001

def test6(t = None):
    """Test mapDepth"""
    t,t2 = (t, translateTree(t)) if t else createTrees()
    D1 = t.mapDepth()
    D2 = t2.mapDepth()
    return set(D1.keys()) == set(D2.keys()) and all([abs(D1[x]-D2[x]) < 0.00001 for x in D1.keys()])

def test7a(t = None):
    """Test mapPath"""
    t,t2 = (t, translateTree(t)) if t else createTrees()
    D1 = t.mapPath()
    try:
        D2 = t2.mapPath({}, [])
    except:
        try:
            D2 = t2.mapPath([], {})
        except:
            D2 = t2.mapPath()

    S1 = set(D1.keys())
    S2 = set(D2.keys())

    score = 0.2 if S1 == S2 else (0.1 if S1 <= S2 else 0) 

    correct = 0.0
    for x in S1:
        if [y.label for y in D1[x]] != [y.label for y in D2[x]]:
            continue
        correct += 1
    score += 0.8*(correct / len(S1))
    
    return score

def test7b(t = None):
    """Test mapPath if the student fogot to add the root to the path"""
    t,t2 = (t, translateTree(t)) if t else createTrees()
    D1 = t.mapPath()
    try:
        D2 = t2.mapPath({}, [])
    except:
        try:
            D2 = t2.mapPath([], {})
        except:
            D2 = t2.mapPath()

    S1 = set(D1.keys())
    S2 = set(D2.keys())

    score = 0.2 if S1 == S2 else (0.1 if S1 <= S2 else 0) 

    correct = 0.0
    for x in S1:
        if [y.label for y in D1[x][:-1]] != [y.label for y in D2[x]]:
            continue
        correct += 1
    score += 0.8*(correct / len(S1))
    
    return 0.95*score

def test7(t = None):
    return max(test7a(t), test7b(t))

def test8a():
    return translateTree(C).MRCA('A','B').label == 'C'

def test8b():
    return translateTree(F).MRCA('A','B').label == 'C'

def test8c():
    return translateTree(F).MRCA('A', 'D').label =='F'


def test8d():
    """Test MRCA"""
    while True:
        t,t2 = createTrees()

        S = karroLeafSet(t)
        if len(S) >= 2:
            break
        
    while True:
        n1 = random.choice(list(S))
        n2 = random.choice(list(S))
        if n1 != n2:
            break

    return (t.MRCA(n1.label,n2.label).label == t2.MRCA(n1.label,n2.label).label)
    

def test8():
    return 0.5*(float(test8a()) + float(test8b()) + float(test8c()))/3 + 0.5*test8d()

def rootByOutgroupTest(t = None, outgroup = None, frac = 0.5, root_label = ""):
    t,t2 = ((t,translateTree(t)) if t else createTrees())

    if not outgroup:
        outgroup = random.choice(list(t.leafSet())).label

    x = t.rootByOutgroup(outgroup, frac, root_label)   # Correct solution
    y1 = t2.rootByOutgroup(outgroup, frac, root_label)   # Student solution if they used frac correctly
    y2 = t2.rootByOutgroup(outgroup, 1-frac, root_label)   # Student solution if they got frac backwards (no penalty)
    return sameTrees(x,y1) or sameTrees(x,y2)

def test9():
    return rootByOutgroupTest(t = PhyloTree_sol.PhyloTree('A', [ (PhyloTree_sol.PhyloTree('B'), 5), (PhyloTree_sol.PhyloTree('C'), 10)  ]), outgroup = 'B')


def test10(t = None):
    return rootByOutgroupTest(t)

def test11(t = None):
    return rootByOutgroupTest(t, root_label = "ABCDEF")
    
def test12(t = None):
    return rootByOutgroupTest(t, frac = 0.333)
    
def test(i, t = None):
    try:
        return float(eval("test" + str(i))())
    except:
        return 0

def avgTest(i, n):
    L = [test(i) for j in range(n)]
    return sum(L)/float(len(L))

def mean(L):
    return float(sum(L))/len(L)

def run_test(n, show_results = True, RNGseed = None):
    if RNGseed:
        random.seed(RNGseed)

    P = [avgTest(i,n) for i in range(1,9)]
    # Tests for variations of rootByOutgroup
    t1 = avgTest(9, n)     # Test last method with a tree of height 2 
    t2 = avgTest(10, n)    # Test last method with frac=0.5, root_label = ""
    t3 = avgTest(11, n)    # Tast last method with frac = 0.5, root_label = "ABCDEF"
    t4 = avgTest(12, n)    # Test last method with frac = 0.333, root_label = ""
    P.append(0.2*t1 + 0.6*(max(t2,t3,t4)) + 0.2*((sum((t2,t3,t4)) - max((t2,t3,t4)))/2))    # 20% for t1, 60% for max(t1,t2,t3), 20% for the other two

    if show_results:
        print("\n".join(["Test %d: %5.3f" % (i+1, v) for i,v in enumerate(P)]) + "\n")

    project_score = mean(P)
    print("SCORE: %5.2f"  % (100*project_score))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Testing dnaSeq project')
    parser.add_argument('-n', type = int, help = "Each result will be averaged over n runs", default = 100)
    parser.add_argument('-s', type = int, help = "Seed for random number generator", default = None)
    parser.add_argument('-d', '--display', action = "store_true", help = "display indvidual test results", default = False)
    args = parser.parse_args()
    run_test(args.n, RNGseed = args.s, show_results = args.display)


