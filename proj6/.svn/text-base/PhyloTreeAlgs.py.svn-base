import sys
sys.path.append("../lib/")
sys.path.append(".")

from PhyloTree import *
from PhyloTree_util import *
import random

def UPGMA(DIST):
    """
    UPGMA Algorithm:
    Input:
    * DIST: A dictinary mapping label pairs to distnce (e.g. for x,y in leafLabels: D[x][y] 
            is the tree dist. from x to y, where x and y are node labels.
    Assumptions: 
    * For every x,y in D.keys(), D[x][y] is defind.
    * DIST[x][x] == 0 and DIST[x][y] == DIST[y][x].
    * DIST must describe a tree with the ultrametric property.
    Return: A PhloyTree object T representing a full, unrooted, binary, ultrametric tree with leaf 
            lables equal to D.keys(), such that the path distance from leaf x to leaf y is D[x][y] for 
            any x,y.
    """
    leafLabels = DIST.keys()
    S = {PhyloTree(label) for label in leafLabels}
    TreeSize = {v:1 for v in S}
    TreeHeight = {v:0 for v in S}
    D = {x:{y:DIST[x.label][y.label] for y in S if y != x} for x in S}

    while len(S) > 1:
        # Find closest nodes
        i = iter(S)
        min_x, min_y = next(i), next(i)
        for x in S:
            for y in S:
                if x != y:
                    if D[x][y] < D[min_x][min_y]:
                        min_x, min_y = x, y

        # Create new node
        new_height = D[min_x][min_y]/2
        C = [(min_x, new_height - TreeHeight[min_x]), (min_y, new_height - TreeHeight[min_y])]
        z = PhyloTree(label = x.label + y.label, children = C)
        TreeHeight[z] = new_height

        # Augment D to reflect distance from z to other nodes
        D[z] = {}
        S.remove(min_x)
        S.remove(min_y)
        for a in S:
            new_dist = (D[min_x][a] + D[min_y][a])/2
            D[a][z] = new_dist
            D[z][a] = new_dist
            
        S.add(z)

    return S.pop() 
        
    
def NJ(DIST):
    """
    UPGMA Algorithm:
    Input:
    * DIST: A dictinary mapping label pairs to distnce (e.g. for x,y in leafLabels: D[x][y] 
            is the tree dist. from x to y, where x and y are node labels.
    Assumptions:
    * For every x,y in leafLabels, D[x][y] is definied.
    * DIST[x][x] == 0 and DIST[x][y] == DIST[y][x].
    Return: A PhloyTree object T representing a binary, unrooted tree with leaf labels equal to
            D.keys(), such that the path distance from leaf x to leaf y is D[x][y] for for
            any x,y.
    Note: While the result is an unrooted binary tree, this will be impossible to exactly represent
          with a PhyloTree object designed to hold a root tree.  Your reslting tree will presumably
          have a root note with three children, and the root node will not necessarily correspond to 
          the actual (unknown) root.
    """
    leafLabels = DIST.keys()
    S = {PhyloTree(label=l) for l in leafLabels}        # Leaf node set
    D = {v:{u:DIST[v.label][u.label] for u in S} for v in S}

    label_num = 1
    while len(S) > 3:
        R = {v:sum([D[v][u] for u in S])/(len(S)-2) for v in S}
        u,v = min([(D[u][v] - (R[u] + R[v]),u,v) for u in S for v in S if u!=v])[1:]  # nodes to merge
        d_uw = 0.5*(D[u][v] + R[u] - R[v])    # Distance from u to new node w
        d_vw = 0.5*(D[v][u] + R[v] - R[u])    # Distance from u to new node w

        S.remove(u)
        S.remove(v)

        w = PhyloTree(label = "L" + str(label_num), children = [(u,d_uw),(v,d_vw)])
        label_num += 1
        D[w] = {w:0}
        for z in S:
            D[w][z] = 0.5*(D[u][z] + D[v][z] - D[u][v])
            D[z][w] = D[w][z]

        S.add(w)
        
    u,v,w = list(S)
    d_ut = 0.5*(D[u][v] + D[u][w] - D[v][w])
    d_vt = 0.5*(D[v][u] + D[v][w] - D[u][w])
    d_wt = 0.5*(D[w][u] + D[w][v] - D[u][v])

    return PhyloTree(children=[(u,d_ut),(v,d_vt),(w,d_wt)])



