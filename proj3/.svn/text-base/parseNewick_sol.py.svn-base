import re
import PhyloTree_sol

# Grammer:
#   TREE     --> SUBTREE;
#   SUBTREE  --> LEAF 5| (INTERNAL) | (INTERNAL)LABEL
#   INTERNAL --> SUBTREE:WEIGHT | SUBTREE:WEIGHT,INTERNAL
#   LEAF     --> LABEL
#   LABEL    --> <string of 1 or more chaacters>
#   WEIGHT   --> <int> | <float>
def parseNewick(s):
    if not s[-1] == ';':
        raise NewickParseError("Newick string missing terminating semi-colon")
    return _parseSubtree(s[:-1])
    

r = re.compile("(\(.*\))?(\w+)?")
subtreeRE = re.compile("^(\((.*)\))?(\w+)?")
def _parseSubtree(s, used = set()):
    if not s.strip():
        raise NewickParseError("Recursed to empty substring.")

    internal_str, label = subtreeRE.search(s).group(2,3)
    if not internal_str:
        return PhyloTree_sol.PhyloTree(label)

    children = []    # Hold the substrees
    level = 0
    last = 0
    for i,c in enumerate(internal_str + ","):
        if c == '(':
            level += 1
        elif c == ')':
            level -= 1
        elif c == ',' and level == 0:
            subtree = internal_str[last:i].strip()
            split = subtree.rfind(":")
            if split == -1:
                raise NewickParseError("Bad INTERNAL: ", subtree)
            children.append( (_parseSubtree(subtree[:split].strip(), used), float(subtree[split+1:].strip())) )
            last = i+1
    return PhyloTree_sol.PhyloTree(label, children)
    

