from re import split

class AlignError(Exception):
    def __init__(self, msg = ""):
        self.msg = str(msg)
    def __str__(self):
        return self.msg

def prettyPrintMatrix(M):
    """
    Pretty-print a matrix to make it easily human readable.
    M must be a list of size n, where each element is a list of m integers.
    """
    l = max([len(str(x)) for row in M for x in row])
    print("\n".join(["".join([("{:<%s}" % (l+2)).format(i) for i in row]) for row in M]))

def readScoringMatrix(scoreMatrixFile):
    """Read in a GenBank-formated scoring matrix and
    return a scoring matrix in dictionary form"""
    with open(scoreMatrixFile) as fp:
        line = fp.readline()
        while line[0] == '#':
            line = fp.readline()

        char_list = split("\s+", line.strip())
        D = {}

        for line in fp:
            arr = split("\s+", line.strip())
            D[arr[0]] = {char_list[i]:float(val)  for i,val in enumerate(arr[1:])}
    return D

def genScoringMatrix(match, mismatch):
    """Generate a DNA scoring matrix with the given mathc and mismatch penalties"""
    return {x:{y:match if x == y else mismatch for y in "ACGT"} for x in "ACGT"}

def scoreAlignment(aln1, aln2, S, gap):
    """Compute the score of an alignment relative to scoring matrix S and gap
    penelty gap (represented as a non-negative number)"""

    if {type(aln1),type(aln2)} != {str}:
        raise AlignError("Aligned sequences are not strings")
    if len(aln1) != len(aln2):
        raise AlignError("Aligned sequences are not of the same length")
    if not (set(aln1) | set(aln2) <= set(S.keys()) | {'-'}):
        raise AlignError("Aligned sequences contain bad characters")
    
    return sum([-gap if '-' in {x,y} else S[x][y]  for x,y in zip(aln1,aln2)])


def scoreAlignmentAffine(aln1, aln2, S, o, c):
    """Compute the score of an alignment relative to scoring matrix S, an open
    gap score of s, and a gap continuations score of c"""
    type = 0
    score = 0
    for a,b in zip(aln1,aln2):
        if a == "-":
            score -= (o if type != 1 else 0) + c
            type = 1
        elif b == "-":
            score -= (o if type != 2 else 0) + c
            type = 2
        else:
            score += S[a][b]
            type = 0
    return score

