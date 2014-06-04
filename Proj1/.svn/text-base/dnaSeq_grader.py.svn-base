import sys
sys.path.append(".")   # Necessary for use in an out-of-date emacs IDE\
import dnaSeq_sol
import dnaSeq
import random
import argparse


def checkType(o, t):
    return type(o) == t if python_version >= 3 else isinstance(o, t)


def randomSeq(n, s = "ACGT"):
    return "".join([random.choice(s) for i in range(n)])

#########################################
# Constructor
def test1():
    """Test that the constructor run on a legit sequence withot exceptions"""
    s = dnaSeq.dnaSeq("ACGTacgtNn")
    return isinstance(s, dnaSeq.dnaSeq)
    
def test2():
    """Testing constructor exception"""
    try:
        dnaSeq.dnaSeq(2)
    except Exception as e:
        return isinstance(e, ValueError)
    return False

def test3():
    """Testing constructor exception"""
    try:
        dnaSeq.dnaSeq("ABCD")
    except Exception as e:
        return type(e) == dnaSeq.DNAError

#########################################
# Test len()
def test4():
    n = random.randint(50,150)
    s = dnaSeq.dnaSeq(randomSeq(n))
    return len(s) == n

#########################################
# Test str()
def test5():
    n = random.randint(50,400)
    s = randomSeq(n)
    return str(dnaSeq.dnaSeq(s)) == s

#########################################
# Test getitem -- index
def test6():
    n = random.randint(50,400)
    s = randomSeq(n)
    seq = dnaSeq.dnaSeq(s)
    i = random.randint(0,n-1)
    j = -1*random.randint(1,n)
    return 0.5*int(seq[i] == s[i]) + 0.5*int(seq[j]==seq[j])
    
#########################################
# Test getitem -- slice
def test7():
    n = random.randint(50,400)
    s = randomSeq(n)
    seq = dnaSeq.dnaSeq(s)
    i = random.randint(0,n//2)
    j = random.randint(n//2+1,n)

    i2 = -1*random.randint(0,n//2)
    j2 = -1*random.randint(n//2+1,n)

    i3 = random.randint(0,n//2)
    j3 = random.randint(n//2+1,n)

    return (1/3.0)*int(str(seq[i:j])==s[i:j]) + (1/3.0)*int(str(seq[i2:j2])==s[i2:j2]) + (1/3.0)*int(str(seq[i3:j3:2])==s[i3:j3:2]) 

def test8():    # Make sure index slice returns a sequence
    s = dnaSeq.dnaSeq("ACGT")[0:2]
    return isinstance(s, dnaSeq.dnaSeq)
    
#########################################
# Test setitem
def test9():
    n = random.randint(50,400)
    s = randomSeq(n)
    seq = dnaSeq.dnaSeq(s)

    i = random.randint(10,n-10)
    s = s[:i] + "a" + s[i+1:]
    seq[i] = "a"

    return str(seq) == s

def test10():
    n = random.randint(50,400)
    s = randomSeq(n)
    seq = dnaSeq.dnaSeq(s)

    i = random.randint(10,n-10)
    j = random.randint(10,n-10)
    if i > j:
        i,j = j,i

    s = s[:i] + "a" + s[j:]
    seq[i:j] = 'a'

    n2 = random.randint(50,400)
    s2 = randomSeq(n)
    seq2 = dnaSeq.dnaSeq(s2)

    i2 = random.randint(10,n2-10)
    j2 = random.randint(10,n2-10)
    if i2 > j2:
        i2,j2 = j2,i2

    s2 = s2[:i2] + "aaa" + s2[j2:]
    seq2[i2:j2] = 'aaa'

    return 0.5*int(str(seq)==s) + 0.5*int(str(seq2)==s2)

#################################
# Test __add__
def test11():
    n = random.randint(50,400)
    s = randomSeq(n)
    seq = dnaSeq.dnaSeq(s)

    n2 = random.randint(50,400)
    s2 = randomSeq(n2)
    seq2 = dnaSeq.dnaSeq(s2)

    return 0.8*int(str(seq+seq2) == s+s2) + 0.2*int(isinstance(seq+seq2, dnaSeq.dnaSeq))

#################################
# Test __contains__
def test12():
    n = random.randint(50,400)
    s = randomSeq(n)
    seq = dnaSeq.dnaSeq(s)

    i = random.randint(0,n)
    j = random.randint(0,n)
    if i > j:
        i,j = j,i

    return str(dnaSeq.dnaSeq(s[i:j])) in seq

def test13():
    n = random.randint(50,400)
    s = randomSeq(n)
    seq = dnaSeq.dnaSeq(s)

    i = random.randint(0,n)
    j = random.randint(0,n)
    if i > j:
        i,j = j,i

    return dnaSeq.dnaSeq(s[i:j]) in seq


################################
# Test complement / reverse / reverse-complement
def test14():
    n = random.randint(50,400)
    s = randomSeq(n)
    seq1 = dnaSeq.dnaSeq(s)   
    seq2 = dnaSeq_sol.dnaSeq(s)

    seq1.complement()
    seq2.complement()
    return str(seq1)==str(seq2)

def test15():
    n = random.randint(50,400)
    s = randomSeq(n)
    seq1 = dnaSeq.dnaSeq(s)   
    seq2 = dnaSeq_sol.dnaSeq(s)

    seq1.reverse()
    seq2.reverse()
    return str(seq1)==str(seq2)

def test16():
    n = random.randint(50,400)
    s = randomSeq(n)
    seq1 = dnaSeq.dnaSeq(s)   
    seq2 = dnaSeq_sol.dnaSeq(s)

    seq1.reverse_complement()
    seq2.reverse_complement()
    return str(seq1)==str(seq2)

############################
# Test readFA 
def test17():
    L = [dnaSeq_sol.dnaSeq(randomSeq(random.randint(300,600))) for i in range(5)]
    for d in L:
        d.info = ">" + "".join([random.choice("ABCDEFH ") for i in range(30)])
    dnaSeq_sol.writeFA(L, "test.fa")
    
    L2 = dnaSeq.readFA("test.fa")
    return all([o1.info.lstrip('>').strip()==o2.info.lstrip('>').strip() and str(o1)==str(o2) for o1,o2 in zip(L,L2)])


############################
# Test writeFA
def test18():
    """Test to make sure it creates a file that is readable"""
    L = [dnaSeq.dnaSeq(randomSeq(random.randint(300,600))) for i in range(5)]
    for d in L:
        d.info = ">" + "".join([random.choice("ABCDEFH ") for i in range(30)])
    dnaSeq.writeFA(L, "test1.fa")
    L2 = dnaSeq.readFA("test1.fa")
    return all([o1.info.lstrip('>').strip()==o2.info.lstrip('>').strip() and str(o1)==str(o2) for o1,o2 in zip(L,L2)])

def test19():
    """Test with differnt column width"""
    colwidth = random.randint(60,100)
    L = [dnaSeq_sol.dnaSeq(randomSeq(random.randint(300,600))) for i in range(5)]
    for d in L:
        d.info = ">" + "".join([random.choice("ABCDEFH ") for i in range(30)])
    dnaSeq_sol.writeFA(L, "test1.fa", col_width = colwidth)

    L2 = []
    for l in L:
        d = dnaSeq.dnaSeq(str(l))
        d.info = l.info
        L2.append(d)

    dnaSeq.writeFA(L2, "test2.fa", col_width = colwidth)

    ## compare files
    for l1,l2 in zip(open("test1.fa"), open("test2.fa")):
        if l1[0] == '>':
            if l1.lstrip('>').strip() != l2.lstrip('>').strip():
                return False
        else:
            if l1 != l2:
                return False

    return True

def test(i):
    try:
        return float(eval("test" + str(i))())
    except:
        return 0
    
def avgTest(i, n):
    L = [test(i) for j in range(n)]
    return sum(L)/len(L)

def mean(L):
    return float(sum(L))/len(L)

def run_test(n, show_results = True, RNGseed = None):
    if RNGseed:
        random.seed(RNGseed)

    tests = [None] + [avgTest(i,n) for i in range(1,20)]

    P = [
         0.7*tests[1] + 0.15*tests[2] + 0.15*tests[3],   # P0: Consturctor test
         tests[4],   # P1: Length test
         tests[5],   # P2: str() test
         0.5*tests[6] + 0.4*tests[7] + 0.1*tests[8], # P3: getitem test
         0.5*tests[9] + 0.5*tests[10], # P4: setitem test
         tests[11], # P5: __add__ test
         0.8*max(tests[12],tests[13]) + 0.2*min(tests[12], tests[13]), # P6: __contains__ test
         tests[14], # P7: complement test
         tests[15], # P8: reverse test
         tests[16], # P9: reverse complement test
         tests[17], # P10: readFA test
         0.7*tests[18] + 0.3*tests[19] # P11: witeFA test
        ]
    if show_results:
        print("\n".join(["P%d: %5.3f" % (i, v) for i,v in enumerate(P)]) + "\n")

    project_score = 100*(0.6*mean(P[:7]) + 0.1*mean(P[7:10]) + 0.3*mean(P[10:]))
    print("SCORE: %5.2f"  % (project_score))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Testing dnaSeq project')
    parser.add_argument('-n', type = int, help = "Each result will be averaged over n runs", default = 100)
    parser.add_argument('-s', type = int, help = "Seed for random number generator", default = None)
    parser.add_argument('-d', '--display', action = "store_true", help = "display indvidual test results", default = False)
    args = parser.parse_args()
    run_test(args.n, RNGseed = args.s, show_results = args.display)
