import argparse
import alignment
import alignment_sol
import alignment_util
import random
import sys


labels = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

D1 = alignment_util.readScoringMatrix("DNA1.txt")
D2 = alignment_util.readScoringMatrix("DNA2.txt")
D3 = alignment_util.genScoringMatrix(10, -5)
D4 = alignment_util.genScoringMatrix(5,0)
B  = alignment_util.readScoringMatrix("Blosum62.txt")

#############################
def test_SW(seq1, seq2, S, g):
    s, a1, a2 = alignment.SmithWaterman(seq1, seq2, S, g)
    s_sol, a1_sol, a2_sol = alignment_sol.SmithWaterman(seq1, seq2, S, g)
    score = 0
    
    # First: test that the function has returned has an optimal alignment
    if alignment_util.scoreAlignment(a1, a2, S, g) == s_sol:
        score += 70

    # Second: test that the function has returned an optimal score
    if s == s_sol:
        score += 20

    # Third: Test that the function has returned the correct score for the alignment
    if alignment_util.scoreAlignment(a1, a2, S, g) == s:
        score += 10

    return score/100.0

def test_sym(seq1, seq2, S, g):
    """Test for symmetric alignments -- partial penalty if it doesn't work when parameters are switched."""
    s1 = test_SW(seq1, seq2, S, g)
    s2 = test_SW(seq2, seq1, S, g)
    return 0.8*max(s1,s2) + 0.2*min(s1,s2)


def test_SWA(seq1, seq2, S, o, c):
    s, a1, a2 = alignment.SmithWatermanAffine(seq1, seq2, S, o, c)
    s_sol, a1_sol, a2_sol = alignment_sol.SmithWatermanAffine(seq1, seq2, S, o, c)
    score = 0
    
    # First: test that the function has returned has an optimal alignment
    if alignment_util.scoreAlignmentAffine(a1, a2, S, o, c) == s_sol:
        score += 70

    # Second: test that the function has returned an optimal score
    if s == s_sol:
        score += 20

    # Third: Test that the function has returned the correct score for the alignment
    if alignment_util.scoreAlignmentAffine(a1, a2, S, o, c) == s:
        score += 10

    return score/100.0

def test_symA(seq1, seq2, S, o, c):
    """Test for symmetric alignments -- partial penalty if it doesn't work when parameters are switched."""
    s1 = test_SWA(seq1, seq2, S, o, c)
    s2 = test_SWA(seq2, seq1, S, o, c)
    return 0.8*max(s1,s2) + 0.2*min(s1,s2)

################################
# First set of tests: make sure this works "globally" -- on sequences
# where the correct local alignment is global.
def test1():
    """Make sure it works in a trivial case."""
    s1 = "CCC"
    s2 = "CCC"
    return test_sym(s1, s2, D3, 3)

def test2():
    """Make sure it works with a mis-match"""
    s1 = "TTT"
    s2 = "TGT"
    return test_sym(s1, s2, D3, 1)

def test3():
    """Test with gap"""
    s1 = "TTTATTT"
    s2 = "TTTTTT"
    return test_sym(s1, s2, D3, 2)

def test4():
    s1 = "TTTAGATTT"
    s2 = "TTTGTTT"
    return test_sym(s1, s2, D3, 2)

def test5():
    s1 = "TTTAAGTTT"
    s2 = "TTTGTTT"
    return test_sym(s1, s2, D3, 3)

#########################
# Second set of tests: make sure this works when correct alignment is
# a local alignment -- discaring at a portion of one of the sequences.
def test6():
    s1 = "GGGTTT"
    s2 = "AAAGGGTTT"
    return test_sym(s1, s2, D3, 1)

def test7():
    s1 = "AAAGGGTTT"
    s2 = "AAAGGG"
    return test_sym(s1, s2, D3, 1)

def test8():
    s1 = "GGG"
    s2 = "AAAGGGTTT"
    return test_sym(s1, s2, D3, 1)

def test9():
    s1 = "TTTGGGTTT"
    s2 = "TTTTTT"
    return test_sym(s1, s2, D3, 10)

def test10():
    s1 = "TTTTTGTTTTT"
    s2 = "CCCTTTTTTTTTTCCC"
    return test_sym(s1, s2, D3, 1)

#########################
# Smith-Waterman Affine: Make sure this works when
# the correct alignment is global.
def test11():
    s1 = "AAA"
    s2 = "ATA"
    return test_symA(s1, s2, D3, 5, 1)

def test12():
    s1 = "AAAAATTTAAAAA"
    s2 = "AAAAAAAAAA"
    return test_symA(s1, s2, D3, 5, 1)

def test13():
    s1 = "TTTTTCACTTTTT"
    s2 = "TTTTTGTTTTT"
    return test_symA(s1, s2, D2, 1, 1)

def test14():
    s1 = "TTTTTCACTTTTT"
    s2 = "TTTTTGTTTTT"
    return test_symA(s1, s2, D2, 2, 1)

def test15():
    s1 = "TTTTTCACTTTTTTTTTT"
    s2 = "TTTTTTTTTACATTTTTT"
    return test_symA(s1, s2, D3, 6, 2)

#########################
# Smith-Waterman Affine: Make sure this works when
# the correct alignment is local.
def test16():
    s1 = "TTTAAA"
    s2 = "AAA"
    return test_symA(s1, s2, D3, 5, 1)

def test17():
    s1 = "AAATTT"
    s2 = "AAA"
    return test_symA(s1, s2, D3, 5, 1)

def test18():
    s1 = "TTTAAATTT"
    s2 = "TTTTTT"
    return test_symA(s1, s2, D3, 10, 10)

def test19():
    s1 = "TTTAAATTT"
    s2 = "CCCGGGCCC"
    return test_symA(s1, s2, D3, 2, 1)

def test20():
    s1 = "CCCTTTTTCACTTTTTTTT"
    s2 = "GGGTTTTTGTTTTTAAA"
    return test_symA(s1, s2, D2, 2, 1)



#########################
# Smith-Waterman Affine: Make sure its
# not just a greedy matching
def test21():
    s1 = "AAATTTCAAA"
    s2 = "AAACAAA"
    return test_symA(s1, s2, D4, 5, 1)

def test22():
    s1 = "AAATTTCAAA"
    s2 = "AAACAA"
    return test_symA(s1, s2, D4, 5, 1)

def test23():
    s1 = "TTTAAATTTCAACCC"
    s2 = "GGGAAACAAA"
    return test_symA(s1, s2, D4, 5, 1)

#####################
# Test each method on a longer alignment
def test24():
    s1 = "A"*100 + "CCC" + "G"*100
    s2 = "A"*95 + "G"*95

    t1 = test_sym(s1, s2, D3, 1)
    t2 = test_symA(s1, s1, D1, 1, 1)   # Yes -- this is on purpose
    t3 = test_symA(s1, s2, D3, 3, 1)
    return (t1 + t2 + t3)/3.0

####################

def test(i):
    try:
        return float(eval("test" + str(i))())
    except:
        return 0

def avgTest(R):
    L = [test(i) for i in R]
    return sum(L)/float(len(L))


def run_test(show_results = True):
    t1 = avgTest(range(1,5))
    t2 = avgTest(range(6,10))
    t3 = avgTest(range(11,15))
    t4 = avgTest(range(16,20))
    t5 = avgTest(range(21,23))
    t6 = test24()
    if show_results:
        print("Test of SW on gloablly alignable sequences (10%s): %5.2f/100" % ("%", 100*t1))
        print("Test of SW on locally alignable sequences (30%s): %5.2f/100" % ("%", 100*t2))
        print("Test of SW-affine on gloablly alignable sequences (10%s): %5.2f/100" % ("%", 100*t3))
        print("Test of SW-affine on locally alignable sequences (30%s): %5.2f/100" % ("%", 100*t4))        
        print("Test of SW-affine for use of greedy matching (10%s): %5.2f/100" % ("%", 100*t5))        
        print("Test of both algorithms on a longer sequence (10%s): %5.2f/100" % ("%", 100*t6))
    score = 0.1*t1 + 0.3*t2 + 0.1*t3 + 0.3*t4 + 0.1*t5 + 0.1*t6
    print("SCORE: %5.2f" % (100*score))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Testing dnaSeq project')
    parser.add_argument('-d', '--display', action = "store_true", help = "display indvidual test results", default = False)
    args = parser.parse_args()
    run_test(show_results = args.display)


