import sys
import unittest
import re
import alignment
import alignment_util


def basicTest(seq1, seq2, S, gap, correct):
    s, a1, a2 = SmithWaterman(seq1, seq2, S, gap)
    return s == scoreAlignment(a1, a2, S, gap) and s == correct


def fullSuite():
    suite = unittest.TestSuite()
    suite.addTest(SmithWatermanGlobalTester('global_test1'))
    suite.addTest(SmithWatermanGlobalTester('global_test2'))
    suite.addTest(SmithWatermanGlobalTester('global_test3'))
    suite.addTest(SmithWatermanGlobalTester('global_test4'))
    suite.addTest(SmithWatermanGlobalTester('global_test5'))
    suite.addTest(SmithWatermanGlobalTester('global_test6'))
    suite.addTest(SmithWatermanGlobalTester('global_test7'))
    suite.addTest(SmithWatermanLocalTester('local_test1'))
    suite.addTest(SmithWatermanLocalTester('local_test2'))
    suite.addTest(SmithWatermanLocalTester('local_test3'))
    suite.addTest(SmithWatermanLocalTester('local_test4'))
    suite.addTest(SmithWatermanLocalTester('local_test5'))
    suite.addTest(SmithWatermanLocalTester('local_test6'))
    suite.addTest(SmithWatermanLocalTester('local_test7'))
    suite.addTest(SmithWatermanAffineTester('test1'))
    suite.addTest(SmithWatermanAffineTester('test2'))
    suite.addTest(SmithWatermanAffineTester('test3'))
    suite.addTest(SmithWatermanAffineTester('test4'))
    suite.addTest(SmithWatermanAffineTester('test5'))
    suite.addTest(SmithWatermanAffineTester('test6'))
    suite.addTest(SmithWatermanAffineTester('test7'))
    suite.addTest(SmithWatermanAffineTester('mcalliTest'))
    return suite


class BasicAlignmentTester(unittest.TestCase):
    def setUp(self):
        self.S1 = alignment_util.readScoringMatrix("DNA1.txt")
        self.S2 = alignment_util.readScoringMatrix("DNA2.txt")
        self.S3 = alignment_util.genScoringMatrix(10,-10)
        self.S4 = alignment_util.readScoringMatrix("Blosum62.txt")

class SmithWatermanGlobalTester(BasicAlignmentTester):
    """In the following tests, the local alignments are the global alignments.  
    These should work if you didn't break anything with your modifications, 
    but do not test that you are computing optimal local alignments."""

    def global_test1(self):
        s1 = "AAA"
        s2 = "AAA"
        S = self.S1
        gap = 1
        correct_score = 3
        s, a1, a2 = alignment.SmithWaterman(s1, s2, S, gap)
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), s, "Returned alignment score does not match actual alignment score")
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), correct_score, "Score of returned alignment is not optimal")


    def global_test2(self):
        s1 = "AAA"
        s2 = "AAA"
        S = self.S2
        gap = 1
        correct_score = 12
        s, a1, a2 = alignment.SmithWaterman(s1, s2, S, gap)
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), s, "Returned alignment score does not match actual alignment score")
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), correct_score, "Score of returned alignment is not optimal")

    def global_test3(self):
        s1 = "AAAA"
        s2 = "AAGA"
        S = self.S2
        gap = 1
        correct_score = 13
        s, a1, a2 = alignment.SmithWaterman(s1, s2, S, gap)
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), s, "Returned alignment score does not match actual alignment score")
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), correct_score, "Score of returned alignment is not optimal")

        s, a1, a2 = alignment.SmithWaterman(s2, s1, S, gap)
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), s, "Returned alignment score does not match actual alignment score")
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), correct_score, "Score of returned alignment is not optimal")

    def global_test4(self):
        s1 = "AATT"
        s2 = "AACTT"
        S = self.S2
        gap = 1
        correct_score = 15
        s, a1, a2 = alignment.SmithWaterman(s1, s2, S, gap)
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), s, "Returned alignment score does not match actual alignment score")
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), correct_score, "Score of returned alignment is not optimal")

        s, a1, a2 = alignment.SmithWaterman(s2, s1, S, gap)
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), s, "Returned alignment score does not match actual alignment score")
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), correct_score, "Score of returned alignment is not optimal")

    def global_test5(self):
        s1 = "AATT"
        s2 = "AACTT"
        S = self.S2
        gap = 2
        correct_score = 14
        s, a1, a2 = alignment.SmithWaterman(s1, s2, S, gap)
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), s, "Returned alignment score does not match actual alignment score")
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), correct_score, "Score of returned alignment is not optimal")

        s, a1, a2 = alignment.SmithWaterman(s2, s1, S, gap)
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), s, "Returned alignment score does not match actual alignment score")
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), correct_score, "Score of returned alignment is not optimal")

    def global_test6(self):
        s1 = "AT"
        s2 = "ACCCT"
        S = self.S3
        gap = 1
        correct_score = 17
        s, a1, a2 = alignment.SmithWaterman(s1, s2, S, gap)
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), s, "Returned alignment score does not match actual alignment score")
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), correct_score, "Score of returned alignment is not optimal")

        s, a1, a2 = alignment.SmithWaterman(s2, s1, S, gap)
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), s, "Returned alignment score does not match actual alignment score")
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), correct_score, "Score of returned alignment is not optimal")

    def global_test7(self):
        s1 = 'GDQSHITZBZYHBPEMKASGGZXACAYVHVNSLYXLFSNZNBIZRRANXYVNKBPAZVGEVMBACMTPQSDSFMFQMBKVCXDVYQBLACIN'
        s2 = 'GDQSHITZBZYHMKASGGZXACAYVHVNBZYCCSLYXLFSNZNBIZRRANXYVNFBPAZVGEVMBACMTPQSDSFMFQMBKVCFZBXDVYQBLACIN'
        S = self.S4
        gap = 5
        correct_score = 375
        s, a1, a2 = alignment.SmithWaterman(s1, s2, S, gap)
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), s, "Returned alignment score does not match actual alignment score")
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), correct_score, "Score of returned alignment is not optimal")

        s, a1, a2 = alignment.SmithWaterman(s2, s1, S, gap)
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), s, "Returned alignment score does not match actual alignment score")
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), correct_score, "Score of returned alignment is not optimal")


class SmithWatermanLocalTester(BasicAlignmentTester):
    def local_test1(self):
        s1 = "CCCAAA"
        s2 = "TTTCCCTTT"
        S = self.S2
        gap = 10
        correct_score = 12
        s, a1, a2 = alignment.SmithWaterman(s1, s2, S, gap)
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), s, "Returned alignment score does not match actual alignment score")
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), correct_score, "Score of returned alignment is not optimal")

        s, a1, a2 = alignment.SmithWaterman(s2, s1, S, gap)
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), s, "Returned alignment score does not match actual alignment score")
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), correct_score, "Score of returned alignment is not optimal")

    def local_test2(self):
        s1 = "TTTCCC"
        s2 = "AAACCCTTT"
        S = self.S2
        gap = 10
        correct_score = 12
        s, a1, a2 = alignment.SmithWaterman(s1, s2, S, gap)
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), s, "Returned alignment score does not match actual alignment score")
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), correct_score, "Score of returned alignment is not optimal")

        s, a1, a2 = alignment.SmithWaterman(s2, s1, S, gap)
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), s, "Returned alignment score does not match actual alignment score")
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), correct_score, "Score of returned alignment is not optimal")


    def local_test3(self):
        s1 = "AAACCC"
        s2 = "TTTCCCTTT"
        S = self.S2
        gap = 10
        correct_score = 12
        s, a1, a2 = alignment.SmithWaterman(s1, s2, S, gap)
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), s, "Returned alignment score does not match actual alignment score")
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), correct_score, "Score of returned alignment is not optimal")

    def local_test4(self):
        s1 = "TTTCCCTTT"
        s2 = "AAACCCAAA"
        S = self.S2
        gap = 10
        correct_score = 12
        s, a1, a2 = alignment.SmithWaterman(s1, s2, S, gap)
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), s, "Returned alignment score does not match actual alignment score")
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), correct_score, "Score of returned alignment is not optimal")

    def local_test5(self):
        s1 = "AAACCCCCACCCCCAAA"
        s2 = "TTTCCCCCCCCCCTTT"
        S = self.S2
        gap = 4
        correct_score = 36
        s, a1, a2 = alignment.SmithWaterman(s1, s2, S, gap)
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), s, "Returned alignment score does not match actual alignment score")
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), correct_score, "Score of returned alignment is not optimal")

        s, a1, a2 = alignment.SmithWaterman(s2, s1, S, gap)
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), s, "Returned alignment score does not match actual alignment score")
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), correct_score, "Score of returned alignment is not optimal")


    def local_test6(self):
        s1 = "CCCCCACCCCC"
        s2 = "TTTCCCCCCCCCCTTT"
        S = self.S2
        gap = 4
        correct_score = 36
        s, a1, a2 = alignment.SmithWaterman(s1, s2, S, gap)
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), s, "Returned alignment score does not match actual alignment score")
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), correct_score, "Score of returned alignment is not optimal")

        s, a1, a2 = alignment.SmithWaterman(s2, s1, S, gap)
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), s, "Returned alignment score does not match actual alignment score")
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), correct_score, "Score of returned alignment is not optimal")

    def local_test7(self):
        s1 = 'WCFFEXTSYRKFNTHPEQRNGDQSHITZBZYHBPEMKASGGZXACAYVHVNSLYXLFSNZNBIZRRANXYVNKBPAZVGEVMBACMTPQSDSFMFQMBKVCXDVYQBLACINNGLTTSYDSYFZTFYBR'
        s2 = 'CPDAZPXGMFLHMZZKKBMTPCDTXPLWKSEAGDQSHITZBZYHMKASGGZXACAYVHVNBZYCCSLYXLFSNZNBIZRRANXYVNFBPAZVGEVMBACMTPQSDSFMFQMBKVCFZBXDVYQBLACINKXDPXMYGY'
        S = self.S4
        gap = 12
        correct_score = 298
        s, a1, a2 = alignment.SmithWaterman(s1, s2, S, gap)
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), s, "Returned alignment score does not match actual alignment score")
        self.assertEqual(alignment_util.scoreAlignment(a1,a2,S,gap), correct_score, "Score of returned alignment is not optimal")

class SmithWatermanAffineTester(BasicAlignmentTester):
    def test1(self):
        s1 = "AATTA"
        s2 = "AAA"
        S = self.S3
        o = 4
        c = 1
        correct_score = 24
        s, a1, a2 = alignment.SmithWatermanAffine(s1, s2, S, o, c)
        self.assertEqual(alignment_util.scoreAlignmentAffine(a1,a2,S,o,c), s, "Returned alignment score does not match actual alignment score")
        self.assertEqual(alignment_util.scoreAlignmentAffine(a1,a2,S,o,c), correct_score, "Score of returned alignment is not optimal")

        s, a1, a2 = alignment.SmithWatermanAffine(s2, s1, S, o, c)
        self.assertEqual(alignment_util.scoreAlignmentAffine(a1,a2,S,o,c), s, "Returned alignment score does not match actual alignment score")
        self.assertEqual(alignment_util.scoreAlignmentAffine(a1,a2,S,o,c), correct_score, "Score of returned alignment is not optimal")

    def test2(self):
        s1 = "AAAAAGTGAAAAA"
        s2 = "AAAAACAAAAA"
        S = self.S2
        o = 1
        c = 1
        correct_score = 37
        s, a1, a2 = alignment.SmithWatermanAffine(s1, s2, S, o, c)
        self.assertEqual(alignment_util.scoreAlignmentAffine(a1,a2,S,o,c), s, "Returned alignment score does not match actual alignment score")
        self.assertEqual(alignment_util.scoreAlignmentAffine(a1,a2,S,o,c), correct_score, "Score of returned alignment is not optimal")

    def test3(self):
        s1 = "AAAAAGTGAAAAA"
        s2 = "AAAAACAAAAA"
        S = self.S2
        o = 2
        c = 1
        correct_score = 35
        s, a1, a2 = alignment.SmithWatermanAffine(s1, s2, S, o, c)
        self.assertEqual(alignment_util.scoreAlignmentAffine(a1,a2,S,o,c), s, "Returned alignment score does not match actual alignment score")
        self.assertEqual(alignment_util.scoreAlignmentAffine(a1,a2,S,o,c), correct_score, "Score of returned alignment is not optimal")

    def test4(self):
        s1 = "AAAAATTTAAAAACCCCC"
        s2 = "AAAAAAAAAATTTCCCCC"
        S = self.S3
        o = 5
        c = 2
        correct_score = 128
        s, a1, a2 = alignment.SmithWatermanAffine(s1, s2, S, o, c)
        self.assertEqual(alignment_util.scoreAlignmentAffine(a1,a2,S,o,c), s, "Returned alignment score does not match actual alignment score")
        self.assertEqual(alignment_util.scoreAlignmentAffine(a1,a2,S,o,c), correct_score, "Score of returned alignment is not optimal")

    def test5(self):
        s1 = "GGGAAAAATTTAAAAACCCCC"
        s2 = "AAAAAAAAAATTTCCCCCGGG"
        S = self.S3
        o = 5
        c = 2
        correct_score = 128
        s, a1, a2 = alignment.SmithWatermanAffine(s1, s2, S, o, c)
        self.assertEqual(alignment_util.scoreAlignmentAffine(a1,a2,S,o,c), s, "Returned alignment score does not match actual alignment score")
        self.assertEqual(alignment_util.scoreAlignmentAffine(a1,a2,S,o,c), correct_score, "Score of returned alignment is not optimal")

    def test6(self):
        s1 = "AAAAATTTAAAAACCCCCGGG"
        s2 = "GGGAAAAAAAAAATTTCCCCC"
        S = self.S3
        o = 5
        c = 2
        correct_score = 128
        s, a1, a2 = alignment.SmithWatermanAffine(s1, s2, S, o, c)
        self.assertEqual(alignment_util.scoreAlignmentAffine(a1,a2,S,o,c), s, "Returned alignment score does not match actual alignment score")
        self.assertEqual(alignment_util.scoreAlignmentAffine(a1,a2,S,o,c), correct_score, "Score of returned alignment is not optimal")

    def test7(self):
        s1 = "CCCAAAAATTTAAAAACCCCCGGG"
        s2 = "GGGAAAAAAAAAATTTCCCCCCCC"
        S = self.S3
        o = 5
        c = 2
        correct_score = 128
        s, a1, a2 = alignment.SmithWatermanAffine(s1, s2, S, o, c)
        self.assertEqual(alignment_util.scoreAlignmentAffine(a1,a2,S,o,c), s, "Returned alignment score does not match actual alignment score")
        self.assertEqual(alignment_util.scoreAlignmentAffine(a1,a2,S,o,c), correct_score, "Score of returned alignment is not optimal")
    
    def mcalliTest(self):
        s1 = 'CCCCCCCCCAAAAAAAAAAAAAAAAAAATAAACCCCCCC'
        s2 = 'CCCCCCCCCTTTTTTTTTTTTTTTTTTTTTTTCCCCCCC'
        S = {x:{y:100*self.S2[x][y] for y in "ACGT"} for x in "ACGT"}
        o = 500
        c = 1
        correct_score = 5462
        
        s, a1, a2 = alignment.SmithWatermanAffine(s1, s2, S, o,c)
        self.assertEqual(alignment_util.scoreAlignmentAffine(a1,a2,S,o,c), s, "Returned alignment score does not match actual alignment")
        self.assertEqual(alignment_util.scoreAlignmentAffine(a1,a2,S,o,c), correct_score, "Score of returned alignment is not optimal")
        
        t, b1, b2 = alignment.SmithWatermanAffine(s2, s1, S, o,c)
        self.assertEqual(alignment_util.scoreAlignmentAffine(b1,b2,S,o,c), t, "Returned alignment score does not match actual alignment")
        self.assertEqual(alignment_util.scoreAlignmentAffine(b1,b2,S,o,c), correct_score, "Score of returned alignment is not optimal")
        

def full_unit_test():
    print("WARNING: This unit tester is only a PARTIAL tester -- it is not comprehensive")
    print("or necessarily even finished.  Do not count on it for the full range of testing")
    print("required to verify your algorithms -- only as a starting place.\n")
    suite = fullSuite()
    unittest.TextTestRunner(verbosity=2).run(suite)

if __name__ == "__main__":
    full_unit_test()

