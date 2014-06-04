import unittest
import random
import filecmp
import dnaSeq

def randomSeq(n, s = "ACGT"):
    return "".join([random.choice(s) for i in range(n)])

def fullSuite():
    suite = unittest.TestSuite()
    suite.addTest(InitTestCase('baseTest'))
    suite.addTest(InitTestCase('valueTest'))
    suite.addTest(InitTestCase('dnaTest'))
    suite.addTest(LenTestCase('lengthTest'))
    suite.addTest(StrTestCase('strTest'))
    suite.addTest(GetitemTestCase('getitemRegTest'))
    suite.addTest(GetitemTestCase('getitemNegTest'))
    suite.addTest(GetitemTestCase('getitemRegSliceTest'))
    suite.addTest(GetitemTestCase('getitemNegSliceTest'))
    suite.addTest(GetitemTestCase('getitemSkipSliceTest'))
    suite.addTest(GetitemTestCase('getitemTypeSliceTest'))
    suite.addTest(SetitemTestCase('setitemRegTest'))
    suite.addTest(SetitemTestCase('setitemNegTest'))
    suite.addTest(SetitemTestCase('setitemRegSliceTest'))
    suite.addTest(SetitemTestCase('setitemNegSliceTest'))
    suite.addTest(SetitemTestCase('setitemMultiSliceTest'))
    suite.addTest(AddTestCase('addTest'))
    suite.addTest(AddTestCase('addTypeTest'))
    suite.addTest(ContainsTestCase('containsTest'))
    suite.addTest(ReadFATestCase('readTest'))
    suite.addTest(WriteFATestCase('writeTest'))
    return suite

class SimpleTestCase(unittest.TestCase):
    def setUp(self):
        self.n = random.randint(30, 400)
        self.s = randomSeq(self.n)

class SimpleSeqTestCase(SimpleTestCase):
    def setUp(self):
        super(SimpleSeqTestCase, self).setUp()
        self.seq = dnaSeq.dnaSeq(self.s)

class SimpleTwoSeqTestCase(SimpleSeqTestCase):
    def setUp(self):
        super(SimpleTwoSeqTestCase, self).setUp()
        self.n2 = random.randint(30, 400)
        self.s2 = randomSeq(self.n2)
        self.seq2 = dnaSeq.dnaSeq(self.s2)

class SimpleItemTestCase(SimpleSeqTestCase):
    def setUp(self):
        super(SimpleItemTestCase, self).setUp()
        self.i = random.randint(0,self.n-1)
        self.j = random.randint(0,self.n-1)
        if self.i > self.j:
            self.i,self.j = self.j,self.i

class SimpleFileTestCase(unittest.TestCase):
    def setUp(self):
        self.seq1 = dnaSeq.dnaSeq("".join(['A' for i in range(80)]))
        self.seq1.info = ">ONE"
        self.seq2 = dnaSeq.dnaSeq("".join(['C' for i in range(80)]))
        self.seq2.info = ">TWO"
        self.seqs = [self.seq1, self.seq2]

class InitTestCase(SimpleTestCase):
    """Test that the constructor runs on a legit sequence without exceptions"""
    def baseTest(self):
        dna = dnaSeq.dnaSeq("ACGTacgtNn")
        self.assertIsInstance(dna, dnaSeq.dnaSeq, "Constructor failed, object is not of type dnaSeq")
    
    """Test that the constructor raises a ValueError properly"""
    def valueTest(self):
        self.assertRaises(ValueError, dnaSeq.dnaSeq, 2)
    
    """Test that the constructor raises a DNAError properly"""
    def dnaTest(self):
        self.assertRaises(dnaSeq.DNAError, dnaSeq.dnaSeq, "ABCD")

class LenTestCase(SimpleSeqTestCase):
    """Test that the length method returns the correct number of bases"""
    def lengthTest(self):
        self.assertEqual(len(self.seq), self.n, "Length of dnaSeq did not match generated sequence")
    
class StrTestCase(SimpleSeqTestCase):
    """Test that the string method returns the correct sequence of bases"""
    def strTest(self):
        self.assertEqual(self.s, str(self.seq), "dnaSeq string did not match generated sequence")

class GetitemTestCase(SimpleItemTestCase):
    """Test that the getitem operator[(i > 0)] returns the correct base"""
    def getitemRegTest(self):
        self.assertEqual(self.seq[self.i], self.s[self.i], "dnaSeq[i] did not match the correct base")
    
    """Test that the getitem operator[(i < 0)] returns the correct base"""
    def getitemNegTest(self):
        self.assertEqual(self.seq[-1*self.i], self.s[-1*self.i], "dnaSeq[-i] did not match the correct base")
    
    """Test that the getitem operator[(i > 0):(j > 0)] returns the correct slice"""
    def getitemRegSliceTest(self):
        self.assertEqual(str(self.seq[self.i:self.j]), self.s[self.i:self.j], "dnaSeq[i:j] did not match the generated slice")
    
    """Test that the getitem operator[(j < 0):(i < 0)] returns the correct slice"""
    def getitemNegSliceTest(self):
        self.assertEqual(str(self.seq[-1*self.j:-1*self.i]), self.s[-1*self.j:-1*self.i], "dnaSeq[-j:-i] did not match the generated slice")
    
    """Test that the getitem operator[(i > 0):(j > 0):2] returns the correct slice"""
    def getitemSkipSliceTest(self):
        self.assertEqual(str(self.seq[self.i:self.j:2]), self.s[self.i:self.j:2], "dnaSeq[i:j:2] did not match the generated slice")
    
    """Test that the getitem operator[i:j] returns a dnaSeq"""
    def getitemTypeSliceTest(self):
        self.assertIsInstance(self.seq[self.i:self.j], dnaSeq.dnaSeq, "dnaSeq[i:j] is not of type dnaSeq")

class SetitemTestCase(SimpleItemTestCase):
    """Test that the setitem operator[(i > 0)] changes the correct base"""
    def setitemRegTest(self):
        self.s = self.s[:self.i] + "n" + self.s[self.i+1:]
        self.seq[self.i] = 'n'
        self.assertEqual(str(self.seq), self.s, "dnaSeq[i] = 'n' failed to change the correct base")
        
    """Test that the setitem operator[(i < 0)] changes the correct base"""
    def setitemNegTest(self):
        self.s = self.s[:-1*self.i] + "n" + self.s[-1*self.i+1:]
        self.seq[-1*self.i] = 'n'
        self.assertEqual(str(self.seq), self.s, "dnaSeq[-i] = 'n' failed to change the correct base")
    
    """Test that the setitem operator[(i > 0):(j > 0)] changes the correct slice"""
    def setitemRegSliceTest(self):
        self.s = self.s[:self.i] + "n" + self.s[self.j:]
        self.seq[self.i:self.j] = 'n'
        self.assertEqual(str(self.seq), self.s, "dnaSeq[i:j] = 'n' failed to change base(s)")
    
        """Test that the setitem operator[(j < 0):(i < 0)] changes the correct slice"""
    def setitemNegSliceTest(self):
        self.s = self.s[:-1*self.j] + "n" + self.s[-1*self.i:]
        self.seq[-1*self.j:-1*self.i] = 'n'
        self.assertEqual(str(self.seq), self.s, "dnaSeq[-j:-i] = 'n' failed to change base(s)")
    
    """Test that the setitem operator[(i > 0):(j > 0)] can add multiple bases"""
    def setitemMultiSliceTest(self):
        self.s = self.s[:self.i] + "nnn" + self.s[self.j:]
        self.seq[self.i:self.j] = 'nnn'
        self.assertEqual(str(self.seq), self.s, "dnaSeq[i:j] = 'nnn' failed to change base(s)")

class AddTestCase(SimpleTwoSeqTestCase):
    """Test that concatenating two dnaSeqs matches concatenating the strings used to initialize said dnaSeqs"""
    def addTest(self):
        self.assertEqual(str(self.seq + self.seq2), self.s + self.s2, "The concatenated dnaSeqs do not match the generated concatenated sequence")
    
    """Test that concatenating two dnaSeqs returns an instance of a dnaSeq"""
    def addTypeTest(self):
        self.assertIsInstance(self.seq + self.seq2, dnaSeq.dnaSeq, "The concatenated sequences are not of type dnaSeq")

class ContainsTestCase(SimpleItemTestCase):
    """Test that the contains operator (in) can locate a substring"""
    def containsTest(self):
        self.assertIn(dnaSeq.dnaSeq(self.s[self.i:self.j]), self.seq, "Could not find a substring of bases in the generated dnaSeq")

class ReadFATestCase(SimpleFileTestCase):
    def readTest(self):
        fileseqs = dnaSeq.readFA("tester.fa")
        self.assertTrue(all([str(fileseqs[i]) == str(self.seqs[i]) and self.seqs[i].info in fileseqs[i].info for i in range(len(self.seqs))]), "Info and sequences read from tester.fa did not match expected results")

class WriteFATestCase(SimpleFileTestCase):
    def writeTest(self):
        dnaSeq.writeFA(self.seqs, "generated.fa", 80)
        with open("tester.fa") as fp1, open("generated.fa") as fp2:
            lines1 = fp1.readlines()
            lines2 = fp2.readlines()
            self.assertTrue(len(lines1)==len(lines2), "The generated file appears to have the wrong number of lines")
            
            for line1,line2 in zip(lines1,lines2):
                self.assertTrue(line1.rstrip() == line2.rstrip(), "The generated file appears to have incorrect lines")


def full_unit_test():
    suite = fullSuite()
    unittest.TextTestRunner(verbosity=2).run(suite)

if __name__ == "__main__":
    full_unit_test()
