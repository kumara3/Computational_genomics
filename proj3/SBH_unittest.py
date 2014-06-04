import sys
import unittest
import re
import random
from SBH import SBH

def fullSuite():
    suite = unittest.TestSuite()
    suite.addTest(SimpleTestCase('basicTest'))
    suite.addTest(SimpleTestCase('repTest'))
    suite.addTest(SimpleTestCase('randomTest'))
    return suite

class SimpleTestCase(unittest.TestCase):
    def setUp(self):
        self.s1 = "ACGT"

    def spectrum(self, s, l, shuffle = True):
        s2 = s+s
        S = [s2[i:i+l] for i in range(len(s))]
        if shuffle:
            random.shuffle(S)
        return S
        
    def TestSBH(self, s, l, shuffle = True):
        S1 = self.spectrum(s, l, shuffle)
        s2 = SBH(S1)
        S2 = self.spectrum(s2, l, shuffle)
        return sorted(S1)==sorted(S2)

    def basicTest(self):
        self.assertTrue(self.TestSBH("ACGT", 3, False))
        self.assertTrue(self.TestSBH("AAAACCCCGGGG", 5, False))
        
    def repTest(self):
        self.assertTrue(self.TestSBH("AAAAAAAACCCCCCCCAAAAAAAAGGGGGGGGAAAAAAAACCCCTTTT", 8))

    def randomTest(self):
        self.assertTrue(self.TestSBH("".join([random.choice("ACGT") for i in range(500)]), 15))

                    
        


    
def full_unit_test():
    suite = fullSuite()
    unittest.TextTestRunner(verbosity=2).run(suite)

if __name__ == "__main__":
    full_unit_test()

