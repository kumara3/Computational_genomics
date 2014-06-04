import sys
import unittest
import re
import PhyloTree

def stripSpace(s):
    return re.sub("\s", "", str(s))

def leafLabels(s):
    return set(x.label.strip() for x in s.leafSet())

def sameTrees(T1, T2):
    """Determined if the trees are the same.
    Assumes every leaf has a (non-empty) unique label"""
    if T1.label != T2.label:
        return False;
    if T1.numChildren() != T2.numChildren():
        return False
    if T1.isLeaf():  
        return True

    C1 = sorted(T1.children, key = lambda x: x[0].label)
    C2 = sorted(T2.children, key = lambda x: x[0].label)
    return all([sameTrees(c1,c2) and abs(float(w1) - float(w2)) < 0.0001 for (c1,w1),(c2,w2) in zip(C1,C2)])
        

def fullSuite():
    suite = unittest.TestSuite()
    suite.addTest(LengthTestCase('lenTest'))
    suite.addTest(StringTestCase('strTest'))
    suite.addTest(CopyTestCase('copyTest'))
    suite.addTest(LeafSetTestCase('leafSetTest'))
    suite.addTest(TreeWeightTestCase('treeWeightTest'))
    suite.addTest(MapDepthTestCase('mapDepthTest'))
    suite.addTest(MapPathTestCase('mapPathTest'))
    suite.addTest(MRCATestCase('mrcaTest'))
    suite.addTest(RootByOutgroupTestCase('rootByOutgroupTest'))
    return suite

class SimpleTestCase(unittest.TestCase):
    def setUp(self):
        self.A = PhyloTree.PhyloTree("A")                                        # A;
        self.B = PhyloTree.PhyloTree("B")                                        # B;
        self.C = PhyloTree.PhyloTree("C", [(self.A,1),(self.B,1)])               # (A:1,B:1):C;
        self.D = PhyloTree.PhyloTree("D")                                        # D;
        self.E = PhyloTree.PhyloTree("E", [(self.D,1.5)])                        # (D:1.5)E;
        self.F = PhyloTree.PhyloTree("F", [(self.E,2), (self.C,3)])              # ((D:1.5)E:2,(A:1,B:1)C:3)F;
        self.G = PhyloTree.PhyloTree("G", [(self.A,10),(self.B,20),(self.E,30)]) # (A:10,B:20,(D:1.5)E:30)G;
        self.U = PhyloTree.PhyloTree(children = [(self.A,27),(self.B,92)])       # (A:27,B:92);

        self.X1 = PhyloTree.PhyloTree('G', children = [(self.B,20), (self.E,30)])
        self.X2 = PhyloTree.PhyloTree(children = [(self.X1,5),(self.A,5)])
        self.X3 = PhyloTree.PhyloTree('R2', children = [(self.X1,8),(self.A,2)])
    
class LengthTestCase(SimpleTestCase):
    def lenTest(self):
        self.assertEqual(len(self.A), 1, 'Length of PhyloTree("A;") did not equal 1.')
        self.assertEqual(len(self.B), 1, 'Length of PhyloTree("B;") did not equal 1.')
        self.assertEqual(len(self.C), 2, 'Length of PhyloTree("(A:1,B:1):C;") did not equal 2.')
        self.assertEqual(len(self.D), 1, 'Length of PhyloTree("D;") did not equal 1.')
        self.assertEqual(len(self.E), 1, 'Length of PhyloTree("(D:1.5)E;") did not equal 1.')
        self.assertEqual(len(self.F), 3, 'Length of PhyloTree("((D:1.5)E:2,(A:1,B:1)C:3)F;") did not equal 3.')
        self.assertEqual(len(self.G), 3, 'Length of PhyloTree("(A:10,B:20,(D:1.5)E:30)G;") did not equal 3.')
        self.assertEqual(len(self.U), 2, 'Length of PhyloTree("(A:27,B:92);") did not equal 2.')
    
class StringTestCase(SimpleTestCase):
    def strTest(self):
        self.assertEqual(stripSpace(self.A), "A;", 'str(A) did not return "A;"')
        self.assertEqual(stripSpace(self.B), "B;", 'str(B) did not return "B;"')
        self.assertEqual(stripSpace(self.C), "(A:1,B:1)C;", 'str(C) did not return "(A:1,B:1)C;"')
        self.assertEqual(stripSpace(self.D), "D;", 'str(D) did not return "D;"')
        self.assertEqual(stripSpace(self.E), "(D:1.5)E;", 'str(E) did not return "(D:1.5)E;"')
        self.assertEqual(stripSpace(self.F), "((D:1.5)E:2,(A:1,B:1)C:3)F;", 'str(F) did not return "((D:1.5)E:2,(A:1,B:1)C:3)F;"')
        self.assertEqual(stripSpace(self.G), "(A:10,B:20,(D:1.5)E:30)G;", 'str(G) did not return "(A:10,B:20,(D:1.5)E:30)G;"')
        self.assertEqual(stripSpace(self.U), "(A:27,B:92);", 'str(U) did not return "(A:27,B:92);"')
    
class CopyTestCase(SimpleTestCase):
    def copyTest(self):
        F2 = self.F.copy()
        F2.children[0][0].children[0][0].label = 'Z';
        self.assertEqual(stripSpace(F2), "((Z:1.5)E:2,(A:1,B:1)C:3)F;", 'A label in the copied tree could not be changed.')
        self.assertEqual(stripSpace(self.D), "D;", 'Changing a label in the copied tree changed a label in another tree.')
    
class LeafSetTestCase(SimpleTestCase):
    def leafSetTest(self):
        self.assertEqual(leafLabels(self.A), {'A'}, "A.leafSet() did not return {'A'}")
        self.assertEqual(leafLabels(self.B), {'B'}, "B.leafSet() did not return {'B')")
        self.assertEqual(leafLabels(self.C), {'A','B'}, "C.leafSet() did not return {'A','B'}")
        self.assertEqual(leafLabels(self.D), {'D'}, "D.leafSet() did not return {'D'}")
        self.assertEqual(leafLabels(self.E), {'D'}, "E.leafSet() did not return {'D'}")
        self.assertEqual(leafLabels(self.F), {'A','B','D'}, "F.leafSet() did not return {'A','B','D'}")
        self.assertEqual(leafLabels(self.G), {'A','B','D'}, "G.leafSet() did not return {'A','B','D'}")
        self.assertEqual(leafLabels(self.U), {'A','B'}, "U.leafSet() did not return {'A','B'}")

class TreeWeightTestCase(SimpleTestCase):
    def treeWeightTest(self):
        self.assertEqual(self.A.treeWeight(), 0, "A.treeWeight() did not match expected value (0)")
        self.assertEqual(self.C.treeWeight(), 2, "C.treeWeight() did not match expected value (2)")
        self.assertLess(abs(self.E.treeWeight() - 1.5), 0.0001, "E.treeWeight() did not fall within acceptable tolerance (1.5 +/- 0.0001)")
        self.assertLess(abs(self.F.treeWeight() - 8.5), 0.0001, "F.treeWeight() did not fall within acceptable tolerance (8.5 +/- 0.0001)")

class MapDepthTestCase(SimpleTestCase):
    def mapDepthTest(self):
        DM = self.F.mapDepth()
        self.assertEqual(DM['A'], 4, "F.mapDepth()->A did not match expected value (4)")
        self.assertEqual(DM['B'], 4, "F.mapDepth()->B did not match expected value (4)")
        self.assertEqual(DM['D'], 3.5, "F.mapDepth()->D did not match expected value (3.5)")
        self.assertEqual(set(DM), set("ABD"), "F.mapDepth() did not return the right set of items, expected {A,B,D}")

class MapPathTestCase(SimpleTestCase):
    def mapPathTest(self):
        DM = self.F.mapPath()
        self.assertEqual(DM['A'], [self.F, self.C, self.A], "F.mapPath()->A did not match expected path (F, C, A)")
        self.assertEqual(DM['D'], [self.F, self.E, self.D], "F.mapPath()->D did not match expected path (F, E, D)")

class MRCATestCase(SimpleTestCase):
    def mrcaTest(self):
        self.assertEqual(self.F.MRCA('A','B'), self.C, "F.MRCA(A,B) did not return expected node (C)")
        self.assertEqual(self.F.MRCA('A','D'), self.F, "F.MRCA(A,D) did not return expected node (F)")
        self.assertEqual(self.G.MRCA('A','D'), self.G, "G.MRCA(A,D) did not return expected node (G)")

class RootByOutgroupTestCase(SimpleTestCase):
    def rootByOutgroupTest(self):
        self.assertTrue(sameTrees(self.G.rootByOutgroup('A'), self.X2), "G.rootByOutgroup('A') did not return the proper tree")
        self.assertTrue(sameTrees(self.G.rootByOutgroup('A', split_fraction = 0.2, root_label = 'R2'), self.X3) or 
                        sameTrees(self.G.rootByOutgroup('A', split_fraction = 0.8, root_label = 'R2'), self.X3), "G.rootByOutgroup('A', fract=0.8, label='R2') did not return the proper tree (X1)")

    
def full_unit_test():
    suite = fullSuite()
    unittest.TextTestRunner(verbosity=2).run(suite)

if __name__ == "__main__":
    full_unit_test()

