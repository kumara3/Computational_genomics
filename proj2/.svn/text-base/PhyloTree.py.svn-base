# Impelementation of the PhyloLab tree: Lab 5.
import copy
class PhyloTree:
    
    """Phylogenetic tree manipulation and Newick conversion"""
    
    #################################################
    # PhyloTree constructor.  You may add to this as needed, but may not change the code
    # provided.
    def __init__(self, label = None, children = None):
        """Contructor: create a root node with the specified label and child set.
        Assumes children is a list of tuples, with each tuple containing a 
        PhyloTree object and a weight."""

        self.label = label if label else ""
        self.children = children if children else []

    ##################################################
    # The following have been implemented for you -- DO NOT CHANGE
    def rootLabel(self):
        """Return the label of the root."""
        return self.label

    def child(self, i):
        """Get the ith child"""
        return self.children[i][0]

    def numChildren(self):
        """Return the number of children"""
        return len(self.children)

    def isLeaf(self):
        """Indicates whether a tree has children"""
        return self.numChildren() == 0

    
    #################################################
    # 1) We will define the length of a tree as the number of leaves in the tree.
    def __len__(self):
        
        """Return the number of leaves in the tree."""
        count = 0
        if  self.isLeaf() == True:
            return 1
        else:
            for  p,q in self.children:
                count = count+len(p)                    
        return count
   #################################################
    # 2) Print the string in Newick format.
    def __str__(self):
        """Return the Newick representation of the tree"""
        if self.isLeaf():
            return self.label+';'
        else:
            return  '"'+'('+", ".join([str(node)+(':')+str(int(weight)) for node, weight in self.children])+')'+self.label+';'+'"'
                
    #################################################
    # 3) Return a deep copy of the tree.  That is, the resulting tree should contain entirely
    #    new node objects (so modifying a node in one tree does not modify nodes in the other), 
    #    but should be isomorphic ("identical") to the source.
    def copy(self):
        """Return a deep copy of the tree"""
        new = []
        if self.isLeaf():
            newobject = PhyloTree(self.label)
            return newobject
        else:
            for k,v in self.children:
                new.append((k.copy(),v))
            new_obj = PhyloTree(self.label,new)
            return new_obj
                
             
        #newobject = copy.deepcopy(self)
##        newobject.label = 'Z'
##        for new_node, new_wt in newobject.children:
##            return  '('+", ".join([str(new_node)+(':')+str(int(new_wt)) for node, weight in self.children])+')'+self.label+';'
##            #print '%c%s%c%d%c%s' %('(',new_node,':',new_wt,')',self.label)
##            ## string formating required
            

    #################################################
    # 4) Return a set of all leaf nodes in the tree.  (The set should contain PhyloTree objects,
    # not just labels.
    def leafSet(self):
        """Return a set of the leaf nodes"""
        leaf_set = set()
        if self.isLeaf():
            leaf_set.add(self)
            return leaf_set
        else:
            for t,wt in self.children:
 #               value = t.leafSet()
                leaf_set = leaf_set.union(t.leafSet())
        return leaf_set
        
    #################################################
    # 5) Return the sum of all edges in the tree.
    def treeWeight(self):
        edge_weight = 0
        if self.isLeaf():
            return edge_weight
        else:
            for k, m in self.children:
                
                edge_weight += k.treeWeight() + m 
        return edge_weight
        """Return the sum of the weight of all edges in the tree"""
            

    #################################################
    # 6) Return a dictionary that maps the labels of the leaves to their dpeth in the tree.
    def mapDepth(self):
        """Return a dictionary mapping leaf labels to their depth in the tree."""
       # depth = 0
        mapping_leaves = dict()
        if self.isLeaf():
            mapping_leaves[self.label] = 0
            return mapping_leaves
        else:
            for x,y in self.children:
                depth = y
                x = PhyloTree(x)
                mapping_leaves[x.label] = depth + x.mapDepth()
        return mapping_leaves
        
    #################################################
    # 7) Return a dictionary mapping each leaf-label of the tree to a list of
    #    the nodes on the path from the root to the leaf (in order).
    def mapPath(self,path_list=[],map_path=dict()):

        if self.isLeaf():
            map_path[self.label] = [self.label]
            path_list.append(self.label)
            return map_path
        else:
            path_list.append(self.label)
            for m,n in self.children:
                L = m.mapPath(path_list)
                #L[L.keys()[0]] = path_list.append(L[L.keys()[0]])
                #map_path[map_path.keys()].append
                map_path.update(L)
        return map_path


    #################################################
    # 8) Return the MRCA (Most Recent Common Ancestor) of the two labels.  That is,
    #    returns the node with the greatest depth that is an ancestor to both leaves. 
    #    D, if provided, is the result of the mapPath call (so that does not need to be 
    #    re-calculated).   Throw a ValueError is one of the labels is not in the tree.
    def MRCA(self, label1, label2,D=None):
    	
   	for x in D.keys():
		if x== label1:
			return True
		elif X==label2:
			return True
			
		else:
			raise ValueError
   	MRCA_list = []
	list_one = D[label1]
	list_two = D[label2]
	for i in reversed(list_one):
		
		for j in reversed(list_two):
			if i.label == j.label:
				return i
	

   	
			


    #################################################
    # 9) Returns a copy of the tree that has been reoriented based on an outgroup.
    #    Input:
    #    * outgroup_label: The label of one leaf that will be a child of the root in the
    #      reoriented tree. 
    #    * split_fraction: The split edge will have 100*split_fraction% of its weight
    #      distriubted to the r->outgroup edge, and the remaining distrubted to the other
    #      new edge.
    #    * root_label: The label of the new root node.
##    def rootByOutgroup(self, outgroup_label, split_fraction = 0.5, root_label = ""):
##        assert 0 < split_fraction < 1
##        pass
#A = PhyloTree( label = 'A', children = [('C',2), ('D',1)]) 
#C = PhyloTree (label = 'C')
#H = PhyloTree(label = 'H', children = [('F',1),('G',2)])
##A = PhyloTree("A")                                        # A;
##B = PhyloTree("B")                                        # B;
##C = PhyloTree("C", [(A,1),(B,1)])
##print C.treeWeight()
