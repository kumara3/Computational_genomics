from PhyloTree import *
import PhyloTree
import PhyloTree_sol
import parseNewick
import PhyloTree_util
from operator import itemgetter
#modified_dict = {}
#used = []
def UPGMA(D):
    #"""
    #UPGMA Algorithm:
    #Input:
    #* DIST: A dictinary mapping label pairs to distnce (e.g. for x,y in leafLabels: D[x][y] 
            #is the tree dist. from x to y, where x and y are node labels.
    #Assumptions: 
    #* For every x,y in D.keys(), D[x][y] is defind.
    #* DIST[x][x] == 0 and DIST[x][y] == DIST[y][x].
    #* DIST must describe a tree with the ultrametric property.
    #Return: A PhloyTree object T representing a full, unrooted, binary, ultrametric tree with leaf 
            #lables equal to D.keys(), such that the path distance from leaf x to leaf y is D[x][y] for 
            #any x,y.
    #"""
    modified_dict = {}
    used = []
    
    
    new_list = []
    node1 = ''
    node2 = ''

    print(D)
    while len(D) > 1:  
        print(D)
        mini = 99999
        children_list = []
        for k in D.keys():
            for z in D.keys():
                if D[k][z] < mini and D[k][z] != 0:
                    mini = D[k][z]
                    node1 = k
                    node2 = z
                    height = D[k][z]/2       
        if node1 not in modified_dict:
            A = PhyloTree.PhyloTree(node1)
            A.h = 0
            modified_dict[node1] = A
        if node2 not in modified_dict:
            B = PhyloTree.PhyloTree(node2)
            B.h = 0
            modified_dict[node2] = B
        distance_1 = height - modified_dict[node1].h
        #print modified_dict[node1].label,distance_1
        distance_2 = height - modified_dict[node2].h
        #print modified_dict[node2].label,distance_2
        child_first = (modified_dict[node1],distance_1)
        child_second = (modified_dict[node2],distance_2)    
        children_list.append(child_first)
        children_list.append(child_second)
        object_new = PhyloTree.PhyloTree(label='', children= children_list)
        object_new.h = height
        modified_dict[node1+node2] = object_new
        
                   
        element = node1+node2
        used.append(node1)
        used.append(node2)
    
        D[element] = {}
        D[element][element] = 0
        for each_element in D.keys():
            if each_element not in used and each_element!= element:
                D[element][each_element] = (D[each_element][node1]+D[each_element][node2])/2
                D[each_element][element] = (D[each_element][node1]+D[each_element][node2])/2

        del D[node1]
        del D[node2]         
        for it in D.keys():
            if it!=element:
                del D[it][node1]
                del D[it][node2]  
    return object_new           

#Z = "(((((B:3,C:3):8,D:11):1,(E:6,F:6):6):4,G:16):2,A:18);"
#T = parseNewick.parseNewick(Z)
#D = PhyloTree_util.createDistMatrix(T)
#D = {'A': {'A': 0, 'D': 36.0, 'G': 36.0, 'CB': 36.0, 'EF': 36.0}, 'D': {'A': 36.0, 'D': 0, 'G': 32.0, 'CB': 22.0, 'EF': 24.0}, 'G': {'A': 36.0, 'D': 32.0, 'G': 0, 'CB': 32.0, 'EF': 32.0}, 'CB': {'A': 36.0, 'D': 22.0, 'G': 32.0, 'CB': 0, 'EF': 24.0}, 'EF': {'A': 36.0, 'CB': 24.0, 'EF': 0, 'D': 24.0, 'G': 32.0}}

#print T
#print len(D)
#print UPGMA(D)

#used_2 = []
def NJ(D):
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
    used_2 = []
    all_node = {}
    for each in D.keys():
        objects = PhyloTree.PhyloTree(each)
        all_node[each] = objects    
    while len(D) > 3:
        all_cluster = {}
        sum_r = 0
        node_1 = ''
        node_2 = ''
    
        for my in D.keys():
            sum_r = 0
            for my_d in D.keys():
                sum_r += D[my][my_d]
            all_cluster[my] = sum_r/(len(D)-2)
                #print all_cluster
       
        
        
    
        
        
        children_list = []
        minimum = 9999
        for ma in D.keys():
            for ja in all_cluster:
                if (D[ma][ja] - all_cluster[ma] -all_cluster[ja]) < minimum and D[ma][ja] != 0:
                    minimum = D[ma][ja]- all_cluster[ma] -all_cluster[ja]
                    node_1 = ma
                    node_2 = ja
        length = (D[node_1][node_2] + all_cluster[node_1]-all_cluster[node_2])/2
        children_1 = (all_node[node_1],length)
        children_2 = (all_node[node_2],D[node_1][node_2] - length)
        children_list.append(children_1)
        children_list.append(children_2)
        #children_list = [(all_node[node_1], length), (all_node[node_2],length)]
        object_1 = PhyloTree.PhyloTree(label='', children = children_list)
        all_node[node_1+node_2] = object_1
        del all_node[node_1]
        del all_node[node_2]
        used_2.append(node_1)
        used_2.append(node_2)
        element = node_1+node_2
        D[element] = {}
        D[element][element] = 0
        for b in D.keys():
            if b not in used_2 and b != element:
                D[element][b] = (D[b][node_1] + D[b][node_2] - D[node_1][node_2])/2
                D[b][element] = (D[b][node_1] + D[b][node_2] - D[node_1][node_2])/2
        del D[node_1]
        del D[node_2]
    
        for fa in D.keys():
            if fa != element:
                del D[fa][node_1]
                del D[fa][node_2]
    #return D

    #C = PhyloTree.PhyloTree(label = '', children = [])
    x = ''
    y = ''
    z = ''
    x = D.keys()[0] 
    y = D.keys()[1]
    z = D.keys()[2]
    
    edge_1 = 0
    edge_2 = 0
    edge_3 = 0
    
    edge_1 = D[x][y] 
    edge_2 = D[x][z]
    edge_3 = D[y][z] 
    
    dist_xr = 0
    dist_yr = 0
    dist_zr = 0
    dist_xr = (edge_2+edge_3 -edge_1)/2
    dist_yr = (edge_2 + edge_1 - edge_3)/2
    dist_zr =  (edge_1 + edge_3 - edge_2)/2
    L = PhyloTree.PhyloTree(" ",[(all_node[x],dist_xr),(all_node[y],dist_yr),(all_node[z],dist_zr)])
    return L
            
         
            
            
            
#Z = "(((((B:3,C:3):8,D:11):1,(E:6,F:6):6):4,G:16):2,A:18);"
#Z = "((A:3,B:4):1,(C:5,D:6):1);" 
#T = parseNewick.parseNewick(Z)
#D = PhyloTree_util.createDistMatrix(T)
#L = UPGMA(D)
#print L
##print D
