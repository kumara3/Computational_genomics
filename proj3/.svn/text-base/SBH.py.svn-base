#def SBH(S):
    # Input: A multi-set S (represented as a list) of length l > 1 strings representing dna sequences.
    # Output: A string s representing a dna sequence.
    # Goal: Spectrum(s,l) == S when s is viewed as a circular genome.
    #   Note: The above should be multi-set equality, not list equality.  (e.g. {"AAA","CCC", "CCC"} == {"CCC", "AAA", "CCC"})

store_fragment = []
def help_in(fragment_list, element_1, element_2,store_fragment):
    if len(store_fragment) > len(fragment_list):
        print store_fragment
        return calculate_overlap(store_fragment)


    for i,itr in enumerate(fragment_list):
            if element_2[0]  ==  itr[0:-1] :
                #print "Sequence not found in " +str(i)+ "round"
                store_fragment.extend([itr[1:]])
                element_2[:] = [itr[1:]]
                return help_in(fragment_list,element_1,element_2,store_fragment)
            else:
                if element_2[0] == str(itr[0]):
                    store_fragment.extend([itr[1:]])
                    element_2[:] = [itr[1:]]
                    return help_in(fragment_list,element_1,element_2,store_fragment)
                    
def calculate_overlap(fragments):
    sequence = str(fragments[0])
    
    for k in range(1,(len(fragments))):
        sequence += (fragments[k][-1])
    result = overlap_helper(sequence)
    print result

def overlap_helper(search_factor):
    for i in range(1,len(search_factor)):
        if search_factor[0] == search_factor[i] and search_factor[i:] == search_factor[:len(search_factor[i:])]:
            return search_factor[:i]
        
        
        
        
def SBH(S):
    fragment_list = list(S)
    first_element = fragment_list[0]
    element_1 = [first_element[0:-1]]
    element_2 = [first_element[1:]]
    store_fragment.extend((element_1[0],element_2[0]))
    return help_in(fragment_list,element_1, element_2,store_fragment)

       
#SBH(S=set(["ATTAC","TACAG","GATTA","ACAGA","CAGAT","TTACA","AGATT"]))
#SBH(S= set(['AC', 'CG', 'GT', 'TA', 'AC']))
#SBH(S= set(['AC', 'CG', 'GT', 'TA', 'AC', 'CCAAAAA', 'CAAAAAA', 'AAAAAAA', 'AAAAAAC', 'AAAAACC', 'AAAACCC', 'AAACCCC', 'AACCCCC', 'ACCCCCC', 'CCCCCCC', 'CCCCCCA', 'CCCCCAA', 'CCCCAAA', 'CCCAAAA', 'CCAAAAA']))


