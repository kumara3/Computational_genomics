import re
import textwrap
BASES = {'A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n'}
COMPLEMENT = {x:y for x,y in zip(['A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n'],
                                  ['T', 'G', 'C', 'A', 'N', 't', 'g', 'c', 'a', 'N'])}
class DNAError(Exception):
    def __init__(self, value = "Bad DNA string"):
        self.value = value

    def __str__(self):
        return str(self.value)


class dnaSeq:
    """Mutable DNA Sequence class."""
    
    # 1) Implement the constructor.  Should take one argument (other than self) and throw a ValueError if
    #    the value is not a string, or a DNAError (see above) if the value is a string but contains bad
    #    bases characters.  (That is, characters that are not in the BASES set.)  Otherwise, store the string
    #    is an attribute in whatever for you see fit.
    #    COMMENT: Since this is intended to be a mutable class, you might want to consider holding
    #    the information in a mutable data strcture -- as opposed to one that you have to copy
    #    every time you want to change the object.
    def __init__(self, s):
        #pass
        self.sequence = s
        
        if type(self.sequence) != str:
            raise ValueError("The sequence is not a string")

        elif not type(self.sequence) == str and not any([bases in BASES for bases in self.sequence]):
            try:
                raise DNAError
            except DNAError,e:
                print e
                
        else:
            self.sequence = [i for i in s]

    def __len__(self):
        return len(self.sequence)
            
        
    def __str__(self):
        return "".join(self.sequence)

#---------------------------------------
    #Implement the __getitem__ function so that we can us the index the operator on a dnaSeq --
    #    both single indexing and splicing.  An index should return a single-character string,
    #    while a slice should return a dnaSeq object.
    #    NOTE: The Python specifications for __getitem__ require that:
    #          1) An illegal index type result in a TypeError being raised.
    #          2) An invalid index value results in a KeyError being raised.
    #    You may assume this for any other object, and implement this for dnaSeq.
    def __getitem__(self, i):
        if type(i) == int:
            if not (-1*len(self.sequence)<=i<len(self.sequence)):
                raise TypeError
            else:
                return self.sequence[i]
        else:
            return dnaSeq(''.join(self.sequence[i]))

#------------------------------------------------------------------------
    #Implement the __setitem__ function so we can use index-assignment (e.g. o[5] = 'A' or O[5:7] = 'AC').
    #    Should throw the same sorts of errors as __getitem__.  Additionally, should throw a ValueError
    #    if the value is not itself either a dnaSeq object or a string, and a DNAError is it is a string
    #    containing bad bases.
    def __setitem__(self, i, value):
            if type(i) == int:
                if not (abs(i)<=len(self.sequence)):
                    raise TypeError
                elif type(i) != int:
                    raise KeyError()
		else:
		    self.sequence[i] = value
            elif not type(value) == str and not isinstance(value,dnaSeq):
                raise ValueError
            elif  not all ([bases in BASES for bases in value]):#self.sequence[i]]):
		raise DNAError
            	   
            elif type(i) == slice:
                value_object = dnaSeq(value)
                self.sequence[i] = value_object.sequence




    #Impelement the __add__ function to enable + to concatenate sequences.
    #    Throw a ValueError if object other is not of class dnaSeq
    #    Make sure it produces a deep copy.

    def __add__(self, other):
        if not type(other) == type(self):
            raise ValueError
        else:
            self.sequence = self.sequence+other.sequence
            return self


    #Implement the necessary method for the "in" operator to to allow operations such as:
    #     s in o (where s is a string or dnaSeq, o is a dnaSeq, and this returns true if either
    #     o contains s as a subsequence.
    #     Raise a ValueError if the query sequence is neither a string nor a dnaSeq.

    def __contains__(self, query):
 
        if type(query) != str and type(query) != type(dnaSeq):
            raise ValueError
        my_match = re.search(str(query),str(self))            

        if my_match:
            return True
        else:
            return False


        
    def complement(self):
        complemented = [COMPLEMENT[z] for z in self.sequence]
        return "".join(complemented)
	
    def reverse(self):
        return "".join(self.sequence[::-1])

    def reverse_complement(self):
        reverse_of_seq = self.sequence[::-1]
        complement_of_reverse_seq = [COMPLEMENT[x] for x in reverse_of_seq]
        return "".join(complement_of_reverse_seq)

sequences = []
def readFA(fasta_file):
	#sequences = []
    	info = ''
	with open(fasta_file) as file:
        	line = file.readline()[:-1]
        	while line:
            		seq = ''
            		if line[0] == '>':
				if info:
					fasta_object = dnaSeq(seq)
					fasta_object.info = line[:-1]
					sequences.append(fasta_object)
                	else:
                    		seq += line

                	line = file.readline()[:-1]
     	return sequences
        

def writeFA(S,output, col_width=60):
    output_file = open(output,'w')
    for k in range(len(S)):
	if S[k].info:
		output_file.write('%s\n' %(S[k].info))
        	output.write("\n".join(S[k].sequence))
    output_file.close()


##                if each[0]:
##                        output.write('%s\n' % each[0])
##                        #for i in xrange (0, len(each[1]), 60):
##                        output.write("\n".join(str(each[1]),col_width))
##                        output.write("\n")
##                else:
##                        print "Description not found"
##        
##        

readFA('sample.fa')            
writeFA(sequences,'example.fa')
#dnaseq_object = dnaSeq("acgtcg")
#dnaseq_object2 = dnaSeq("aa")
