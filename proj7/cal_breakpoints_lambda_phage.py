##Calculation of breakpoints in lambda phage###

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-seqfile', type=argparse.FileType('r'))  #### -s = input as most likle path###
#parser.add_argument('-outfile', type=argparse.FileType('w'))
results = parser.parse_args()

n = 0
m = 1
High_GC_followed_Low_GC = []
Low_GC_followed_High_GC = []
breakpoints = []
for each in results.seqfile:
    if each[n] == 'H' and each[m] == 'L':
        High_GC_followed_Low_GC.append(len(each[0:n]))
    elif each[n] == 'L' and each[m] == 'H':
        Low_GC_followed_High_GC.append(len(each[0:n]))
    breakpoint = [High_GC_followed_Low_GC]+[Low_GC_followed_High_GC]
#print breakpoint

with open('breakpoints_from_lambda.txt', 'w') as fp:
    fp.write(' ,'.join([str(i) for i in breakpoint]))
    fp.write("\n")
    fp.write("The first coordinate in the list corresponds to High GC content followed by Low GC content\n")
    fp.write("The second coordinate in the list corresponds to Low GC content followed by high GC content\n")    
        
        
        
    
