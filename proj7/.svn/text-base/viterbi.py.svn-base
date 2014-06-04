import numpy as np
import argparse
import re


states = []
#observation = set()
observation = []

start_probability = {}  
T = {}
O = {}
### Using the argsparse to feed option for command line####
parser = argparse.ArgumentParser()
parser.add_argument('-infile', type=argparse.FileType('r'))
parser.add_argument('-seqfile', type=argparse.FileType('r'))
#parser.add_argument('-outfile', type=argparse.FileType('w'))
results = parser.parse_args()
try:
 
    sequence = ''
    if results.seqfile:
        for each in results.seqfile:
            my_match = re.search(r'^>.*',each)
            if my_match:
                print "IGNORE FIRST LINE"
            else:
                sequence += each.replace('\n','')
                observation = list(sequence)
            
    for each in results.infile:
        match_states = re.search(r'^states\s\=\s(.*)',each)
        match_observation = re.search(r'^observation\s\=\s(.*)',each)
        match_start_probability = re.search(r'^start_probability\s\=\s(.*)',each)
        match_transition_probability = re.search(r'^Transition_probability\s\=\s(.*)',each)
        match_emission_probability = re.search(r'^emission_probability\s\=\s(.*)',each)

        if match_states:
            states = match_states.group(1)
            states = list(states)
        if match_observation:
            observation = match_observation.group(1)
            observation = observation.split(',')
        if match_start_probability:
            start_probability = match_start_probability.group(1)
        if  match_transition_probability:
            T = match_transition_probability.group(1)
        if  match_emission_probability:
            O = match_emission_probability.group(1)
except IOError,msg:
    parser.error("Invalid Input file")
    


def most_likely_path(observation,states,T,O,start_probability):
    n = len(observation)
    m = len(states)
    score = 0
    V= np.zeros((n,m))
    backtrace = np.zeros((n,m))
    maximum_score = 0
    max_index = (0,0)
 
    for no,j in enumerate(states):                                              ###Filling the first row of (n,m) matrix###
        
        V[0][no] = np.log((eval(start_probability)[j])) * (eval(O)[j][observation[0]])
        backtrace[0][no] = -1

        
    for o in range(1,n):                                                                ####Filling the viterbi matrix####
        for j in states:
            max_list = []
            increase = 0

            for k in range(0,m):
                trans_state = [s for s in states[0+increase:1+increase]]
                max_list.append((V[o-1,k]*eval(T)[trans_state[0]][j],k))
                increase += 1
                max_tup_list = tuple(max_list)
               # print max_tup_list
            V[o][k],backtrace[o][k] = max(max_tup_list)
            V[o][k] = np.log(V[o][k] * eval(O)[j][observation[o]])
    k = 0
    l = 0
    for obs in range(n-1,n):
        maximum = -10 #V[n-1][0]
        for sts in range(1,m):
            if V[obs][sts] > maximum:
                maximum = V[obs][sts]
            max_index = (obs,sts)
            k = obs
            l = sts
            
    best_path = []
    most_likely_best_path = ""
    index = int(backtrace[k][l])
    best_path.append(states[index])
    
    while backtrace[k][l] != -1:    ### Trace back of pointer matrix to get the most likely path###(-1 denotes the constant values in first row, 0 if the max value came from state1, 1 if max value came from state2 and so on###
        k = k-1
        l = int(backtrace[k][l])
        best_path.append(states[l])
    best_path = best_path[::-1]
    most_likely_best_path = "".join(best_path)
    
    
    print "VITERBI MATRIX = ",V,"\n"
    print "POINTER MATRIX =", backtrace,"\n"
    #print "MOST LIKELY PATH =",most_likely_path,"\n"
    print "SCORE =",score
    with open('output.txt', 'w') as fp:
        fp.write(most_likely_best_path)


   

    

    
   


    
           
                
 
              

 














##states = ['H', 'F']
## 
##observation = ('normal', 'cold', 'dizzy')
## 
##start_probability = {'H': 0.6, 'F': 0.4}
 
##T = {
##   'Healthy' : {'Healthy': 0.7, 'Fever': 0.3},
##   'Fever' : {'Healthy': 0.4, 'Fever': 0.6},
##   }
## 
##O = {
##   'Healthy' : {'normal': 0.5, 'cold': 0.4, 'dizzy': 0.1},
##   'Fever' : {'normal': 0.1, 'cold': 0.3, 'dizzy': 0.6},
##   }


##T = {
##   'H' : {'H': 0.5,'F':0},
##   'F' : {'H': 1.0,'F':0}
##   }
## 
##O = {
##   'H' : {'cold':0.5, 'dizzy':0.5,'normal':0},
##   'F' : {'normal': 1.0,'dizzy':0,'cold':0},
##   }



##states = ['l', 'o']
##start_probability = {'l': 0.9, 'o': 0.1}
##O = {
##    'l': {000:0.2, 001:0.0, 010:0.4, 011:0.2, 100:0.1, 101:0.0, 110:0.1, 111:0.0},
##    'o': {000:0.0, 001:0.05, 010:0.1, 011:0.25, 100:0.05, 101:0.0, 110:0.25, 111:0.3},
##}
##
##observation =(010, 011, 011, 010, 111, 111, 011, 011, 111, 010)
##T = {
##    'l': {'l':0.33, 'o':0.67},
##    'o': {'l':0.40, 'o':0.60},
##}
 



    
##states = ['s1','s2','s3']
##O = { 's1':{'x':0.5,'y':0.5,'z':0}, 's2':{'x':0,'y':0.5,'z':0.5}, 's3':{'x':0.5,'y':0,'z':0.5}}
##states = ['M','N','Q']
##observation = ("xyz")
##
##
##
##start_probability = {'M':0.5,'N':0.5,'Q':0}
####
####T = { 's1':{'s1':0,'s2':0.33,'s3':0.66}, 's2':{'s1':0.33,'s2':0,'s3':0.66}, 's3':{'s1':0.33,'s2':0.33,'s3':0.33}}
##
##
##T = {'M': {'N': 1.0}, 'N': {'Q': 1.0}}
##O = {'Q': {'z': 1.0}, 'M': {'x': 1.0}, 'N': {'y': 1.0}}
print most_likely_path(observation,states,T,O,start_probability)
