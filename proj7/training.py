import numpy as np
import argparse
import re

### Using the argsparse to feed option for command line####
parser = argparse.ArgumentParser()
parser.add_argument('-infile', type=argparse.FileType('r'))    ### -i input as most likely path informtion###
parser.add_argument('-seqfile', type=argparse.FileType('r'))  #### -s = input as sequence information path for lambda phage###
parser.add_argument('-obsfile', type=argparse.FileType('r'))  #### -o = input as sequence of observation###3
results = parser.parse_args()
ST = ''  ###Most likely path from viterbi algorithm###
OT = ''  ###Seq of observations#####

try:
    if results.obsfile:
        for ob in results.obsfile:
            match_observation = re.search(r'^observation\s\=\s(.*)',ob)
            if match_observation:
                OT = match_observation.group(1)
                OT = OT.split(',')
 
    sequence = ''
    if results.seqfile:
        for each in results.seqfile:
            my_match = re.search(r'^>.*',each)
            if my_match:
                print "IGNORE FIRST LINE"
            else:
                sequence += each.replace('\n','')
                OT = list(sequence)   ###OT is the observation###
    if results.infile:
        for each in results.infile:
            ST = list(each)               ###ST is the most likely path####
except IOError,msg:
    parser.error("Invalid Input file")

transition_dict = {}
transition_probability = {}
count_dict = {}
n = 1
m=2
for i in ST:
    fill_list = []
    
    if i not in count_dict:
        count_dict[i] = 1
    else:
        count_dict[i] += 1
    if n < len(ST):
   
    
        
        key = i+ST[n]
        n += 1
        if key not in transition_dict:
            transition_dict[key] = 1
        else:
            transition_dict[key] += 1
##print transition_dict
        
#transition_probability.setdefault(i,{}).setdefault(ST[j],count_dict[i])
##print transition_dict

new_dict = {}
SH = ST[:]
for i in ST:
    for j in SH:
        new_dict.setdefault(i,{}).setdefault(j,0)
for each in new_dict:
    for i in new_dict[each]:
        key = each+i
        if key in transition_dict:
            new_dict[each][i] = float(transition_dict[key])/(count_dict[each])
print new_dict ###New transition probability dictionary#####
    
###### update the emission probability######
n = 0
count_states_list = []
count_states_dict = {}
count_obser_dict = {}
zipped = zip(ST,OT)
new_observation_probability = {}
for each in zipped:
    if each[0] not in count_states_dict:
        count_states_dict[each[0]] = [each[1]]
    else:
        count_states_dict[each[0]].append(each[1]) 

for a in count_states_dict:
    if a not in count_obser_dict:
        count_obser_dict[a] = dict(zip(count_states_dict[a],map(count_states_dict[a].count,count_states_dict[a])))
##print count_states_dict

##print count_states_dict
emission_probability={}
for i in ST:
    for j in OT:
        emission_probability.setdefault(i,{}).setdefault(j,0)
##print count_obser_dict
##print count_dict
for i in emission_probability:
    for j in emission_probability[i]:
        if j in count_obser_dict[i]:
            emission_probability[i][j] = float(count_obser_dict[i][j])/(count_dict[i])
print emission_probability ####New emission probability dictionary#####

    
        



        

    
    
         



