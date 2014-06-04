## The script uses Entrez module of biopython to call eUtils program (efetch) and downloads the protein sequences in fasta format.
# Input file should be given as a list of protein primary ID's.

from Bio import Entrez
from Bio import SeqIO
import re
import os

my_email = "ashwani.vit@gmail.com"
fh = open('/home/kumara3/CoV.ids','r')
output_file = "/home/kumara3/protein_fastaseq.fa"
id_dict = []

for each in fh:
    my_match = re.match(r'^\D.*\t(.*)', each,re.M)
    ids = my_match.group(1)
    id_dict.append(ids)
    protein_id = ','.join(id_dict)
if not os.path.isfile(output_file):
    file_handle = Entrez.efetch(db="protein", id = protein_id, rettype="fasta", retmode = "text")
    output_file_handle = open(output_file,'w')
    for eachline in file_handle:
    	output_file_handle.write(eachline)
    #output_file_handle.write(file_handle.read())
output_file_handle.close()
file_handle.close()
    

