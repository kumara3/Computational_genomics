

README

Source Species: Palm civet (accession number : AAV49723.1)
Evidences: From the phylogenetic tree generated, protein sequences of spike gene from palm civet is found to be closest to human sars coronavirus as compared to other species in the tree. To answer the open question of 'from where the virus originated' it appears to be a result of cross jump between either feline,porcine or canine.(from the tree SARS genome appears to be a mixture from different animal coronavirus like cats, dogs, pigs and does not show close relatedness to other humans coronavirus).It could be predicted that SARS would show signs and symptoms of bronchitis, peritonitis.

Code explanation:
download_fasta_proteinseq.py
This scripts download the fasta file.

CoV_dist.py
This script reads the downloaded fasta file, draws the alignment by applying neighbor joining methods, calculates the distance matrix and finally generates the newick string to be used for generating phylogenetic tree.
Function readfasta()
Reads the fasta file. Creates a dictionary called store_seq with sequence description as key and protein sequences as value.
Function find_distance()
Input = sequence_dictionary (dictionary returned by the readfasta function.)
output = distance matrix.
Alignment by SmithWatermanAffine.
The function calls SmithWatermanAffine to calculate the alignment between a pair of sequence. To give the input sequences as parameter to SmithwatermanAffine, I am using two dictionary namely sequence_dictionary and a copied version called copied_dict. The alignment output is returned in a tuple which is then converted to a list called output. Using element two and three of the output list, I am calculating the p distance (abbreviated as p_distance) and then actual distances to be put in distance matrix (abbreviated as evol_distance).The distances are fed into a nested dictionary called distance_matrix which is then returned.

Pickling distance matrix and newick object.

The distance matrix has been pickled in file CoV.dist
PS:The newick object is pickled in file CoV.tree. The newick object which I got had all the sequence description as label. While trying to feed the newick string in newick viewer, it was showing error as “the object should not have more than 50 characters”. So I deleted the description and preserved accession number as label. Tree generated has accession number as label.Below is the newick string which I used to generate tree.







(NP_937950.1:0.0513315984227,(AAL57308.1:0.00384472619636,AAK83356.1:0.0013086683124)L9:0.0365486889083,(((YP_209233.1:0.0597521797342,NP_045300.1:0.0623080355177)L6:0.143390077517,((AAP41037.1:0.00889859625493,AAV49723.1:0.00396000963094)L1:0.887101627804,(NP_040831.1:0.876753716618,(BAA06805.1:0.484724073671,(NP_598310.1:0.60101666848,S41453:0.391388342075)L2:0.0407311861752)L3:0.555861116874)L4:0.600794856232)L5:0.680202292055)L7:0.176507455018,AAL80031.1:0.0920156251066)L8:0.073270850806);



