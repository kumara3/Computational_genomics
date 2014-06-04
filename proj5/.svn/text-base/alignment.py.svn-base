import  alignment_util
import numpy as np

def NeedlemanWunsch(seq1, seq2, S, gap):
    """
    Input:
    * seq1, seq2: Two *strings* (not dnaSeqs) to be aligned.
      -- May not necessarily represent a dna sequence.
    * S: a scoring dictionary such that for any characters a,b: S[a][b] is the miss-match penalty for aligning
         an a to a b.
      -- You may assume that every character contained in seq1 or seq2 is defined in S.
    * gap: The gap penalty, listed as a non-negative number.  (e.g. If gap=4, we will deduct 4 points for each -.)
    Output: A tuple containing three elements: 
    1) The optimal global alignment score.
    2) The corresponding alignment string of seq1.  (That is, a copy of seq1 string with "-" characters added corresponding to an optimal alignment.
    3) The corresponding alignment string of seq2.
    """
    assert gap > 0, "Gap should be a positive number"
    seq1 = " "+seq1
    seq2 =  " "+seq2
    n, m = len(seq1), len(seq2)

    M = [[ float(-gap*j) for j in range(0, m) ]]    # Score matrix
    P = [[0] + [ 2 for i in range(1,m) ]]   # Pointer matrix -- record whether the score vame from the top-left (0), top (1), or left (2).
    for i in range(1, n):
        M.append([0]*m)
        P.append([-1]*m)
        M[i][0] = float(-gap*i)
        P[i][0] = 1
        
        for j in range(1,m):
            t1 = M[i-1][j-1] + S[seq1[i]][seq2[j]]
            t2 = M[i-1][j] - gap
            t3 = M[i][j-1] - gap
            M[i][j], P[i][j] = max((t1,0), (t2,1), (t3,2))  # max will use the tuple with the largest first element
            
            if M[i][j] < 0:
                M[i][j] == 0
    i, j = n-1, m-1
    s1, s2 = "", ""
    while i > 0 or j > 0:
        if P[i][j] == 0:
            s1 += seq1[i]
            s2 += seq2[j]
            i -= 1
            j -= 1
        elif P[i][j] == 1:
            s1 += seq1[i]
            s2 += '-'
            i -= 1
        elif P[i][j] == 2:
            s1 += '-'
            s2 += seq2[j]
            j -= 1
        else:
            assert(False)
    
    return M[n-1][m-1], s1[::-1], s2[::-1]

#########S = alignment_util.readScoringMatrix("/Users/ashwanikumar/Documents/project_5/Blosum62.txt")
#########S = alignment_util.readScoringMatrix("/Users/ashwanikumar/Documents/project_5/PAM250.txt")
######S = alignment_util.readScoringMatrix("/Users/ashwanikumar/Documents/project_5/DNA1.txt")
#########S = alignment_util.readScoringMatrix("/Users/ashwanikumar/Documents/project_5/DNA2.txt")
#########print NeedlemanWunsch("MNALSDRT", "MGSDRTTET", S, 12)
######print NeedlemanWunsch("AATT", "AACTT",S,1)
########
def SmithWaterman(seq1, seq2, S, gap):
    """
   Input:
    * seq1, seq2: Two *strings* (not dnaSeqs) to be aligned.
      -- May not necessarily represent a dna sequence.
    * S: a scoring dictionary such that for any characters a,b: S[a][b] is the miss-match penalty for aligning
         an a to a b.
      -- You may assume that every character contained in seq1 or seq2 is defined in S.
    * gap: The gap penalty, listed as a non-negative number.  (e.g. If gap=4, we will deduct 4 points for each -.)
    Output: A tuple containing three elements: 
    1) The optimal local alignment score.
    2) The corresponding alignment string using a substring of seq1.  (That is, a copy of a substring of seq1 string with "-" characters added corresponding to an optimal alignment.)
    3) The corresponding alignment string using a substring of seq2.  (That is, a copy of a substring of seq2 string with "-" characters added corresponding to an optimal alignment.)
    """
    assert gap > 0, "Gap penalties should be a positive number"
    seq1 = " " + seq1
    seq2 = " " + seq2
    n, m = len(seq1), len(seq2)

    M = [[ float(0*j) for j in range(0, m) ]]    # Score matrix
    P = [[0] + [ 2 for i in range(1,m) ]]   # Pointer matrix -- record whether the score vame from the top-left (0), top (1), or left (2).
    maximum_value = 0
    backtracing_start  = 0
    for i in range(1, n):
        M.append([0]*m)
        P.append([-1]*m)
        M[i][0] = float(0*i)
        P[i][0] = 1
        for j in range(1,m):
            t1 = M[i-1][j-1] + S[seq1[i]][seq2[j]]
            t2 = M[i-1][j] - gap
            t3 = M[i][j-1] - gap
            M[i][j], P[i][j] = max((t1,0), (t2,1), (t3,2),(0,-1))  # max will use the tuple with the largest first element.
            if M[i][j] > maximum_value:
                i_max, j_max = 0,0
                maximum_value = M[i][j]
                i_max = i
                j_max = j
                backtracing = M[i][j]
    i, j = i_max, j_max
    s1, s2 = "", ""
    while M[i][j] > 0:
        if P[i][j] == 0:
            s1 += seq1[i]
            s2 += seq2[j]
            i -= 1
            j -= 1
        elif P[i][j] == 1:
            s1 += seq1[i]
            s2 += '-'
            i -= 1
        elif P[i][j] == 2:
            s1 += '-'
            s2 += seq2[j]
            j -= 1
        else:
            assert(False)

    return M[i_max][j_max], s1[::-1], s2[::-1]
#######S = alignment_util.readScoringMatrix("/Users/ashwanikumar/Documents/project_5/Blosum62.txt")
#####S = alignment_util.readScoringMatrix("/Users/ashwanikumar/Documents/project_5/PAM250.txt")
#####S = alignment_util.readScoringMatrix("/Users/ashwanikumar/Documents/project_5/DNA1.txt")
#####S = alignment_util.readScoringMatrix("/Users/ashwanikumar/Documents/project_5/DNA2.txt")
#####print SmithWaterman("AATT", "AACTT",S,1)
#####print SmithWaterman("MNALSDRT", "MGSDRTTET", S, 12)


def SmithWatermanAffine(seq1, seq2, S, o, c):
    """
   Input:
    * seq1, seq2: Two *strings* (not dnaSeqs) to be aligned.
      -- May not necessarily represent a dna sequence.
    * S: a scoring dictionary such that for any characters a,b: S[a][b] is the miss-match penalty for aligning
         an a to a b.
      -- You may assume that every character contained in seq1 or seq2 is defined in S.
    * o: The gap opening penalty: the penalty assessed for each segment of gaps in a sequence.
    * c: The gap continuation penalty: the individual penalty for each gap symbol in the alignment.
    Example: Given sequences AAATTTAAAAAA and AAAAAAGGGGAAA, with S[A][A] = 10, S[C][C] = 10, S[a][b] = 0 for a!=b, o = 5, and c=1, the
    optimal alignment would be:
               AAATTTAAA----AAA
               AAA---AAAGGGGAAA
    with a score of 30 points for the matches, -5 for each of the initial gaps, and -1 for each gap continuation.
    30 - 5 - 3 + 30 - 5 - 4 + 30 = 73.

    Output: A tuple containing three elements: 
    1) The optimal (affgine) local alignment score.
    2) The corresponding alignment string using a substring of seq1.  (That is, a copy of a substring of seq1 string with "-" characters added corresponding to an optimal alignment.)
    3) The corresponding alignment string using a substring of seq2.  (That is, a copy of a substring of seq2 string with "-" characters added corresponding to an optimal alignment.)

    """
    assert o > 0 and c > 0, "Gap penalties should be a positive numbers"
    m = len(seq1)+1
    n = len(seq2)+1
    M = np.zeros((m, n))
    P = np.zeros((m, n))
    Q = np.zeros((m, n))
    pointer_M = np.zeros((m,n))
    pointer_P = np.zeros((m,n))
    pointer_Q = np.zeros((m,n))
    pointer_MH = np.zeros((m,n))
    pointer_MU = np.zeros((m,n))
    pointer_QU = np.zeros((m,n))
    pointer_PH = np.zeros((m,n))
    
    for k in range(0,m):
        pointer_M[k][0] = 1
        pointer_P[k][0] = 1
        pointer_Q[k][0] = 1
        pointer_MH[k][0] = 1
        pointer_MU[k][0] = 1
        pointer_QU[k][0] = 1
        pointer_PH[k][0] = 1
        
    for z in range(0,n):
        pointer_M[0][z] = 2
        pointer_P[0][z] = 2
        pointer_Q[0][z] = 2
        pointer_MU[0][z] = 2
        pointer_MH[k][0] = 2
        pointer_QU[k][0] = 2
        pointer_PH[k][0] = 2
        

    pointer_M[0][0] = 0
    pointer_P[0][0] = 0
    pointer_Q[0][0] = 0
    pointer_MU[0][z] = 2
    pointer_MH[k][0] = 2
    pointer_QU[k][0] = 0
    pointer_PH[k][0] = 0
    
    for i in range(0,m):
        M[i][0] = 0
        P[i][0] = 0
        Q[i][0] = 0
    for j in range(0,n):
        M[0][j] = 0
        P[0][j] = 0
        Q[0][j] = 0
    max_score = 0
    i_max = 0
    j_max = 0
    M_score = 0
    P_score = 0
    Q_score = 0
    for i in range(1,m):
        for j in range(1,n):            
            ##calculation of maximum score for P[i][j] (This matrix handles gap in seq 1) and assignment of pointers##
 
            MP_GOP = M[i][j-1] - o - c
            P_Extension = P[i][j-1] - c
            P[i][j] = max(MP_GOP, P_Extension,0)
            if P[i][j] == P_Extension:
                pointer_PH[i][j-1] = 1
            if P[i][j] == MP_GOP:
                pointer_MH[i][j-1] = 1
           
            ##calculation of maximum score for Q[i][j] (This matrix handles gap in seq 2) and assignment of pointers ###

            MQ_GOP = M[i-1][j] - o - c
            Q_Extension = Q[i-1][j] - c
            Q[i][j] = max(MQ_GOP,Q_Extension,0)
            if Q[i][j] == Q_Extension:
                pointer_QU[i-1][j] = 2
            if Q[i][j] == MQ_GOP:
                pointer_MU[i-1][j] = 2

                ###calculation of maximum score for M[i][j] and the assignment of pointers##

            M_score = S[seq1[i-1]][seq2[j-1]] + M[i-1][j-1]
            P_score = S[seq1[i-1]][seq2[j-1]] + P[i-1][j-1]
            Q_score = S[seq1[i-1]][seq2[j-1]] + Q[i-1][j-1]
            M[i][j] = max(M_score,P_score,Q_score,0)
            if M[i][j] == M_score: 
                pointer_M[i-1][j-1] = 6
            elif M[i][j] == P_score:
                pointer_P[i-1][j-1] = 6
            elif M[i][j] == Q_score: 
                pointer_Q[i-1][j-1] = 6
            
            if M[i][j] >= max_score:
                i_max, j_max = 0,0
                max_score = M[i][j]
                i_max = i
                j_max = j
    i,j = i_max,j_max


##    print M
##    print P
##    print Q
##    print pointer_M
##    print pointer_P
##    print pointer_Q
##    print pointer_MU
##    print pointer_MH
    s1, s2 = "", ""
    #print i,j
    while M[i][j] > 0:
        
        if pointer_M[i-1][j-1] == 6:
            s1 += seq1[i-1]
            s2 += seq2[j-1]
            i -= 1
            j -= 1
            
        elif pointer_Q[i-1][j-1] == 6:
            s1 += seq1[i-1]
            s2 += seq2[j-1]
            i -= 1
            j -= 1
            while pointer_QU[i-1][j] == 2:
                s2 += '-'
                s1 += seq1[i]
                i -= 1               
            else:
                if pointer_MU[i-1][j] == 2:
                    s1 += seq1[i]
                    s2 += '-'
                    i -= 1
                          
        elif pointer_P[i-1][j-1] == 6:
            s1 += seq1[i-1]
            s2 += seq2[j-1]
            i -= 1
            j -= 1
            while pointer_PH[i][j-1] == 1:
                s2 += seq2[j]
                s1 += '-'
                j -= 1
            else:
                if pointer_MH[i][j-1] == 1:
                    s1 += '-'
                    s2 += seq2[j]
                    j -= 1
        else:
            print "I am done"
                 
    return M[i_max][j_max], s1[::-1], s2[::-1]

                        
#S = alignment_util.readScoringMatrix("/Users/ashwanikumar/Documents/project_5/Blosum62.txt")
#S = alignment_util.readScoringMatrix("/Users/ashwanikumar/Documents/project_5/DNA2.txt")
#S = {'A': {'A': 10, 'C': -10, 'T': -10, 'G': -10}, 'C': {'A': -10, 'C': 10, 'T': -10, 'G': -10}, 'T': {'A': -10, 'C': -10, 'T': 10, 'G': -10}, 'G': {'A': -10, 'C': -10, 'T': -10, 'G': 10}}
#A = "MERPEPELIRQSWRAVSRSPLEHGTVLFARLFALEPDLLPLFQYNCRQFSSPEDCLSSPEFLDHIRKVMLVIDAAVTNVEDLSSLEEYLASLGRKHRAVGVKLSSFSTVGESLLYMLEKCLGPAFTPATRAAWSQLYGAVVQAMSRGWDGE"
#B = "MGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKSEDEMKASEDLKKHGATVLTALGGILKKKGHHEAEIKPLAQSHATKHKIPVKYLEFISECIIQVLQSKHPGDFGADAQGAMNKALELFRKDMASNYKELGFQG"
#print SmithWatermanAffine("AAAAAGTGAAAAA","AAAAACAAAAA" ,S, 1, 1)
#print SmithWatermanAffine("AATTA","AAA" ,S, 4, 1)
#print SmithWatermanAffine("MQRGGDEQ","MQRGGGGGDEQ" ,S, 11, 1)
#print SmithWatermanAffine("MPRSFLVRKPSDPRRKPNYSELQDACVEFTFQQPYDQAHLLAAIPPPEVLNPAASLPTLIWDSLLVPQVRPVAWATLPLRESPKAVELTSLSDEDSGKSSQPPSPPSPASSFSSTSASS","MPRSFLVRKPSDPRRKPNYSELQEVLNPAASLPTLIWDSLLVPQVRPVAWATLPLRESPKAVELTSLSDEDSGKSSQPPSPPSPAPSSFSSTSASS",S, 11, 1)
#A = "MADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTIDFPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVDEMIREADIDGDGQVNYEEFVQMMTAK"
#B = "MSEEIITPVYCTGVSAQVQKQRARELGLGRHENAIKYLGQDYEQLRVRCLQSGTLFRDEAFPPVPQSLGYKDLGPNSSKTYGIKWKRPTELLSNPQFIVDGATRTDICQGALGDCWLLAAIASLTLNDTLLHRVVPHGQSFQNGYAGIFHFQLWQFGEWVDVVVDDLLPIKDGKLVFVHSAEGNEFWSALLEKAYAKVNGSYEALSGGSTSEGFEDFTGGVTEWYELRKAPSDLYQIILKALERGSLLGCSIDISSVLDMEAITFKKLVKGHAYSVTGAKQVNYRGQVVSLIRMRNPWGEVEWTGAWSDSSSEWNNVDPYERDQLRVKMEDGEFWMSFRDFMREFTRLEICNLTPDALKSRTIRKWNTTLYEGTWRRGSTAGGCRNYPATFWVNPQFKIRLDETDDPDDYGDRESGCSFVLALMQKHRRRERRFGRDMETIGFAVYEVPPELVGQPAVHLKRDFFLANASRARSEQFINLREVSTRFRLPPGEYVVVPSTFEPNKEGDFVLRFFSEKSAGTVELDDQIQANLPDEQVLSEEEIDENFKALFRQLAGEDMEISVKELRTILNRIISKHKDLRTKGFSLESCRSMVNLMDRDGNGKLGLVEFNILWNRIRNYLSIFRKFDLDKSGSMSAYEMRMAIESAGFKLNKKLYELIITRYSEPDLAVDFDNFVCCLVRLETMFRFFKTLDTDLDGVVTFDLFKWLQLTMFA"
#print SmithWatermanAffine(A,B ,S, 11, 1)
   










    
