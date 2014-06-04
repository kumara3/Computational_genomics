##def NeedlemanWunsch(seq1, seq2, S, gap):
##    """
##    Input:
##    * seq1, seq2: Two *strings* (not dnaSeqs) to be aligned.
##      -- May not necessarily represent a dna sequence.
##    * S: a scoring dictionary such that gor any characters a,b: S[a][b] is the miss-match penalty for aligning
##         an a to a b.
##      -- You may assume that every character contained in seq1 or seq2 is defined in S.
##    * gap: The gap score, listed as a non-negative number.  (e.g. If gap=4, we will deduct 4 points for each -.)
##    Output: A tuple containing four elements: 
##    1) The optimal global alignment score.
##    3) The corresponding alignment string of seq1.  (That is, a copy of seq1 string with "-" characters added corresponding to an optimal alignment.
##    5) The corresponding alignment string of seq2.
##    """
##    assert gap > 0, "Gap should be a positive number"
##    seq1 = " " + seq1
##    seq2 = " " + seq2
##    n, m = len(seq1), len(seq2)
##
##    M = [[ float(-gap*j) for j in range(0, m) ]]    # Score matrix
##    P = [[0] + [ 2 for i in range(1,m) ]]   # Pointer matrix -- record whether the score vame from the top-left (0), top (1), or left (2).
##    for i in range(1, n):
##        M.append([0]*m)
##        P.append([-1]*m)
##        M[i][0] = float(-gap*i)
##        P[i][0] = 1
##        for j in range(1,m):
##            t1 = M[i-1][j-1] + S[seq1[i]][seq2[j]]
##            t2 = M[i-1][j] - gap
##            t3 = M[i][j-1] - gap
##            M[i][j], P[i][j] = max((t1,0), (t2,1), (t3,2))  # max will use the tuple with the largest first element
##
##    i, j = n-1, m-1
##    s1, s2 = "", ""
##    while i > 0 or j > 0:
##        if P[i][j] == 0:
##            s1 += seq1[i]
##            s2 += seq2[j]
##            i -= 1
##            j -= 1
##        elif P[i][j] == 1:
##            s1 += seq1[i]
##            s2 += '-'
##            i -= 1
##        elif P[i][j] == 2:
##            s1 += '-'
##            s2 += seq2[j]
##            j -= 1
##        else:
##            assert(False)
##    
##    return M[n-1][m-1], s1[::-1], s2[::-1]
##
##def SmithWaterman(seq1, seq2, S, gap):
##    """
##   Input:
##    * seq1, seq2: Two *strings* (not dnaSeqs) to be aligned.
##      -- May not necessarily represent a dna sequence.
##    * S: a scoring dictionary such that gor any characters a,b: S[a][b] is the miss-match penalty for aligning
##         an a to a b.
##      -- You may assume that every character contained in seq1 or seq2 is defined in S.
##    * gap: The gap score, listed as a non-negative number.  (e.g. If gap=4, we will deduct 4 points for each -.)
##    Output: A tuple containing four elements: 
##    1) The optimal local alignment score.
##    3) The corresponding alignment string using a substring of seq1.  (That is, a copy of seq1 string with "-" characters added corresponding to an optimal alignment.)
##    5) The corresponding alignment string using a substring of seq2.  (That is, a copy of seq2 string with "-" characters added corresponding to an optimal alignment.)
##    """
##    assert gap > 0, "Gap should be a positive number"
##    seq1 = " " + seq1
##    seq2 = " " + seq2
##    n, m = len(seq1), len(seq2)
##
##    M = [[0]*m]    # Score matrix
##    P = [[ -1 for i in range(0,m) ]]   # Pointer matrix -- record whether the score vame from the top-left (0), top (1), or left (2).
##    max_index = (0,0)
##    max_score = 0
##    for i in range(1, n):
##        M.append([0]*m)
##        P.append([-1]*m)
##        P[i][0] = -1
##        for j in range(1,m):
##            t1 = M[i-1][j-1] + S[seq1[i]][seq2[j]]
##            t2 = M[i-1][j] - gap
##            t3 = M[i][j-1] - gap
##            M[i][j], P[i][j] = max((t1,0), (t2,1), (t3,2), (0,-1))  # max will use the tuple with the largest first element
##            if M[i][j] > max_score:
##                max_score = M[i][j]
##                max_index = (i,j)
##
##    i, j = max_index
##    s1, s2 = "", ""
##    while P[i][j] != -1:
##        if P[i][j] == 0:
##            s1 += seq1[i]
##            s2 += seq2[j]
##            i -= 1
##            j -= 1
##        elif P[i][j] == 1:
##            s1 += seq1[i]
##            s2 += '-'
##            i -= 1
##        elif P[i][j] == 2:
##            s1 += '-'
##            s2 += seq2[j]
##            j -= 1
##        else:
##            assert(False)
##    
##    return M[max_index[0]][max_index[1]], s1[::-1], s2[::-1]
##
##
##def NeedlemanWunschAffine(seq1, seq2, S, s, c):
##    """
##   Input:
##    * seq1, seq2: Two *strings* (not dnaSeqs) to be aligned.
##      -- May not necessarily represent a dna sequence.
##    * S: a scoring dictionary such that for any characters a,b: S[a][b] is the miss-match penalty for aligning
##         an a to a b.
##      -- You may assume that every character contained in seq1 or seq2 is defined in S.
##    * s: The gap *start* score: the first - in any continuis sequence of gaps will be penalized at this value.
##    * c: The gap *continuation* score: the i > 1 gap - in any continuous sequence of gaps will be penalized at this score.
##    Example: Given sequences AAATTTAAAAAA and AAAAAAGGGGAAA, with S[A][A] = 10, S[C][C] = 10, S[a][b] = 0 for a!=b, s = 5, and c=1, the
##    optimal alignment would be:
##               AAATTTAAA----AAA
##               AAA---AAAGGGGAAA
##    with a score of 30 points for the matches, -5 for each of the initial gaps, and -1 for each gap continuation.
##    30 - (5+1+1) + 30 - (5+1+1+1) + 30 = 75.
##
##    Output: A tuple containing four elements: 
##    1) The optimal global alignment score.
##    3) The corresponding alignment string seq1.  (That is, a copy of seq1 string with "-" characters added corresponding to an optimal alignment.)
##    5) The corresponding alignment string seq2.  (That is, a copy of seq2 string with "-" characters added corresponding to an optimal alignment.)
##    """
##    assert s > 0 and c > 0, "Gap scores should be a positive numbers"
##    min_val = min([v for row in S.values() for v in row.values()] + [-c,-s])
##    min_score = min_val*(max(len(seq1),len(seq2))+2)
##
##    seq1 = " " + seq1
##    seq2 = " " + seq2
##    n, m = len(seq1), len(seq2)
##
##    M = [ # This contains the three score-holding matrices
##         [[0] + [min_score]*(m-1)],    # "Match" matrix
##         [[min_score] + [-s - j*c for j in range(1,m)]],   # "Insert" matrix
##         [[min_score] + [-s - j*c for j in range(1,m)]]    # "Delete" matrix
##        ]
##
##    P = [ # This contains the three point-back matrics.  P[k][i][j] = (m,d), where a is matrix a b is direction
##          [[(None,None)]*m],
##          [[(None,None)] + [(0,2)] + [(1,2)]*(m-2)],
##          [[(None,None)] + [(0,2)] + [(2,2)]*(m-2)]        
##        ]
##    # P[k][i][j] = (a,b) where:
##    # * a is the matrix (0, 1, or 2)
##    # * b is the direction: 0 = diagonal, 1 = up, 2 = left
##
##    for i in range(1,n):
##        M[0].append([min_score] + [None]*(m-1))
##        M[1].append([-s - i*c] + [None]*(m-1))
##        M[2].append([-s - i*c] + [None]*(m-1))
##    
##        P[0].append([(None,None)] + [None]*(m-1))
##        P[1].append([(0 if i==1 else 1, 1)] + [None]*(m-1))
##        P[2].append([(0 if i==1 else 2, 1)] + [None]*(m-1))
##
##        for j in range(1,m):
##            match_score = S[seq1[i]][seq2[j]]
##
##            t1 = M[0][i-1][j-1] + match_score
##            t2 = M[1][i-1][j-1] + match_score
##            t3 = M[2][i-1][j-1] + match_score
##            M[0][i][j], P[0][i][j] = max( (t1, (0,0)), (t2, (1,0)), (t3, (2,0)) )
##
##            t1 = M[0][i-1][j] - s - c
##            t2 = M[1][i-1][j] - c
##            t3 = M[2][i-1][j] - s - c
##            M[1][i][j], P[1][i][j] = max( (t1, (0,1)), (t2, (1,1)), (t3, (2,1)) )
##
##            t1 = M[0][i][j-1] - s - c
##            t2 = M[1][i][j-1] - s - c
##            t3 = M[2][i][j-1] - c
##            M[2][i][j], P[2][i][j] = max( (t1, (0,2)), (t2, (1,2)), (t3, (2,2)) )
##
##    s1, s2 = "", ""
##    k, i, j = 0, n-1, m-1
##
##    while i > 0 or j > 0:
##        p,d = P[k][i][j]
##
##        if d == 0:
##            s1 += seq1[i]
##            s2 += seq2[j]
##            i -= 1
##            j -= 1
##        elif d == 1:
##            s1 += seq1[i]
##            s2 += '-'
##            i -= 1
##        else:   # d == 2
##            s1 += '-'
##            s2 += seq2[j]
##            j -= 1
##        k = p
##    return M[0][n-1][m-1], s1[::-1], s2[::-1]
##    
##from alignment_util import *
##S1 = readScoringMatrix("DNA1.txt")
##S2 = readScoringMatrix("DNA2.txt")
##S3 = {x:{y:10 if x == y else -10 for y in "ACGT"} for x in "ACGT"}
###M = NeedlemanWunschAffine("ATCG", "CG", S3, 10, 5)
###M = NeedlemanWunschAffine("ATCG", "CG", S3, 1, 1)
##    
        

def SmithWatermanAffine(seq1, seq2, S, s, c):
    """
   Input:
    * seq1, seq2: Two *strings* (not dnaSeqs) to be aligned.
      -- May not necessarily represent a dna sequence.
    * S: a scoring dictionary such that gor any characters a,b: S[a][b] is the miss-match penalty for aligning
         an a to a b.
      -- You may assume that every character contained in seq1 or seq2 is defined in S.
    * s: The gap *start* score: the first - in any continuis sequence of gaps will be penalized at this value.
    * c: The gap *continuation* score: the i > 1 gap - in any continuous sequence of gaps will be penalized at this score.
    Example: Given sequences AAATTTAAAAAA and AAAAAAGGGGAAA, with S[A][A] = 10, S[C][C] = 10, S[a][b] = 0 for a!=b, s = 5, and c=1, the
    optimal alignment would be:
               AAATTTAAA----AAA
               AAA---AAAGGGGAAA
    with a score of 30 points for the matches, -5 for each of the initial gaps, and -1 for each gap continuation.
    30 - (5+1+1) + 30 - (5+1+1+1) + 30 = 75.

    Output: A tuple containing four elements: 
    1) The optimal local alignment score.
    3) The corresponding alignment string using a substring of seq1.  (That is, a copy of seq1 string with "-" characters added corresponding to an optimal alignment.)
    5) The corresponding alignment string using a substring of seq2.  (That is, a copy of seq2 string with "-" characters added corresponding to an optimal alignment.)
    """
    assert s > 0 and c > 0, "Gap scores should be a positive numbers"

    seq1 = " " + seq1
    seq2 = " " + seq2
    n, m = len(seq1), len(seq2)

    M = [ # This contains the three score-holding matrices
         [[0]*m],    # "Match" matrix
         [[0]*m],   # "Insert" matrix
         [[0]*m]    # "Delete" matrix
        ]

    P = [ # This contains the three point-back matrics.  P[k][i][j] = (m,d), where a is matrix a b is direction
          [[(-1,-1)]*m],
          [[(-1,-1)]*m],
          [[(-1,-1)]*m]
        ]
    # P[k][i][j] = (a,b) where:
    # * a is the matrix (0, 1, or 2)
    # * b is the direction: 0 = diagonal, 1 = up, 2 = left, -1 = stop

    best_point = (-1, -1, -1)
    for i in range(1,n):
        M[0].append([0] + [None]*(m-1))
        M[1].append([0] + [None]*(m-1))
        M[2].append([0] + [None]*(m-1))
    
        P[0].append([(-1,-1)] + [None]*(m-1))
        P[1].append([(-1,-1)] + [None]*(m-1))
        P[2].append([(-1,-1)] + [None]*(m-1))

        for j in range(1,m):
            match_score = S[seq1[i]][seq2[j]]

            t1 = M[0][i-1][j-1] + match_score
            t2 = M[1][i-1][j-1] + match_score
            t3 = M[2][i-1][j-1] + match_score
            M[0][i][j], P[0][i][j] = max( (t1, (0,0)), (t2, (1,0)), (t3, (2,0)), (0,(-1,-1)) )
            best_point = max(best_point, (M[0][i][j], i, j))

            t1 = M[0][i-1][j] - s - c
            t2 = M[1][i-1][j] - c
            t3 = M[2][i-1][j] - s - c
            M[1][i][j], P[1][i][j] = max( (t1, (0,1)), (t2, (1,1)), (t3, (2,1)), (0,(-1,-1)) )

            t1 = M[0][i][j-1] - s - c
            t2 = M[1][i][j-1] - s - c
            t3 = M[2][i][j-1] - c
            M[2][i][j], P[2][i][j] = max( (t1, (0,2)), (t2, (1,2)), (t3, (2,2)), (0,(-1,-1)) )

    s1, s2 = "", ""

    k, i, j = (0,) + best_point[1:]
    while P[k][i][j][1] != -1:
        p,d = P[k][i][j]

        if d == 0:
            s1 += seq1[i]
            s2 += seq2[j]
            i -= 1
            j -= 1
        elif d == 1:
            s1 += seq1[i]
            s2 += '-'
            i -= 1
        else:   # d == 2
            s1 += '-'
            s2 += seq2[j]
            j -= 1
        k = p
    return best_point[0], s1[::-1], s2[::-1]
            

        
            
                     
                     

