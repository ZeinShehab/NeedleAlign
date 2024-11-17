# python3

import sys
import numpy as np
from parser import Parser

INF = 1000000

class NeedleSolver:
    def __init__(self):
        self.read_input_files()
        
        self.scoring_matrix_file = "ednafull.txt"
        self.gap_opening_penalty = 10
        self.gap_extension_penalty = 0.5
    
        self.find_best_alignment()

            
    def find_best_alignment(self):
        print("[*] Beginning allignment...\n")
        
        n = len(self.seq_list1)
        
        for seq1_idx in range(n):
            seq1 = self.seq_list1[seq1_idx]
            
            best_score = 0
            best_seq2 = ""
            best_seq1_idx = -1
            best_seq2_idx = -1
            results = ()
            
            m = len(self.seq_list2)
            for seq2_idx in range(m):
                seq2 = self.seq_list2[seq2_idx]
                
                score, backtrack, i, j = self.LCSBackTrack(seq1, seq2)
                
                if score > best_score:
                    best_score = score
                    results = (backtrack, i, j)
                    best_seq2 = seq2
                    best_seq2_idx = seq2_idx
                    
            backtrack, i, j = results
            s1, s2 = self.OutputLCS(backtrack, seq1, best_seq2, i, j)
            
            print(f"Best score: {best_score}")
            print(f"Line: {seq1_idx:<5} | {s1}")
            print(f"Line: {best_seq2_idx:<5} | {s2}")            
            print("\n")

    def read_input_files(self):
        if len(sys.argv) < 3:
            print("Expected 2 arguments got 1. Please enter name of files")
            sys.exit(1)
        
        parser = Parser()
    
        temp_seq_list1 = parser.parse_file(sys.argv[1])
        temp_seq_list2 = parser.parse_file(sys.argv[2])
        
        if len(temp_seq_list1[0]) > len(temp_seq_list2[0]):
            self.seq_list1 = temp_seq_list1
            self.seq_list2 = temp_seq_list2
        else:
            self.seq_list1 = temp_seq_list2
            self.seq_list2 = temp_seq_list1
    
    def scoring_matrix(self):
        sMatrixTxt = open(self.scoring_matrix_file, "r").read()
        sMatrixList = sMatrixTxt.strip().split('\n')
        aaList = sMatrixList[0].split()
        sMatrix = dict()
        
        for aa in aaList:
            sMatrix[aa] = dict()
        
        for i in range(1, len(aaList) + 1):
            currRow = sMatrixList[i].split()
            for j in range(len(aaList)):
                sMatrix[currRow[0]][aaList[j]] = int(currRow[j + 1])
        
        return sMatrix

    def LCSBackTrack(self, v, w):
        sMatrix = self.scoring_matrix()
        n = len(v)
        m = len(w)
        sLower = np.matrix(-INF * np.ones((n+1)*(m+1), dtype = np.int64).reshape((n+1, m+1)))
        sMiddle = np.matrix(-INF * np.ones((n+1)*(m+1), dtype = np.int64).reshape((n+1, m+1)))
        sUpper = np.matrix(-INF * np.ones((n+1)*(m+1), dtype = np.int64).reshape((n+1, m+1)))
        s = [sLower, sMiddle, sUpper]
        backtrackLower = np.matrix(np.zeros((n+1)*(m+1), dtype = np.int64).reshape((n+1, m+1)))
        backtrackMiddle = np.matrix(np.zeros((n+1)*(m+1), dtype = np.int64).reshape((n+1, m+1)))
        backtrackUpper = np.matrix(np.zeros((n+1)*(m+1), dtype = np.int64).reshape((n+1, m+1)))
        backtrack = [backtrackLower, backtrackMiddle, backtrackUpper]
        s[0][0, 0] = 0
        s[1][0, 0] = 0
        s[2][0, 0] = 0
        s[0][1, 0] = -self.gap_opening_penalty
        s[1][1, 0] = -self.gap_opening_penalty
        for i in range(2, n+1):
            s[0][i, 0] = s[0][i-1, 0] - self.gap_extension_penalty
            s[1][i, 0] = s[0][i, 0]
        s[2][0, 1] = -self.gap_opening_penalty
        s[1][0, 1] = -self.gap_opening_penalty
        for j in range(2, m+1):
            s[2][0, j] = s[2][0, j-1] - self.gap_extension_penalty
            s[1][0, j] = s[2][0, j]
        for i in range(1, n+1):
            for j in range(1, m+1):
                score1 = s[0][i-1, j] - self.gap_extension_penalty
                score2 = s[1][i-1, j] - self.gap_opening_penalty
                lowerScore = max(score1, score2)
                s[0][i, j] = lowerScore
                if lowerScore == score2:
                    backtrack[0][i, j] = 1 # From middle
                score1 = s[2][i, j-1] - self.gap_extension_penalty
                score2 = s[1][i, j-1] - self.gap_opening_penalty
                upperScore = max(score1, score2)
                s[2][i, j] = upperScore
                if upperScore == score2:
                    backtrack[2][i, j] = 1 # From middle
                score3 = s[1][i-1, j-1] + sMatrix[v[i-1]][w[j-1]]
                middleScore = max(lowerScore, score3, upperScore)
                s[1][i, j] = middleScore
                if middleScore == score3:
                    backtrack[1][i, j] = 1 # From middle
                elif middleScore == upperScore:
                    backtrack[1][i, j] = 2 # From upper
        return s[1][n, m], backtrack, n, m                       


    def OutputLCS(self, backtrack, v, w, i, j):
        s1 = ''
        s2 = ''
        level = 1 # middle
        while i > 0 and j > 0:
            if 0 == level:
                if 1 == backtrack[0][i, j]:
                    level = 1
                i -= 1
                s1 = v[i] + s1
                s2 = '-' + s2
                continue
            elif 2 == level:
                if 1 == backtrack[2][i, j]:
                    level = 1
                j -= 1
                s1 = '-' + s1
                s2 = w[j] + s2
                continue
            else:
                if 1 == backtrack[1][i, j]:
                    i -= 1
                    j -= 1
                    # try:
                        # if (i)
                        
                    s1 = v[i] + s1
                    s2 = w[j] + s2
                    continue
                    # except IndexError:
                    #     print(f"i={i} | j={j}")
                    #     print(f"v={v}\nw={w}")
                    #     sys.exit(1)
                else:
                    level = backtrack[1][i, j]
        return s1, s2
    
    def computeCLS(self, seq1, seq2):
        maxSocre, backtrack, i, j = self.LCSBackTrack(seq1, seq2)
        s1, s2 = self.OutputLCS(backtrack, seq1, seq2, i, j)
        return maxSocre, s1, s2

if __name__ == "__main__":
    NeedleSolver()