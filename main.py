# python3

import sys
import numpy as np
from parser import Parser
from tqdm import tqdm

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

        node_names_1 = self.node_list1.keys()
        node_names_2 = self.node_list2.keys()

        for node1 in node_names_1:
            for node2 in node_names_2:
                seq1 = self.node_list1[node1]
                seq2 = self.node_list2[node2]

                best_score, s1, s2, gaps, matches, mismatches = self.computeCLS(seq1, seq2)

                max_seq_len = gaps + matches + mismatches
                max_name_len = int(max(len(node1), len(node2)))

                gaps_percentage = (gaps/max_seq_len) * 100
                match_percentage = (matches/max_seq_len) * 100
                mismatch_percentage = (mismatches/max_seq_len) * 100

                print(f"Best score: {best_score}")
                print(f"Gaps: {gaps} | {gaps_percentage:.2f} %")
                print(f"Matches: {matches} | {match_percentage:.2f} %")
                print(f"Mismatches: {mismatches} | {mismatch_percentage:.2f} %")

                self.write_output("output.txt", node1, node2, s1, s2)

    def write_output(self, filename, node1, node2, s1, s2):
        with open(filename, 'a+') as file:
            max_name_len = int(max(len(node1), len(node2)))
            align_string = ""

            for i in range(len(s1)):
                if s1[i] == s2[i]:
                    align_string += "|"
                elif s1[i] == "-" or s2[i] == "-":
                    align_string += " "
                else:
                    align_string += "."

            output_s1 = [s1[i:i+50] for i in range(0, len(s1), 50)]
            output_s2 = [s2[i:i+50] for i in range(0, len(s2), 50)]
            output_align = [align_string[i:i+50] for i in range(0, len(align_string), 50)]

            padding = " " * (max_name_len + 1)

            for i in range(len(output_s1)):
                file.write(f"{node1:<{max_name_len}}  {output_s1[i]:<51}\n")
                file.write(f"{padding:<{len(padding)}} {output_align[i]:<51}\n")
                file.write(f"{node2:<{max_name_len}}  {output_s2[i]:<51}\n")
                file.write("\n")
                file.write("\n")

    def read_input_files(self):
        if len(sys.argv) < 3:
            print("Expected 2 arguments got 1. Please enter name of files")
            sys.exit(1)

        parser = Parser()

        self.node_list1 = parser.parse_nosplit(sys.argv[1])
        self.node_list2 = parser.parse_nosplit(sys.argv[2])

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
        shape = (n + 1, m + 1)

        s = [np.full(shape, -INF, dtype=np.int64) for _ in range(3)]
        backtrack = [np.zeros(shape, dtype=np.int64) for _ in range(3)]

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

        for i in tqdm(range(1, n+1)):
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
        gaps = 0
        matches = 0
        mismatches = 0
        level = 1 # middle
        while i > 0 and j > 0:
            if 0 == level:
                if 1 == backtrack[0][i, j]:
                    level = 1
                i -= 1
                s1 = v[i] + s1
                s2 = '-' + s2
                gaps += 1
                continue
            elif 2 == level:
                if 1 == backtrack[2][i, j]:
                    level = 1
                j -= 1
                s1 = '-' + s1
                s2 = w[j] + s2
                gaps += 1
                continue
            else:
                if 1 == backtrack[1][i, j]:
                    i -= 1
                    j -= 1
                    if v[i] == w[j]:
                        s1 = v[i] + s1
                        s2 = w[j] + s2
                        matches += 1
                    else:
                        s1 = v[i] + s1
                        s2 = w[j] + s2
                        mismatches += 1
                    continue
                else:
                    level = backtrack[1][i, j]
        return s1, s2, gaps, matches, mismatches

    def computeCLS(self, seq1, seq2):

        if len(seq1) > len(seq2):
            s = seq1
            t = seq2
        else:
            s = seq2
            t = seq1

        maxScore, backtrack, i, j = self.LCSBackTrack(s, t)
        s1, s2, gaps, matches, mismatches = self.OutputLCS(backtrack, s, t, i, j)
        return maxScore, s1, s2, gaps, matches, mismatches

if __name__ == "__main__":
    NeedleSolver()