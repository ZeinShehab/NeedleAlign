#!/usr/bin/env python
"""
The Needleman-Wunsch Algorithm
==============================

This is a dynamic programming algorithm for finding the optimal alignment of
two strings.

Example
-------

    >>> x = "GATTACA"
    >>> y = "GCATGCU"
    >>> print(nw(x, y))
    G-ATTACA
    GCA-TGCU


LICENSE

This is free and unencumbered software released into the public domain.
Anyone is free to copy, modify, publish, use, compile, sell, or
distribute this software, either in source code form or as a compiled
binary, for any purpose, commercial or non-commercial, and by any
means.

In jurisdictions that recognize copyright laws, the author or authors
of this software dedicate any and all copyright interest in the
software to the public domain. We make this dedication for the benefit
of the public at large and to the detriment of our heirs and
successors. We intend this dedication to be an overt act of
relinquishment in perpetuity of all present and future rights to this
software under copyright law.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

For more information, please refer to <http://unlicense.org/>
"""

import numpy as np

from parser import Parser

def scoring_matrix():
    sMatrixTxt = open("ednafull.txt", "r").read()
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


def nw(x, y, scoring_matrix, gap_opening=10, gap_extension=0.5):
    nx = len(x)
    ny = len(y)

    # Optimal score at each possible pair of characters.
    F = np.zeros((nx + 1, ny + 1))
    G_x = np.zeros((nx + 1, ny + 1))  # Tracks gaps in x (row gaps)
    G_y = np.zeros((nx + 1, ny + 1))  # Tracks gaps in y (column gaps)

    # Initialize first row and column with gap penalties
    for i in range(1, nx + 1):
        F[i, 0] = -gap_opening - (i - 1) * gap_extension
    for j in range(1, ny + 1):
        F[0, j] = -gap_opening - (j - 1) * gap_extension

    # Pointers to trace through an optimal alignment
    P = np.zeros((nx + 1, ny + 1))
    P[:, 0] = 3  # Vertical gaps
    P[0, :] = 4  # Horizontal gaps

    # Temporary scores
    t = np.zeros(3)

    # Fill in the scoring matrix
    for i in range(1, nx + 1):
        for j in range(1, ny + 1):
            # Match or mismatch score from scoring matrix
            t[0] = F[i - 1, j - 1] + scoring_matrix[x[i - 1]][y[j - 1]]
            # Gap penalties
            t[1] = max(F[i - 1, j] - gap_opening, G_x[i - 1, j] - gap_extension)  # Vertical gap
            t[2] = max(F[i, j - 1] - gap_opening, G_y[i, j - 1] - gap_extension)  # Horizontal gap

            # Update F, G_x, G_y with max scores
            F[i, j] = max(t)
            G_x[i, j] = t[1]
            G_y[i, j] = t[2]

            # Update pointer matrix
            if t[0] == F[i, j]:
                P[i, j] += 2  # Diagonal (match/mismatch)
            if t[1] == F[i, j]:
                P[i, j] += 3  # Vertical gap
            if t[2] == F[i, j]:
                P[i, j] += 4  # Horizontal gap

    # Trace through an optimal alignment
    i = nx
    j = ny
    rx = []
    ry = []
    while i > 0 or j > 0:
        if P[i, j] in [2, 5, 6, 9]:
            rx.append(x[i - 1])
            ry.append(y[j - 1])
            i -= 1
            j -= 1
        elif P[i, j] in [3, 5, 7, 9]:
            rx.append(x[i - 1])
            ry.append('-')
            i -= 1
        elif P[i, j] in [4, 6, 7, 9]:
            rx.append('-')
            ry.append(y[j - 1])
            j -= 1

    # Reverse the strings
    rx = ''.join(rx)[::-1]
    ry = ''.join(ry)[::-1]
    return rx, ry
    # return '\n'.join([rx, ry])


# def nw(x, y, match = 1, mismatch = 1, gap = 1):
#     nx = len(x)
#     ny = len(y)
#     # Optimal score at each possible pair of characters.
#     F = np.zeros((nx + 1, ny + 1))
#     F[:,0] = np.linspace(0, -nx * gap, nx + 1)
#     F[0,:] = np.linspace(0, -ny * gap, ny + 1)
#     # Pointers to trace through an optimal aligment.
#     P = np.zeros((nx + 1, ny + 1))
#     P[:,0] = 3
#     P[0,:] = 4
#     # Temporary scores.
#     t = np.zeros(3)
#     for i in range(nx):
#         for j in range(ny):
#             if x[i] == y[j]:
#                 t[0] = F[i,j] + match
#             else:
#                 t[0] = F[i,j] - mismatch
#             t[1] = F[i,j+1] - gap
#             t[2] = F[i+1,j] - gap
#             tmax = np.max(t)
#             F[i+1,j+1] = tmax
#             if t[0] == tmax:
#                 P[i+1,j+1] += 2
#             if t[1] == tmax:
#                 P[i+1,j+1] += 3
#             if t[2] == tmax:
#                 P[i+1,j+1] += 4
#     # Trace through an optimal alignment.
#     i = nx
#     j = ny
#     rx = []
#     ry = []
#     while i > 0 or j > 0:
#         if P[i,j] in [2, 5, 6, 9]:
#             rx.append(x[i-1])
#             ry.append(y[j-1])
#             i -= 1
#             j -= 1
#         elif P[i,j] in [3, 5, 7, 9]:
#             rx.append(x[i-1])
#             ry.append('-')
#             i -= 1
#         elif P[i,j] in [4, 6, 7, 9]:
#             rx.append('-')
#             ry.append(y[j-1])
#             j -= 1
#     # Reverse the strings.
#     rx = ''.join(rx)[::-1]
#     ry = ''.join(ry)[::-1]
#     return '\n'.join([rx, ry])

scoring_matrix = scoring_matrix()

parser = Parser()
seq2 = parser.parse_nosplit("input1.txt")
seq1 = parser.parse_nosplit("input2.txt")

names1 = seq1.keys()
names2 = seq2.keys()

# for name1 in names1:
#     for name2 in names2:
#         x = seq1[name1]
#         y = seq2[name2]

#         s1, s2 = nw(x, y, scoring_matrix)

#         for i in range(0, len(s1), 50):
#             print(s1[i:i+50])
#             print(s2[i:i+50])
#             print("\n")

# x = "ATATTAGGTTTTTACCTACCCAGGAAAAGCCAACCAACCTCGATCTCTTGTAGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTAGCTGTCGCTCGGCTGCATGCCTAGTGCACCTACGCAGTATAAACAATAATAAATTTTACTGTC"
# y = "AGCAGTA"


print(nw('GATTACA', 'GCATGCT', scoring_matrix))