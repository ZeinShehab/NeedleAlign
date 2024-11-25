import numpy as np
from matrix import initialize_matrices


# Define the Needleman-Wunsch alignment function
def needleman_wunsch(seq1, seq2, match_score, mismatch_score, gap_open, gap_extend):
    # Calculate the dimensions of the matrix based on the sequence lengths
    rows, cols = len(seq1) + 1, len(seq2) + 1
    # Initialize the score matrices for alignment
    S, Ix, Iy = initialize_matrices(rows, cols, gap_open, gap_extend)
    # Fill the score matrices with scores for matching, mismatching, and gaps
    for i in range(1, rows):
        for j in range(1, cols):
            # Calculate match or mismatch score
            match = S[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score)
            # Calculate the score for inserting a gap in seq1
            Ix[i][j] = max(Ix[i - 1][j] + gap_extend, S[i - 1][j] + gap_open)
            # Calculate the score for inserting a gap in seq2
            Iy[i][j] = max(Iy[i][j - 1] + gap_extend, S[i][j - 1] + gap_open)
            # Choose the maximum score from match, insert gap in seq1, or insert gap in seq2
            S[i][j] = max(match, Ix[i][j], Iy[i][j])
    # Return the filled matrices and the final alignment score from the bottom-right corner
    return S, Ix, Iy, S[rows - 1][cols - 1]


# Define the traceback function to find the best alignment path
def trace_back(S, seq1, seq2, match_score, mismatch_score, gap_extend):
    alignment_a, alignment_b = '', ''
    i, j = len(seq1), len(seq2)
    # Trace back from the bottom-right corner of the matrix to the top-left
    while i > 0 and j > 0:
        # Check if the current cell was filled from a match or mismatch
        if S[i][j] == S[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score):
            alignment_a = seq1[i - 1] + alignment_a
            alignment_b = seq2[j - 1] + alignment_b
            i -= 1
            j -= 1
        # Check if the current cell was filled from inserting a gap in seq2
        elif S[i][j] == S[i][j - 1] + gap_extend:
            alignment_a = '-' + alignment_a
            alignment_b = seq2[j - 1] + alignment_b
            j -= 1
        # Otherwise, it must have been filled from inserting a gap in seq1
        else:
            alignment_a = seq1[i - 1] + alignment_a
            alignment_b = '-' + alignment_b
            i -= 1
    # Fill in the remaining gaps if the beginning of the matrix has not been reached
    while i > 0:
        alignment_a = seq1[i - 1] + alignment_a
        alignment_b = '-' + alignment_b
        i -= 1
    while j > 0:
        alignment_a = '-' + alignment_a
        alignment_b = seq2[j - 1] + alignment_b
        j -= 1
    # Return the final aligned sequences
    return alignment_a, alignment_b
