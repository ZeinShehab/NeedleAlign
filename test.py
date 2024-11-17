from Bio import pairwise2
from Bio.pairwise2 import format_alignment




alignments = pairwise2.align.globalms(v, w, match_score, mismatch_penalty, gap_open_penalty, gap_extend_penalty)
for alignment in alignments:
    print(format_alignment(*alignment))
