#!/usr/bin/python3
seq1 = "CAATTTCCACAC"
import sys
if len(sys.argv) == 2:
    seq1 = sys.argv[1]

print(seq1)

seq2 = ""
mapper = {"A": "T",
            "C": "G",
            "G": "C",
            "T": "A",
            }

# for c in seq1:
#     seq2 += c.lower()
seq1 = seq1[::-1]
for c in seq1:
    seq2 += mapper[c]

print(seq2)
