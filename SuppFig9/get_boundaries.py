import os, sys
import numpy as np
import matplotlib.pyplot as plt

chezod_dict = {}

fin = open('final_plddt_scores/MFIB_CheZOD_final_pLDDT_scores.txt', 'r')
for line in fin:
    new = line.strip().split()
    if new[0] != 'UniProt_ID':
        if int(new[2]) == 0:
            uniprot = str(new[0])
            res = int(new[1])
            if uniprot in chezod_dict.keys():
                chezod_dict[uniprot].append(res)
            else:
                chezod_dict[uniprot] = []
                chezod_dict[uniprot].append(res)
fin.close()

# https://stackoverflow.com/questions/55211695/find-consecutive-integers-in-a-list
from itertools import groupby

fout = open('new_CheZOD_boundaries.txt', 'w')
for key, val in chezod_dict.items():

    # Enumerate and get differences between counter—integer pairs
    # Group by differences (consecutive integers have equal differences)  
    gb = groupby(enumerate(val), key=lambda x: x[0] - x[1])

    # Repack elements from each group into list
    all_groups = ([i[1] for i in g] for _, g in gb)

    # Filter out one element lists
    final = list(filter(lambda x: len(x) > 1, all_groups))
    for i,res_range in enumerate(final):
        fout.write('%s\t%i\t%i\n' % (key, final[i][0], final[i][-1]))
        
fout.close()
    