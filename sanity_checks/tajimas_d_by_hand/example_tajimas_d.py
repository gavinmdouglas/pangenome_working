#!/usr/bin/python3

import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, '/home/gdouglas/scripts/pangenome_working')

from functions.pop_gen import (num_pairwise_diff_bases,
                               tajimas_d_and_wattersons_theta)

from statistics import mean

num_diff1 = [0, 5, 5, 5, 5, 5, 5, 0, 0, 0]
num_diff2 = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2]


print(mean(num_diff1))
print(tajimas_d_and_wattersons_theta(mean(num_diff1), 5, 5))

print("\n\n\n")


print(mean(num_diff2))
print(tajimas_d_and_wattersons_theta(mean(num_diff2), 5, 5))