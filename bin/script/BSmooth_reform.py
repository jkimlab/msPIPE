#! /usr/bin/env python3

import sys


inF = sys.argv[1]


for line in open(inF, 'r'):
    col = line.strip().split("\t")
    new_line = col[:3] + [col[15]] + [col[6]] + [col[14]]
    print("\t".join(new_line))
    

