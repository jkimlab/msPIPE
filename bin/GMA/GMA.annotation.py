#! /usr/bin/env python3

import sys
import re
'''
import subprocess as sub
import os

abs_path = os.path.split(os.path.abspath(sys.argv[0]))[0]
mspipeD = "/".join(abs_path.split("/")[:-2])
libD = mspipeD + "/lib"

sys.path.append(libD)

from func import *
from methClass import *
'''
try:
    gtfF = sys.argv[1]
    outD = sys.argv[2] if sys.argv[2][-1] != "/" else sys.argv[2][:-1]
except:
    print( f'GMA.annotation.py [gtf file] [out dir] ' ,file=sys.stderr)
    sys.exit()


def field(line) :
    col = line.strip(";").split(";")
    dic = {}
    for e in col:
        e = e.strip()
        E = e.split(' ')
        dic[E[0].strip()] = E[1].replace('"',"")
    return dic


exonF = open(outD + "/exon.bed", 'w')
geneF = open(outD + "/gene.bed", 'w')


for l in open(gtfF , 'r'):
    e = l.rstrip().split("\t")
    ch = e[0]
    start = int(e[3]) -1
    end = e[4]

    if l[0] == "#" : continue
    opt = field("".join(e[8:]))
    # filtering
    if e[2] == 'transcript':
        gene_id = opt['gene_id']

        gene_line= [ch, str(start), end,gene_id, e[5], e[6], e[1],e[2],e[7], e[8]]
        geneF.write('\t'.join(gene_line) + '\n')
    
    elif e[2] == 'exon':
        gene_id = opt['gene_id']
        exon_line = [ch, str(start), end,gene_id, e[5], e[6], e[1],e[2],e[7], e[8]]
        exonF.write('\t'.join(exon_line) + '\n')

'''
for ch in sorted(geneD):
    for start in sorted(geneD[ch]):
        geneF.write('\t'.join(geneD[ch][start]) + '\n')

for ch in sorted(exonD):
    for start in sorted(exonD[ch]):
        exonF.write('\t'.join(exonD[ch][start]) + '\n')
'''
geneF.close()
exonF.close()


        











