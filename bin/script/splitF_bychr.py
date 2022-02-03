#! /usr/bin/env python3

import sys

try :
    inF = sys.argv[1]
    outPATH = sys.argv[2]
except IndexError:
    print('     "program.py [input file] [path/prefixed name]"')
    sys.exit()


bfkey = ""
i =0

for l in open(inF, 'r'):
    e = l.rstrip().split('\t')
    key = e[0]
    if bfkey=='' :
        outF = f'{outPATH}.{key}.txt' 
        O = open(outF, "w")
    elif bfkey != key :
        O.close()
        outF = f'{outPATH}.{key}.txt' 
        O = open(outF, "w")
    
    bfkey = key
    
    O.write(l)
    

    

