#! /usr/bin/env python

import sys


if __name__ == '__main__':
    
    bs2_result = sys.argv[1]
    out = sys.argv[2]
    
    lines = {}
    lines['CG'] = []
    lines['CHG'] = []
    lines['CHH'] = []
    for l in open(bs2_result, 'r'):
        ch, nt, pos, context, di_context, level, c_count, ct_count = l.rstrip().split("\t")
        lines[context].append( f'{ch}\t{pos}\t{pos}\t{float(level)*100}\t{c_count}\t{int(ct_count)-int(c_count)}\n' )

    O_CpG = open( out + '_CpG.cov.txt', 'w')
    O_CpG.write( "".join(lines['CG']) )
    O_CpG.close()
    
    O_CHG = open( out + '_CHG.cov.txt', 'w')
    O_CHG.write( "".join(lines['CHG']) )
    O_CHG.close()

    
    O_CHH = open( out + '_CHH.cov.txt', 'w')
    O_CHH.write( "".join(lines['CHH']) )
    O_CHH.close()


