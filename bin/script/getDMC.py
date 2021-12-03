#! /usr/bin/env python3

import sys

try :
    methF= sys.argv[1]
    methdiffF= sys.argv[2]
    q_cutoff = float(sys.argv[3])
    treatments = sys.argv[4]
    outF = sys.argv[5]
except:
    print('[program] [methF] [rawDMBsF] [q_cutoff] [treatment ( ex) 0,0,1,1)] [outF]')
    sys.exit()

try:
    con_name = sys.argv[6]
    case_name = sys.argv[7]
except:
    con_name = 'control'
    case_name =  'case'

qval_diffD = {}
levelD = {}

for l in open(methdiffF , 'r'):
    e = l.rstrip().replace('"','').split(' ')
    if e[0] == 'chr' : continue
    ch, start, end = e[1:4]
    start = int(start)
    end = int(end)
    if not ch in qval_diffD:
        qval_diffD[ch] = {}
        levelD[ch] = {}
    if not start in qval_diffD[ch]:
        qval_diffD[ch][start] = {}
        levelD[ch][start] = {}
    if not end in qval_diffD[ch][start]:
        qval_diffD[ch][start][end] = ( float(0), float(0) )
        levelD[ch][start][end] = []

    qval_diffD[ch][start][end] = ( float(e[6]), float(e[7]) )

O =open(outF , 'w')
treatment = treatments.split(',')
for l in open(methF, 'r'):
    e = l.rstrip().replace('"','').split(' ')
    #print(e)
    if e[0] == 'chr':
        sample_num = (len(e) -4)/3
        if sample_num != float(len(treatment)):
            sys.stderr.write(f" sample not matched in treatment({treatments}) and {methF}\n")
            sys.exit()
        else:
            header = f'chr\tstart\tend\tmethyl(per)_{con_name}\tmethyl(per)_{case_name}\tq_val\tmeth_diff\n'
            O.write(header)
    else :
        ch, start, end, strand= e[1:5]
        start, end = int(start), int(end)
        
        qval , diff = qval_diffD[ch][start][end]
        if qval > q_cutoff: continue

        total_control, methyl_control, total_case,  methyl_case = 0,0,0,0
        CpG = list(map(lambda x:int(x) ,e[5:]))
        for i in range(0,int(sample_num)):
            if treatment[i] == '0':
                total_control += CpG[i*3]
                methyl_control += CpG[i*3+1]
            if treatment[i] == '1':
                total_case += CpG[i*3]
                methyl_case += CpG[i*3+1]
        
        p_methyl_control =(methyl_control/total_control)*100
        p_methyl_case =(methyl_case/total_case)*100
        line = f'{ch}\t{start}\t{end}\t{total_control}\t{methyl_control}\t{p_methyl_control}\t{total_case}\t{methyl_case}\t{p_methyl_case}\n'
        sys.stderr.write(line)
        
        levelD[ch][start][end] = [strand, str(p_methyl_control), str(p_methyl_case), qval, diff]

path = outF.split('/')
reformF = '/'.join(path[:-1]) + '/reform.' + path[-1]
O2 = open(reformF, 'w' )
header = f'chr\tstart\tend\tmethyl(per)_{case_name}\tmethyl(per)_{con_name}\n'
O2.write(header)
for ch in sorted(levelD.keys()):
    for start in sorted(levelD[ch].keys()):
        for end in sorted(levelD[ch][start].keys()):
            line = f'{ch}\t{int(start)-1}\t{end}\t'
            
            if levelD[ch][start][end] == []:
                continue
            else:
                strand, con, case, qval, diff = levelD[ch][start][end]
                line1 =line +  f'{con}\t{case}\t{qval}\t{diff}\n'
                O.write(line1)
                line2 = line + f'{case}\t{con}\t{qval}\t{diff}\n'
                O2.write(line2)



O.close()
O2.close()
            

