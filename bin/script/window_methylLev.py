#! /usr/bin/env python3

import sys
import os
import subprocess as sub

try :
	winsize = int(sys.argv[1])
	IdeoF = sys.argv[2]
	sample_name = sys.argv[3]
	methylCallF = sys.argv[4]
	outF = sys.argv[5]

except IndexError:
    print('     "window_methylLev.py [window size] [Ideogram.bed] [sample_name] [methyl call file] [output]"')
    sys.exit()

    

#Anal_D = os.path.join( mspipe_outputpath,"Analysis")
#IdeoF = os.path.join(Anal_D, "Ideogram.bed")
#vis_paramF = os.path.join(Anal_D, "vis_params.txt")


## READ Ideogram
windowD = {}

for l in open(IdeoF , 'r'):
	(ch,chr_start,chr_end)= l.rstrip().split('\t')[:3]
	
	windowD[ch] = []

	chr_start = int(chr_start)
	chr_end = int(chr_end)
	start = chr_start
	size = winsize

	for i in range(chr_start, chr_end//winsize) :
		windowD[ch].append( (start,start+size) )
		start = start + size

	size = chr_end%winsize
	windowD[ch].append( (start,start+size) )


## READ PARAMS

	
if True:
	O = open(outF, 'w')
	win_i =0
	window_lev_list= []
	ch = ""
	for l in open(methylCallF , "r"):
		col = l.rstrip().split("\t")
	
		bf_chr = ch
		ch = col[0]
		
		pos,all,methyl = map(int, col[1:4])
		pos -=1
		level = (methyl/all)*100

		if ch != bf_chr:
			if bf_chr != "":
				avg_lev = sum(window_lev_list)/len(window_lev_list)
				O.write( f"{bf_chr}\t{win_start}\t{win_end}\t{avg_lev:.3f}\n" )
				window_lev_list = []
			
			window_list = windowD[ch]
			win_i =0
			win_start,win_end = window_list[win_i]
			window_lev_list= []
	
		if pos<win_end and pos >= win_start:
			window_lev_list.append(level)

		else:
			if len(window_lev_list)>0:
				avg_lev = sum(window_lev_list)/len(window_lev_list)
				O.write( f"{ch}\t{win_start}\t{win_end}\t{avg_lev:.3f}\n" )
			
				window_lev_list = []
				win_i +=1
				win_start,win_end = window_list[win_i]
			
			while pos >= win_end:
				line = f"{ch}\t{win_start}\t{win_end}\tNaN\n"
				O.write(line)	
				win_i +=1
				win_start,win_end = window_list[win_i]
		
			window_lev_list.append(level)

	if len(window_lev_list)>0:
		avg_lev = sum(window_lev_list)/len(window_lev_list)
		O.write( f"{ch}\t{win_start}\t{win_end}\t{avg_lev:.3f}\n" )
	
	O.close()
    

