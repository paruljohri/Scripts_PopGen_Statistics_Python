#Basic stats, sliding window, SLIm ms output
#adding divergence to this
#python statistics_slidingwindow_pylibseq_general.py -winSize 200 -stepSize 200 -regionLen 2341 -sim_num 21
#
from __future__ import print_function
import libsequence
import sys
import pandas
import math
import argparse


#parsing user given constants
parser = argparse.ArgumentParser(description='Information about number of sliding windows and step size')
parser.add_argument('-winSize', dest = 'winSize', action='store', nargs = 1, type = int, help = 'size of each sliding window in bp')#500 bp for small, 5000bp for big
parser.add_argument('-stepSize', dest = 'stepSize', action='store', nargs = 1, type = int, help = 'size of step size in bp')#250 bp for small, 5000 bp for big
parser.add_argument('-filename', dest = 'filename', action='store', nargs = 1, type = str, help = 'full filename like sim1.ms')
parser.add_argument('-regionLen', dest = 'regionLen', action='store', nargs = 1, type = int, help = 'length in bp of region simulated')#Length of coding region simulated
parser.add_argument('-folder', dest = 'folder', action='store', nargs = 1, type = str, help = 'path to folder')
args = parser.parse_args()
chr_len =  args.regionLen[0]
win_size = args.winSize[0]/float(chr_len)
step_size = args.stepSize[0]/float(chr_len)
folder = args.folder[0]
filename = args.filename[0]


result = open(folder + "/" + filename.split(".")[0] + "_" + str(args.winSize[0]) +	".stats", 'w+')
result.write("simID" + '\t' + "posn" + '\t' + "S" + '\t' + "thetapi" + '\t' + "thetaw" + '\t' + "thetah" + '\t' + "hprime" + '\t' + "tajimasd" +  '\t' + "numSing" + '\t' + "hapdiv" + '\t' + "rsq" + '\t' + "D" + '\t' + "Dprime" + '\t' + "div" + '\n')

#go through the ms file and divergence file if it exists:
f_ms = open(folder + "/" + filename, 'r')
S = get_S(f_ms)
l_Pos = [] #list of positions of SNPs
l_Genos = [] #list of alleles
d_tmp = {}
for line in f_ms:
	line1 = line.strip('\n')
	if "positions" in line1:
		line2 = line1.split()
		i = 0
		for x in line2:
			if "position" not in x:
				l_Pos.append(float(x))
				d_tmp[str(i)] = ""
				i = i + 1
	elif "//" not in line and "segsites" not in line:
		#print (d_tmp)
		i = 0
		while i < len(line1):
			d_tmp[str(i)] = d_tmp[str(i)] + line1[i]
			i = i + 1
#print (d_tmp)
l_data = []
i = 0
while i < len(l_Pos):
	l_Genos.append(d_tmp[str(i)])
	t_tmp = (l_Pos[i], d_tmp[str(i)])
	l_data.append(t_tmp)
	i = i + 1
#print (l_Pos)
#print (l_Genos)


#assign object
sd = libsequence.SimData(l_data)

#define sliding windows:
w = libsequence.Windows(sd,window_size=win_size,step_len=step_size,starting_pos=0.0,ending_pos=1.0)
num_win = len(w)

#calculate summary statistic in sliding window:
numsim=1
print ("calculating stats in windows")
win_name = 1
for i in range(len(w)):
	wi = w[i]
	#print (wi)
	pswi = libsequence.PolySIM(wi)
	result.write("sim" + str(numsim) + '\t' + str(win_name) + '\t' + str(pswi.numpoly()) + '\t' + str(pswi.thetapi()) + '\t' + str(pswi.thetaw()) + '\t' + str(pswi.thetah()) + '\t' + str(pswi.hprime()) + '\t' + str(pswi.tajimasd()) + '\t' + str(pswi.numexternalmutations()) + '\t' + str(pswi.hapdiv()) + '\t')
	#read data to calculate LD based stats:
			
	if len(wi.pos()) >= 5: #These are pairwise stats. If only 1 site exists, it'll show an error.
		#print (i)
		LD_tmp = libsequence.ld(wi)
		LDstats = pandas.DataFrame(LD_tmp)
		meanrsq = sum(LDstats['rsq'])/len(LDstats['rsq'])
		meanD = sum(LDstats['D'])/len(LDstats['D'])
		meanDprime = sum(LDstats['Dprime'])/len(LDstats['Dprime'])
		result.write(str(meanrsq) + '\t' + str(meanD) + '\t' + str(meanDprime) + '\n')
	else:
		result.write("NA" + '\t' + "NA" + '\t' + "NA" + '\n')
	win_name = win_name + 1

print ("Finished")






