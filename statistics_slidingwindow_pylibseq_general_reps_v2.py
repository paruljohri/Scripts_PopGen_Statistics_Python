#Slightly modified but better:
#Basic stats, sliding window, uses SLiM ms output
#No divergence here, only polymorphism-based stats from a single population
#This program finds every *.ms file in the input_folder and calculates stats for all the files. 
#A single output file is created in output_folder with the name output_prefix + "_" + window_size + ".stats".
#How to run:
#python statistics_slidingwindow_pylibseq_general_reps_v2.py -fileExt ".ms" -winSize 500 -stepSize 500 -regionLen 10000 -input_folder /path/to/folder/with/ms/files -output_folder /path/to/folder/with/results -output_prefix name_of_result_files

from __future__ import print_function
import libsequence
import sys
import pandas
import math
import argparse
import os

#parsing user given constants
parser = argparse.ArgumentParser(description='Information about number of sliding windows and step size')
parser.add_argument('-fileExt', dest = 'fileExt', action='store', nargs = 1, type = str, help = 'example: _masked.ms')
parser.add_argument('-winSize', dest = 'winSize', action='store', nargs = 1, default = 100, type = int, help = 'size of each sliding window in bp')#500 bp for small, 5000bp for big
parser.add_argument('-stepSize', dest = 'stepSize', action='store', nargs = 1, default = 100, type = int, help = 'size of step size in bp')#250 bp for small, 5000 bp for big
parser.add_argument('-regionLen', dest = 'regionLen', action='store', nargs = 1, type = int, help = 'length in bp of region simulated')#Length of coding region simulated
parser.add_argument('-input_folder', dest = 'input_folder', action='store', nargs = 1, type = str, help = 'full path to folder with .ms files')
parser.add_argument('-output_folder', dest = 'output_folder', action='store', nargs = 1, type = str, help = 'full path to folder where you want to write the output')
parser.add_argument('-output_prefix', dest = 'output_prefix', action='store', nargs = 1, type = str, help = 'full path to output file')
args = parser.parse_args()
file_ext = args.fileExt[0]
chr_len =  args.regionLen[0]
win_size = args.winSize[0]/float(chr_len)
step_size = args.stepSize[0]/float(chr_len)
infolder = args.input_folder[0]
outfolder = args.output_folder[0]
prefix = args.output_prefix[0]

def read_ms_file(f_MS):
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
    return(l_data)

result =  open(outfolder + "/" + prefix + "_" + str(args.winSize[0]) +  ".stats", 'w+')
result.write("filename" + '\t' + "posn" + '\t' + "S" + '\t' + "thetapi" + '\t' + "thetaw" + '\t' + "thetah" + '\t' + "hprime" + '\t' + "tajimasd" +  '\t' + "numSing" + '\t' + "hapdiv" + '\t' + "rsq" + '\t' + "D" + '\t' + "Dprime" + '\t' + "BriansD" + '\n')

#go through all simulation replicates and read data into pylibseq format
#addin the option of ignoring some files if they don't exist
os.system("ls " + infolder + "/*" + file_ext + " > " + outfolder + "/" + prefix + ".list")


f_list = open(outfolder + "/" + prefix + ".list", 'r')
numsim = 1
s_absent = 0
for Aline in f_list:
    Aline1 = Aline.strip('\n')
    f_name = Aline1#.split(".")[0]
    f_name_small = Aline1.split("/").pop()
    print ("Reading file:" + Aline1)
    #try:
    if numsim > 0:
        f_ms = open(f_name, 'r')
        l_data = read_ms_file(f_ms)
        f_ms.close()

		#assign object
        sd = libsequence.SimData(l_data)

		#define sliding windows:
        w = libsequence.Windows(sd,window_size=win_size,step_len=step_size,starting_pos=0.0,ending_pos=1.0)
        num_win = len(w)

		#calculate summary statistic in sliding window:
        print ("calculating stats in windows")
        win_name = 1
        for i in range(len(w)):
            wi = w[i]
            #print (wi)
            pswi = libsequence.PolySIM(wi)
            result.write(f_name_small + '\t' + str(win_name) + '\t' + str(pswi.numpoly()) + '\t' + str(pswi.thetapi()) + '\t' + str(pswi.thetaw()) + '\t' + str(pswi.thetah()) + '\t' + str(pswi.hprime()) + '\t' + str(pswi.tajimasd()) + '\t' + str(pswi.numexternalmutations()) + '\t' + str(pswi.hapdiv()) + '\t')
			#read data to calculate LD based stats:
            #These are pairwise stats. If only 1 site exists, it'll show an error.
            LD_tmp = libsequence.ld(wi)
            LDstats = pandas.DataFrame(LD_tmp)
            #print(len(LDstats))
            if len(LDstats) > 0:
                meanrsq = sum(LDstats['rsq'])/len(LDstats['rsq'])
                meanD = sum(LDstats['D'])/len(LDstats['D'])
                meanDprime = sum(LDstats['Dprime'])/len(LDstats['Dprime'])
                result.write(str(meanrsq) + '\t' + str(meanD) + '\t' + str(meanDprime) + '\t')
            else:
                result.write("NA" + '\t' + "NA" + '\t' + "NA" + '\t') 
            #Brian's D:                                                         
            if pswi.thetaw() > 0:                                               
                result.write(str(1.0 - (pswi.thetapi()/pswi.thetaw())) + '\n')  
            else:                                                               
                result.write("NA" + '\n')

            win_name = win_name + 1
    #except:
    else:
        s_absent = s_absent + 1
        print ("This file does not exist or cannot be read or is empty")
	
    numsim = numsim + 1

print ("Number of files not read:" + '\t' + str(s_absent))
print ("Finished")






