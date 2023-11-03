#Only LD stats summarizing w.r.t positions:
#sliding window, SLIm ms output
#No divergence here

from __future__ import print_function
import libsequence
import sys
import pandas
import math
import numpy
import argparse
import os

#parsing user given constants
parser = argparse.ArgumentParser(description='Information about number of sliding windows and step size')
parser.add_argument('-fileExt', dest = 'fileExt', action='store', nargs = 1, type = str, help = 'example: _masked.ms')
parser.add_argument('-winSize', dest = 'winSize', action='store', nargs = 1, default = 100, type = int, help = 'size of each sliding window in bp')#500 bp for small, 5000bp for big
parser.add_argument('-stepSize', dest = 'stepSize', action='store', nargs = 1, default = 100, type = int, help = 'size of step size in bp')#250 bp for small, 5000 bp for big
parser.add_argument('-binSize', dest = 'binSize', action='store', nargs = 1, default = 50, type = int, help = 'bin size in bp')#50 bp
parser.add_argument('-regionLen', dest = 'regionLen', action='store', nargs = 1, type = int, help = 'length in bp of region simulated')#Length of coding region simulated
parser.add_argument('-input_folder', dest = 'input_folder', action='store', nargs = 1, type = str, help = 'full path to folder with .ms files')
parser.add_argument('-output_folder', dest = 'output_folder', action='store', nargs = 1, type = str, help = 'full path to folder where you want to write the output')
parser.add_argument('-output_prefix', dest = 'output_prefix', action='store', nargs = 1, type = str, help = 'full path to output file')
args = parser.parse_args()
file_ext = args.fileExt[0]
chr_len =  args.regionLen[0]
bin_size = args.binSize[0]
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

def calculate_mean(l_values):
    if len(l_values) > 0:
        return(str(round(numpy.mean(l_values), 3)))
    else:
        return("NA")
def refresh_dict(BIN_SIZE, TOT_LEN):
    d_LD_BINNED = {}
    d_LD_BINNED['rsq'], d_LD_BINNED['D'], d_LD_BINNED['Dprime'] = {}, {}, {}
    i = 1
    while i <= TOT_LEN:
        BIN = str(i) + "-" + str(i + BIN_SIZE - 1)
        #print(BIN)
        d_LD_BINNED['rsq'][BIN] = []
        d_LD_BINNED['D'][BIN] = []
        d_LD_BINNED['Dprime'][BIN] = []
        i = i + BIN_SIZE
    #print(d_LD_BINNED)
    return(d_LD_BINNED)

def find_my_bin(DIST, BIN_SIZE, TOT_LEN):
    if DIST == 0:
        DIST = 1 #make it so, 0 doesn't mean anything.
    #print(DIST)
    BIN = ""
    i = 1
    while i <= TOT_LEN:
        bin_start = i
        bin_end = bin_start + int(BIN_SIZE) - 1
        if i <= DIST:
            #print("yes 1")
            #print (bin_start)
            #print (bin_end)
            if int(DIST) <= bin_end:
                #print ("yes 2")
                BIN = str(bin_start) + "-" + str(bin_end)
        i = i + BIN_SIZE
    if BIN == "":
        print ("error finding bin for " + str(DIST))
    else:
        #print(BIN)
        return(BIN)

def bin_LD_stats_by_distance(l_LD_stats, BIN_SIZE, d_LD_BINNED):
    for d_LD in l_LD_stats:
        DIST = round((float(d_LD['j'])*chr_len) - (float(d_LD['i'])*chr_len))
        #print(DIST)
        d_LD_BINNED['rsq'][find_my_bin(DIST, int(BIN_SIZE), int(args.winSize[0]))].append(d_LD['rsq'])
        d_LD_BINNED['D'][find_my_bin(DIST, int(BIN_SIZE), int(args.winSize[0]))].append(d_LD['D'])
        d_LD_BINNED['Dprime'][find_my_bin(DIST, int(BIN_SIZE), int(args.winSize[0]))].append(d_LD['Dprime'])
    return()

result =  open(outfolder + "/" + prefix + "_" + str(args.winSize[0]) +  ".LDstats", 'w+')
result.write("filename" + '\t' + "bin" + '\t' + "mean_distance" + '\t' + "rsq" + '\t' + "D" + '\t' + "Dprime" + '\n')

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
        
        #set up a new dict to store LD stats from all windows:
        d_LD_binned = refresh_dict(bin_size, int(args.winSize[0]))
		
        #calculate LD summary statistic in sliding window:
        print ("calculating stats in windows")
        win_name = 1
        for i in range(len(w)):
            wi = w[i]
			#read data to calculate LD based stats:
            #These are pairwise stats. If only 1 site exists, it'll show an error.
            LD_tmp = libsequence.ld(wi)
            #print(LD_tmp)
            bin_LD_stats_by_distance(LD_tmp, bin_size, d_LD_binned)
            win_name = win_name + 1
        
        #write it in a file:
        for s_bin in d_LD_binned['rsq'].keys():
            print(s_bin)
            mean_dist = (int(s_bin.split('-')[1]) + int(s_bin.split('-')[0])) / 2.0
            result.write(f_name_small + '\t' + s_bin + '\t' + str(mean_dist) + '\t' + calculate_mean(d_LD_binned['rsq'][s_bin]) + '\t' + calculate_mean(d_LD_binned['D'][s_bin]) + '\t' + calculate_mean(d_LD_binned['Dprime'][s_bin]) + '\n')
    #except:
    else:
        s_absent = s_absent + 1
        print ("This file does not exist or cannot be read or is empty")
	
    numsim = numsim + 1

result.close()
print ("Number of files not read:" + '\t' + str(s_absent))
print ("Finished")






