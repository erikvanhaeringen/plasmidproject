#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
STEP 2: ALIGNMENT
-----------------
This step will perform the match of the sequences.
As input, it will take the files produced in "Step 1", and it will generate a
tab delimited text file as output.
This file will contain the sample id and sample sequence in the first 2 columns,
than the plasmid types in the number of columns equal to the number of extracted 
plasmid types, and finally the resistance genes in the rest of columns, again as 
many as the number of resistance genes found.
This step makes use of pooled multiprocessing for the function that matches the 
samples to the plasmid and ARG types, and runs as many parallel processes as 
there are cores in the system. 
"""

from __future__ import division
import sys   
import multiprocessing
from multiprocessing import Pool
from functools import partial
import pandas as pd
import matchSequences #the script that matches if a plasmids type and arg type are present in a sample sequence

def main(SampleDf, SampleFileNames, SampleCount, PlasmidDf, ArgDf, UserOutputFile, RelativePath):
    PathAddition=""
    if RelativePath:
        PathAddition = RelativePath
    
    TempDataSet = [] #create new dataSet to hold the results of the sequence matching
    print "\n[ MATCHING SEQUENCES ]"
    print "using %i cores for paralel processing\n" % multiprocessing.cpu_count()

#    print "seq",SampleDf[1][0][0]
#    print "name",SampleDf[0][0][0]
#    print "tname",type(SampleDf[0][0][0])
#    print "file",SampleDf[0][0][0]
#    print "tfile",type(SampleDf[0][0][0])
#    print "sample",SampleDf[0][0][1]
#    print "tsample",type(SampleDf[0][0][1])
#    print "lsample",len(SampleDf[1][0])
#    print dsjkh
    #if __name__ == '__main__':
    for SampleFile in range(len(SampleDf)):
        print "Sample file:", SampleFileNames[SampleFile] 
        p = Pool(multiprocessing.cpu_count())  #create a pool with the amount of processes based on the number of cpus 
        func = partial(matchSequences.worker,PlasmidDf,ArgDf,SampleDf[0][SampleFile],SampleDf[1][SampleFile])
        ProgressCounter = 0.0
        for i in p.imap_unordered(func, range(len(SampleDf[0][SampleFile]))):
            ProgressCounter += 1
            sys.stderr.write('\rprogress: %i%%' % ((ProgressCounter/SampleCount)*100))            
            TempDataSet.append(i)
        p.close()
        p.join()
        
    #To generate a tab delimited text file of the output, add argument -o <FILE_NAME> when calling the main master script
    if UserOutputFile != None:
        OutputFile = PathAddition+"results/"+UserOutputFile   
        with open(OutputFile, 'w') as file:
            file.writelines('\t'.join(i) + '\n' for i in TempDataSet)
        print "\nExported dataset to", OutputFile

    ResultDf = pd.DataFrame(TempDataSet)

    return ResultDf       