#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
STEP 2: ALIGNMENT
-----------------

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