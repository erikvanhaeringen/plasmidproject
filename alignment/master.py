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
    #Defines a relative path addition for use with the central coordinating master script
    PathAddition=""
    if RelativePath:
        PathAddition = RelativePath
    
    TempDataSet = [] #create new temporary dataSet to hold the results of the sequence matching
    print "\n[ MATCHING SEQUENCES ]"
    print "using %i cores for paralel processing\n" % multiprocessing.cpu_count()

    for SampleFile in range(len(SampleDf)): #for each sample file do:
        print "Sample file:", SampleFileNames[SampleFile] #print the current sample file to screen
        p = Pool(multiprocessing.cpu_count())  #create a pool with the amount of processes based on the number of cpus 
        func = partial(matchSequences.worker,PlasmidDf,ArgDf,SampleDf[0][SampleFile],SampleDf[1][SampleFile]) #create a partial function call to include more than the iterated variable in the function call
        ProgressCounter = 0.0 #create a progress counter set to 0
        for i in p.imap_unordered(func, range(len(SampleDf[0][SampleFile]))): #for every sample in sample file create a separate task and add to pool
            ProgressCounter += 1 #increase the progress counter
            sys.stderr.write('\rprogress: %i%%' % ((ProgressCounter/SampleCount)*100)) #show progress in percentage processed of current sample file 
            TempDataSet.append(i) #add the output data from the match to the temporary dataset
        p.close()
        p.join() #wait for all tasks in pool to complete
        
    #To generate a tab delimited text file of the output, add argument -o <FILE_NAME> when calling the main master script
    if UserOutputFile != None:
        OutputFile = PathAddition+"results/"+UserOutputFile   
        with open(OutputFile, 'w') as file:
            file.writelines('\t'.join(i) + '\n' for i in TempDataSet) #join all the list elements together with tabs in between and a line break on the end of every line (row)
        print "\nExported dataset to", OutputFile

    ResultDf = pd.DataFrame(TempDataSet) #store temporary data set with matching data to a pandas dataframe

    return ResultDf       