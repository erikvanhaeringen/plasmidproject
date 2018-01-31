#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from __future__ import division
import sys

#import numpy as np
#import re
#import os
#            
#def parseSequences(Regexp, Filename, RemoveLineEnding):
#    Data = np.fromregex(Filename, Regexp,[('id', object), ('sequence', object)])
#    if RemoveLineEnding == True:
#        for i in range(len(Data)):
#            #Data[i][1].replace('\n', '').replace('\r', '') 
#            Data[i][1] = (re.sub(r"[\r\n\s]+",r"",Data[i][1]))
#            Data[i][1] = Data[i][1].upper()
#    return Data        
#
##*************** IMPORT SAMPLE SEQUENCES ***************
#Temp=[]
#SampleFileNames=[]
#SamplePath="../setup/data/test/"
#for File in os.listdir(SamplePath):
#    if File.endswith(".fa"):
#        SampleFileNames.append(File)
#        Regexp = r">(.+?)[\r\n]+([atgcATGC\n\r]+)"
#        Temp.append(parseSequences(Regexp, SamplePath+File,True))   
#SampleData = np.array(Temp)
#
##*************** IMPORT ARG SEQUENCES ***************
#Temp=[]
#ArgPath="../setup/results/arg/"
#for File in os.listdir(ArgPath):
#    if File.endswith(".fsa"):
#        Regexp = r">(.+?)[\r\n]+([atgcATGC\n\r]+)"
#        Temp.append(parseSequences(Regexp, ArgPath+File,True))   
#ArgData = np.array(Temp)
#
##*************** IMPORT PLASMID SEQUENCES ***************
#Temp=[]
#PlasmidPath="../setup/results/plasmid/"
#for File in os.listdir(PlasmidPath):
#    if File.endswith(".fsa"):
#        Regexp = r">(.+?)[\r\n]+([atgcATGC\n\r]+)"
#        Temp.append(parseSequences(Regexp, PlasmidPath+File,True))   
#PlasmidData = np.array(Temp)
#
##*************** DATASET CONSTRUCTION ***************
##Count available sequences
#PlasmidCount = 0
#for File in PlasmidData:
#    PlasmidCount = PlasmidCount + len(File)
#ArgCount = 0
#for File in ArgData:
#    ArgCount = ArgCount + len(File)
#SampleCount = 0
#for File in SampleData:
#    SampleCount = SampleCount + len(File)
#print "\n| IMPORTED SEQUENCES |\nSamples:%i\nPlasmids:%i\nArg:%i" % (SampleCount, PlasmidCount, ArgCount)
#
#
#n = 2+PlasmidCount+ArgCount
#OutputDB = np.empty(shape=[SampleCount, n],dtype=object)
#
#for SampleFile in SampleData:
#    for Sample  in range(len(SampleFile)):
#        OutputDB[Sample][0]=SampleFile[Sample][0]
#        OutputDB[Sample][1]=SampleFile[Sample][1]


#*************** SEQUENCE MATCHING ***************

#print "\n| PLASMID MATCHES |"
#for SampleFile in SampleData:
#    sCounter = 0
#    for Sample in SampleFile:
#        for PlasmidFile in PlasmidData:
#            pCounter = 0
#            for Plasmid in PlasmidFile:
#                if Plasmid[1] in Sample[1]:
#                    print "Sample:\t%s\tPlasmid:\t%s" % (Sample[0], Plasmid[0])
#                    OutputDB[sCounter][2+pCounter] = 1
#                pCounter += 1
#        sCounter += 1
#                    
#print "\n| ARG MATCHES |"
#for SampleFile in SampleData:
#    for Sample in SampleFile:
#        for ArgFile in ArgData:
#            for Arg in ArgFile:
#                if Arg[1] in Sample[1]:  
#                    print "Sample:\t%s\tArg:\t%s" % (Sample[0], Arg[0])
                 
                    
import multiprocessing
from multiprocessing import Pool
from functools import partial
import matchSequences #the script that contains the function for determining which plasmids and arg are present in a sample sequence


def main(SampleData, SampleFileNames, SampleCount, PlasmidData, ArgData, RelativePath):
    PathAddition=""
    if RelativePath:
        PathAddition = RelativePath
        
    DataSet = [] #create new dataSet to hold the results of the sequence matching
    print "\n[ MATCHING SEQUENCES ]"
    print "using %i cores for paralel processing\n" % multiprocessing.cpu_count()
    
    #if __name__ == '__main__':
    FileProgress = 0
    for SampleFile in SampleData:
        print "Sample file:", SampleFileNames[FileProgress] 
        p = Pool(multiprocessing.cpu_count())  #create a pool with the amount of processes based on the number of cpus 
        func = partial(matchSequences.worker,PlasmidData,ArgData)
        ProgressCounter = 0.0
        for i in p.imap_unordered(func, SampleFile):
            ProgressCounter += 1
            sys.stderr.write('\rprogress: %i%%' % ((ProgressCounter/SampleCount)*100))            
            DataSet.append(i)
        p.close()
        p.join()
        FileProgress += 1
    
    
    outputFile = PathAddition+"results/DataSet.txt"
    with open(outputFile, 'w') as file:
        file.writelines('\t'.join(i) + '\n' for i in DataSet)
    print "\nExported dataset to", outputFile
    
    Header = []
    Header.append("SampleID")
    Header.append("SampleSeq")   
    
    for PlasmidFile in PlasmidData:
        for Plasmid in PlasmidFile:
            Header.append(Plasmid[0].strip())
 
    for ArgFile in ArgData:
        for Arg in ArgFile:
            Header.append(Arg[0].strip())
    
    outputFile = PathAddition+"results/Header.txt"
    with open(outputFile, 'w') as file:
        file.writelines('\t'+ i for i in Header)
    print "Exported header to", outputFile 

    return DataSet       