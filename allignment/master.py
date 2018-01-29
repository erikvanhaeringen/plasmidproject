#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import re
import os

def matchSequences(Seq1, Seq2):
    if Seq1 in Seq2:
            return True
            
def parseSequences(Regexp, Filename, RemoveLineEnding):
    Data = np.fromregex(Filename, Regexp,[('id', object), ('sequence', object)])
    if RemoveLineEnding == True:
        for i in range(len(Data)):
            #Data[i][1].replace('\n', '').replace('\r', '') 
            Data[i][1] = (re.sub(r"[\r\n\s]+",r"",Data[i][1]))
            Data[i][1] = Data[i][1].upper()
    return Data        




#IMPORT SAMPLE SEQUENCES
temp=[]
SamplePath="../setup/data/test/"
for file in os.listdir(SamplePath):
    if file.endswith(".fa"):
        Regexp = r">(.+?)[\r\n]+([atgcATGC\n\r]+)"
        temp.append(parseSequences(Regexp, SamplePath+file,True))   
SampleData = np.array(temp)

#IMPORT ARG SEQUENCES
temp=[]
ArgPath="../setup/results/arg/"
for file in os.listdir(ArgPath):
    if file.endswith(".fsa"):
        Regexp = r">(.+?)[\r\n]+([atgcATGC\n\r]+)"
        temp.append(parseSequences(Regexp, ArgPath+file,True))   
ArgData = np.array(temp)

#IMPORT PLASMID SEQUENCES
temp=[]
PlasmidPath="../setup/results/plasmid/"
for file in os.listdir(PlasmidPath):
    if file.endswith(".fsa"):
        Regexp = r">(.+?)[\r\n]+([atgcATGC\n\r]+)"
        temp.append(parseSequences(Regexp, PlasmidPath+file,True))   
PlasmidData = np.array(temp)





print "PLASMID MATCHES"
for SampleFile in SampleData:
    for Sample in SampleFile:
        for PlasmidFile in PlasmidData:
            for Plasmid in PlasmidFile:
                if Plasmid[1] in Sample[1]:
                    print "Sample:\t%s\tPlasmid:\t%s" % (Sample[0], Plasmid[0])
print "ARG MATCHES"
for SampleFile in SampleData:
    for Sample in SampleFile:
        for ArgFile in ArgData:
            for Arg in ArgFile:
                if Arg[1] in Sample[1]:  
                    print "Sample:\t%s\tArg:\t%s" % (Sample[0], Arg[0])
            