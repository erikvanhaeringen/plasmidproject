#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
MATCHING FUNCTION
-----------------
Checks if one sequence string is present (as a whole) in another sequence string
and returns the result as a list containing sample id, sample sequence and a row 
of binomial values.
"""

def worker(PlasmidDf, ArgDf, SampleID, SampleSeq, SampleNumber):
    DataRow = []
    DataRow.append(SampleID[SampleNumber]) #add sample id to the output vector
    DataRow.append(SampleSeq[SampleNumber]) #add sample sequence to the output vector
    for PlasmidFile in range(len(PlasmidDf)):
        for Plasmid in range(len(PlasmidDf[1][PlasmidFile])):
            if PlasmidDf[1][PlasmidFile][Plasmid] in SampleSeq[SampleNumber]:  
                #print "Found plasmid type!\t", PlasmidDf[0][PlasmidFile][Plasmid], "for", SampleID[SampleNumber]
                DataRow.append("1")
            else:
                DataRow.append("0")
    for ArgFile in range(len(ArgDf)):
        for Arg in range(len(ArgDf[1][ArgFile])):
            if ArgDf[1][ArgFile][Arg] in SampleSeq[SampleNumber]: 
                #print "Found ARG type!\t", ArgDf[0][ArgFile][Arg], "for", SampleID[SampleNumber]  
                DataRow.append("1")
            else:
                DataRow.append("0")   
    return DataRow