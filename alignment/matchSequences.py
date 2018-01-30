#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 22:34:35 2018

@author: erik
"""
#import os
#from multiprocessing import Queue 

def worker(PlasmidData, ArgData, Sample):
    #print 'Worker %s started on %s' % (os.getpid(),Sample[0])
    DataRow = []
    DataRow.append(Sample[0]) #add sample id to the output vector
    DataRow.append(Sample[1]) #add sample sequence to the output vector
    for PlasmidFile in PlasmidData:
        for Plasmid in PlasmidFile:
            if Plasmid[1] in Sample[1]:
                #print "Sample:\t%s\tPlasmid:\t%s" % (Sample[0], Plasmid[0])
                #q.put("Sample:\t%s\tPlasmid:\t%s" % (Sample[0], Plasmid[0]))                
                DataRow.append("1")
            else:
                DataRow.append("0")
    for ArgFile in ArgData:
        for Arg in ArgFile:
            if Arg[1] in Sample[1]:
                #print "Sample:\t%s\tArg:\t%s" % (Sample[0], Arg[0])
                #q.put("Sample:\t%s\tArg:\t%s" % (Sample[0], Arg[0]))                
                DataRow.append("1")
            else:
                DataRow.append("0")   
    return DataRow
    #return "done"