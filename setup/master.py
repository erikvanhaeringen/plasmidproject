#!usr/bin/env python

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import urllib
import os
from zipfile import ZipFile
import numpy as np
import re

#used this function by ekhumoro found on 'https://stackoverflow.com/questions/8689938/extract-files-from-zip-without-keep-the-top-level-folder-with-python-zipfile'
def get_members(zip):
    parts = []
    # get all the path prefixes
    for name in zip.namelist():
        # only check files (not directories)
        if not name.endswith('/'):
            # keep list of path elements (minus filename)
            parts.append(name.split('/')[:-1])
    # now find the common path prefix (if any)
    prefix = os.path.commonprefix(parts)
    if prefix:
        # re-join the path elements
        prefix = '/'.join(prefix) + '/'
    # get the length of the common prefix
    offset = len(prefix)
    # now re-set the filenames
    for zipinfo in zip.infolist():
        name = zipinfo.filename
        # only check files (not directories)
        if len(name) > offset:
            # remove the common prefix
            zipinfo.filename = name[offset:]
            yield zipinfo

def parseSequences(Regexp, Filename, RemoveLineEnding):
    print "\treformatting ", Filename
    Data = np.fromregex(Filename, Regexp,[('id', object), ('sequence', object)])
    if RemoveLineEnding == True:
        for i in range(len(Data)):
            #Data[i][1].replace('\n', '').replace('\r', '') 
            Data[i][1] = (re.sub(r"[\r\n\s]+",r"",Data[i][1]))
            Data[i][1] = Data[i][1].upper()
    return Data   
    

def main(RelativePath):
    PathAddition=""
    if RelativePath:
        PathAddition = RelativePath
    #DOWNLOAD PLASMID SEQUENCE FILE
    print "\n[ Updating plasmid data ]"
    PlasmidFinder = PathAddition+"PlasmidFinder.zip"
    print "\tdownloading plasmid data from", PlasmidFinder
    urlPlasmid = 'https://bitbucket.org/genomicepidemiology/plasmidfinder_db/get/1689b7aaee28.zip'
    urllib.urlretrieve(urlPlasmid, PlasmidFinder)
    
    #extract contens of the first folder to results folder
    archive = ZipFile(PlasmidFinder)
    print "\textracting", PlasmidFinder
    archive.extractall(PathAddition+"results/plasmid", get_members(archive)) 
    print "\tremoving", PlasmidFinder
    os.remove(PlasmidFinder)
    
    #DOWNLOAD ARG SEQUENCE FILE
    print "\n[ Updating arg data ]"
    ArgFinder = PathAddition+"ArgFinder.zip"
    print "\tdownloading arg data from", ArgFinder
    urlArg = 'https://bitbucket.org/genomicepidemiology/resfinder_db/get/64794c304abc.zip'
    urllib.urlretrieve(urlArg, ArgFinder)
    
    #extract contens of the first folder to results folder            
    archive = ZipFile(ArgFinder)
    print "\textracting", ArgFinder
    archive.extractall(PathAddition+"results/arg", get_members(archive)) 
    print "\tremoving", ArgFinder
    os.remove(ArgFinder)
    
    
    
    
    #REFORMAT DATA AND ADD TO DATASET
    print "\n[ Importing data ]"
    #*************** IMPORT SAMPLE SEQUENCES ***************
    print "\timporting samples"
    Temp=[]
    SampleFileNames=[]
    SamplePath=PathAddition+"data/test/"
    for File in os.listdir(SamplePath):
        if File.endswith(".fa"):
            SampleFileNames.append(File)
            Regexp = r">(.+?)[\r\n]+([atgcATGC\n\r]+)"
            Temp.append(parseSequences(Regexp, SamplePath+File,True))   
    SampleData = np.array(Temp)
    
    #*************** IMPORT ARG SEQUENCES ***************
    print "\timporting args"
    Temp=[]
    ArgPath=PathAddition+"results/arg/"
    for File in os.listdir(ArgPath):
        if File.endswith(".fsa"):
            Regexp = r">(.+?)[\r\n]+([atgcATGC\n\r]+)"
            Temp.append(parseSequences(Regexp, ArgPath+File,True))   
    ArgData = np.array(Temp)
    
    #*************** IMPORT PLASMID SEQUENCES ***************
    print "\timporting plasmids.."
    Temp=[]
    PlasmidPath=PathAddition+"results/plasmid/"
    for File in os.listdir(PlasmidPath):
        if File.endswith(".fsa"):
            Regexp = r">(.+?)[\r\n]+([atgcATGC\n\r]+)"
            Temp.append(parseSequences(Regexp, PlasmidPath+File,True))   
    PlasmidData = np.array(Temp)
    
    #*************** DATASET CONSTRUCTION ***************
    #Count available sequences
    PlasmidCount = 0
    for File in PlasmidData:
        PlasmidCount = PlasmidCount + len(File)
    ArgCount = 0
    for File in ArgData:
        ArgCount = ArgCount + len(File)
    SampleCount = 0
    for File in SampleData:
        SampleCount = SampleCount + len(File)
    print "\n| IMPORTED SEQUENCES |\nSamples:%i\nPlasmids:%i\nArg:%i" % (SampleCount, PlasmidCount, ArgCount)
    
    return SampleData, SampleFileNames, SampleCount, PlasmidData, ArgData
