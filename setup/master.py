#!usr/bin/env python

# -*- coding: utf-8 -*-
"""
STEP 1: SETUP
-------------
In this script the data is setup for further analysis.
The data consists of plasmid sample sequences (sample), plasmid type sequences 
(plasmid) and resistance gene sequences (arg). Plasmid type and resistance gene 
sequences are first downloaded from the repository and extracted. Then sample, 
plasmid and arg headers and sequences are extracted from the files, reformatted 
with regular expressions and organised in a dataset.
"""

import urllib
import os
from zipfile import ZipFile
import re
import pandas as pd

#This is a function by ekhumoro retrieved from 'https://stackoverflow.com/questions/8689938/extract-files-from-zip-without-keep-the-top-level-folder-with-python-zipfile'
#The function takes a zipfile object, extracts the names of files and removes the first level folder structure to also extract the files in subdirectories to the specified output directory
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

#Function that opens a file containing sequences and performs a regexp search to extract all identifiers (header) and sequences 
def parseSequences(Regexp, Filename, RemoveLineEnding):
    print "\treformatting ", Filename
    File = open(Filename,"r") #open file by provided path
    Data=[[],[]] #make empty list to hold extracted data
    for i in re.findall(Regexp,File.read()): #read file and apply a search with the provided regular expression
        Data[0].append(i[0]) #store identifier string
        Data[1].append(i[1]) #store sequence string
    if RemoveLineEnding == True: #if option for line-ending removal is selected
        for i in range(len(Data[1])): #do for all sequence strings..
            Data[1][i] = re.sub(r"[\r\n\s]+",r"",Data[1][i]) #replace all line breaks and white space characters with nothing using regular expressions
            Data[1][i] = Data[1][i].upper() #convert all letters to upper case
    return Data
    

def main(RelativePath):
    #Defines a relative path addition for use with the central coordinating master script
    PathAddition=""
    if RelativePath:
        PathAddition = RelativePath
     
    #*************** RETRIEVE AND EXTRACT DATA *************** 
    #DOWNLOAD PLASMID SEQUENCE FILE
    print "\n[ Updating plasmid data ]"
    PlasmidFinder = PathAddition+"PlasmidFinder.zip" #destination name and path for the zip file
    print "\tdownloading plasmid data from", PlasmidFinder
    urlPlasmid = 'https://bitbucket.org/genomicepidemiology/plasmidfinder_db/get/1689b7aaee28.zip' #url to download
    urllib.urlretrieve(urlPlasmid, PlasmidFinder) #download from defined url to destination
    
    #EXTRACT CONTENTS OF ZIP FOLDER TO RESULTS FOLDER
    archive = ZipFile(PlasmidFinder) #read the archive to a zifile object
    print "\textracting", PlasmidFinder
    archive.extractall(PathAddition+"results/plasmid", get_members(archive)) #send zipfile object to get_members function and extract the resulting filelist in results folder 
    print "\tremoving", PlasmidFinder
    os.remove(PlasmidFinder) #remove the zip file
    
    #DOWNLOAD ARG SEQUENCE FILE
    print "\n[ Updating arg data ]"
    ArgFinder = PathAddition+"ArgFinder.zip" #destination name and path for the zip file
    print "\tdownloading arg data from", ArgFinder
    urlArg = 'https://bitbucket.org/genomicepidemiology/resfinder_db/get/64794c304abc.zip' #url to download
    urllib.urlretrieve(urlArg, ArgFinder) #download from defined url to destination
    
    #EXTRACT CONTENTS OF ZIP FOLDER TO RESULTS FOLDER      
    archive = ZipFile(ArgFinder) #read the archive to a zifile object
    print "\textracting", ArgFinder
    archive.extractall(PathAddition+"results/arg", get_members(archive)) #send zipfile object to get_members function and extract the resulting filelist in results folder  
    print "\tremoving", ArgFinder
    os.remove(ArgFinder) #remove the zip file
    

    #*************** REFORMAT DATA AND ADD TO DATASET ***************
    print "\n[ Importing data ]"
    
    #IMPORT SAMPLE SEQUENCES
    print "Importing samples"
    Temp=[] #temporary list to store parsed sequence data 
    SampleFileNames=[] #empty list to store all procesed sample files
    SamplePath=PathAddition+"data/samples/" #path where the samples are imported from
    for File in os.listdir(SamplePath): #for each file in destination..
        if File.endswith(".fa"): #providing it is a '.fa' file do:
            SampleFileNames.append(File) #add sample file name to file list
            Regexp = r">(.+?)[\r\n]+([atgcATGC\n\r]+)" #defines the regular expression used to extract the data (should contain 2 captures: ID & Sequence)
            Temp.append(parseSequences(Regexp, SamplePath+File,True)) #parse sequence data with regular expression and store result in the tempory list           
    SampleDf = pd.DataFrame(Temp) #store the tempory list in a pandas datasframe
    
    #IMPORT ARG SEQUENCES
    print "Importing args"
    Temp=[] #temporary list to store parsed sequence data 
    ArgPath=PathAddition+"results/arg/" #path where the plasmid types are imported from
    for File in os.listdir(ArgPath): #for each file in destination..
        if File.endswith(".fsa"): #providing it is a '.fsa' file do:
            Regexp = r">(.+?)[\r\n]+([atgcATGC\n\r]+)" #defines the regular expression used to extract the data (should contain 2 captures: ID & Sequence)
            Temp.append(parseSequences(Regexp, ArgPath+File,True)) #parse sequence data with regular expression and store result in the tempory list              
    ArgDf = pd.DataFrame(Temp) #store the tempory list in a pandas datasframe 
    
    #IMPORT PLASMID SEQUENCES
    print "Importing plasmids.."
    Temp=[] #temporary list to store parsed sequence data 
    PlasmidPath=PathAddition+"results/plasmid/" #path where the arg types are imported from
    for File in os.listdir(PlasmidPath): #for each file in destination..
        if File.endswith(".fsa"): #providing it is a '.fsa' file do:
            Regexp = r">(.+?)[\r\n]+([atgcATGC\n\r]+)" #defines the regular expression used to extract the data (should contain 2 captures: ID & Sequence)
            Temp.append(parseSequences(Regexp, PlasmidPath+File,True)) #parse sequence data with regular expression and store result in the tempory list   
    PlasmidDf = pd.DataFrame(Temp) #store the tempory list in a pandas datasframe
    
    #Count available sequences and print to screen
    SampleCount=sum(len(SampleDf[0][File]) for File in range(len(SampleDf[0])))
    PlasmidCount=sum(len(PlasmidDf[0][File]) for File in range(len(PlasmidDf[0])))        
    ArgCount=sum(len(ArgDf[0][File]) for File in range(len(ArgDf[0])))        
    print "\n| IMPORTED SEQUENCES |\nSamples:%i\nPlasmids:%i\nArg:%i" % (SampleCount, PlasmidCount, ArgCount)
    
    return SampleDf,SampleFileNames,SampleCount,PlasmidDf,PlasmidCount,ArgDf,ArgCount


if __name__ == "__main__":
    main("")