#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
STEP 3: STATISTICS
------------------
tranform dataset
export to r
do statistics
return output

"""
import re
import pandas as pd
import subprocess

def main(ResultDf,PlasmidCount,ArgCount,PlasmidDf,ArgDf,RelativePath):
    #if called from the main master script add the relative path to all used file paths
    PathAddition=""
    if RelativePath:
        PathAddition = RelativePath
     
    print "\n\n[ Preparing data for statistical analysis ]"
    #build a subdataset that shows the time a specific was found arg per plasmid type, with the last row showing the arg count with no plasmid type
    

    print "parsing plasmid type IDs"
    #generate a list of plasmid types with the subtype information removed for further analysis
    PlasmidIDs = []

    #iterate through the plasmid type sequences and perform a search and replace with regular expressions (re)
    for PlasmidFile in range(len(PlasmidDf)): 
        for Plasmid in range(len(PlasmidDf[0][PlasmidFile])):
            PlasmidIDs.append(re.sub(r"^([A-Za-z0-9/]+).*",r"\1",PlasmidDf[0][PlasmidFile][Plasmid]))
    PlasmidIDs.append("NoType") #add a 'no type' row the plasmid types list to store all resistance matches without a plasmid type match
 
    print "parsing resistance type IDs"
    #generate a list of arg types with the subtype information removed for further analysis  
    ArgIDs = []
    for ArgFile in range(len(ArgDf)):
        for Arg in range(len(ArgDf[0][ArgFile])):
            TempReplace = re.sub(r"([A-Za-z0-9()]+).+",r"\1",ArgDf[0][ArgFile][Arg])
            ArgIDs.append(re.sub(r"(.+)\w\d$",r"\1",TempReplace))
            
    print "create new data set"
    #Create empty dataset with plasmid type as row index and arg type as column headers            
    SubDataSet = pd.DataFrame(columns=(ArgIDs)) #sets column names
    
    #loop that fills the data set with 0's for the number of plasmid types (including no type) times the number of resistance gene types
    for PlasmidNumber in range(PlasmidCount+1): #(+1 for notype row)
        SubDataSet.loc[PlasmidNumber] = [0]*ArgCount
    SubDataSet.index = PlasmidIDs #sets row names based on the shortened plasmid type names 
    
    print "transpose data to new data set"
    #fill the new data set with data from the output of the analysis step
    #using a loop that iterates through the samples 
    for Sample in range(len(ResultDf)):
        #and per sample checks if any plasmid type was found
        #if not it adds the antibiotic types count to the 'no type' row
        if sum([int(i) for i in ResultDf.loc[Sample][2:(PlasmidCount + 2)]]) == 0: #no plasmid type identified in sample          
            for ArgNumber in range(ArgCount):
                if ResultDf.loc[Sample][2+PlasmidCount+ArgNumber]=="1": #uses string type as the output data set is all strings
                    SubDataSet.loc[PlasmidIDs[PlasmidCount],ArgIDs[ArgNumber]]+=1 #add to new dataset
        #if so it adds the antibiotic types count to the row of the plasmid type(s) that were found
        else:
            for Plasmid in len(ResultDf.loc[Sample][2:(PlasmidCount+2)]):
                if ResultDf.loc[Sample][2+Plasmid] == "1":
                    for ArgNumber in range(ArgCount):
                        if ResultDf.loc[Sample][2+PlasmidCount+ArgNumber]=="1": #uses string type as the output data set is all strings
                            SubDataSet.loc[PlasmidIDs[Plasmid],ArgIDs[ArgNumber]]+=1 #add to new dataset


    #export the subdataset for r to a tab delimited file with column and row names included
    outputFile = PathAddition+"results/SubDataSet.txt"
    with open(outputFile, 'w') as file:
        file.writelines('\t'+ Arg for Arg in ArgIDs) #header with the column names
        file.writelines("\n") #line break for the header line
        #iterate through each item in the data set and add it to a temporary string, separated by tabs
        for Plasmid in range(PlasmidCount):
            TempStr=""
            for Arg in range(ArgCount):
                #determin if it is the first item in the row, if so also add the row name
                if TempStr=="": 
                    TempStr=PlasmidIDs[Plasmid] + "\t" + str(SubDataSet.loc[PlasmidIDs[Plasmid],ArgIDs[Arg]]) #add row name (plasmid type) and item to temp string preceded by a tab
                else:
                    TempStr = TempStr + "\t" + str(SubDataSet.loc[PlasmidIDs[Plasmid],ArgIDs[Arg]])
            file.writelines(TempStr + "\n") #write the temporary string to the file plus a line break and continue to the next line (plasmid type)
    print "Exported transposed dataset for r to", outputFile #show user where the file is saved to
    
    print "\n[ STATISTICAL ANALYSIS WITH R SCRIPT ]"
    RScriptPath = PathAddition+"/scripts/summarize.r"
    print "Running", RScriptPath
    subprocess.call (["/usr/bin/Rscript", "--vanilla", RScriptPath])
    
    
    return SubDataSet
    
