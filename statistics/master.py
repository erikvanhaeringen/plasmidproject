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

def main(ResultDf,PlasmidCount,ArgCount,PlasmidDf,ArgDf,RelativePath):

    PathAddition=""
    if RelativePath:
        PathAddition = RelativePath
     
    print "\n\n[ Preparing data for statistical analysis ]"
    #build a subdataset that shows the time a specific was found arg per plasmid type, with the last row showing the arg count with no plasmid type
    
    #generate a list of plasmid types with the subtype information removed for further analysis
    print "parsing plasmid type IDs"
    PlasmidIDs = []

    for PlasmidFile in range(len(PlasmidDf)):
        for Plasmid in range(len(PlasmidDf[0][PlasmidFile])):
            PlasmidIDs.append(re.sub(r"^([A-Za-z0-9/]+).*",r"\1",PlasmidDf[0][PlasmidFile][Plasmid]))
    PlasmidIDs.append("NoType")

    #generate a list of arg types with the subtype information removed for further analysis   
    print "parsing resistance type IDs"
    ArgIDs = []
    for ArgFile in range(len(ArgDf)):
        for Arg in range(len(ArgDf[0][ArgFile])):
            TempReplace = re.sub(r"([A-Za-z0-9()]+).+",r"\1",ArgDf[0][ArgFile][Arg])
            ArgIDs.append(re.sub(r"(.+)\w\d$",r"\1",TempReplace))
            
    print "create new data set"
    #Create empty dataset with plasmid type as row index and arg type as column headers            
    SubDataSet = pd.DataFrame(columns=(ArgIDs)) #sets column names
    for PlasmidNumber in range(PlasmidCount+1): #(+1 for notype row)
        SubDataSet.loc[PlasmidNumber] = [0]*ArgCount
    SubDataSet.index = PlasmidIDs #sets row names
    
    print "transpose data to new data set"
    for Sample in range(len(ResultDf)):
        if sum([int(i) for i in ResultDf.loc[Sample][2:(PlasmidCount + 2)]]) == 0: #no plasmid type identified in sample          
            for ArgNumber in range(ArgCount):
                if ResultDf.loc[Sample][2+PlasmidCount+ArgNumber]=="1":
                    SubDataSet.loc[PlasmidIDs[PlasmidCount],ArgIDs[ArgNumber]]+=1
        else:
            for Plasmid in len(ResultDf.loc[Sample][2:(PlasmidCount+2)]):
                if ResultDf.loc[Sample][2+Plasmid] == "1":
                    for ArgNumber in range(ArgCount):
                        if ResultDf.loc[Sample][2+PlasmidCount+ArgNumber]=="1":
                            SubDataSet.loc[PlasmidIDs[Plasmid],ArgIDs[ArgNumber]]+=1


    #export the subdataset for r to a tab delimited file with column and row names included
    outputFile = PathAddition+"results/SubDataSet.txt"
    with open(outputFile, 'w') as file:
        file.writelines('\t'+ Arg for Arg in ArgIDs)
        file.writelines("\n")
        PlasmidCounter = 0
        for index, row in SubDataSet.iterrows():
            TempStr=""
            for Arg in row.iteritems():
                if TempStr=="":
                    TempStr=PlasmidIDs[PlasmidCounter] + "\t" + str(Arg[1])
                else:
                    TempStr = TempStr + "\t" + str(Arg[1])
            file.writelines(TempStr + "\n")
            PlasmidCounter+=1
    print "Exported plasmid versus arg dataset to", outputFile   
    
    return SubDataSet

#import re
#
#def main(DataSet,PlasmidCount,ArgCount,PlasmidData,ArgData,RelativePath):
#    PathAddition=""
#    if RelativePath:
#        PathAddition = RelativePath
#        
#    #build a subdataset that shows the time a specific was found arg per plasmid type, with the last row showing the arg count with no plasmid type
#    SubDataSet = [[0 for col in range(2264)] for row in range(271)]
#    for Row in DataSet:
#        if sum([int(i) for i in Row[2:(PlasmidCount + 2)]]) == 0: #no plasmid type identified in sample
#            for ArgNumber in range(ArgCount):
#                if Row[(2+PlasmidCount+ArgNumber)]=="1":
#                    SubDataSet[270][ArgNumber]+=1
#        else:
#            PlasmidNumber = 0
#            for Plasmid in Row[2:(PlasmidCount + 2)]:
#                if Plasmid == "1":
#                    for ArgNumber in range(ArgCount):
#                        if Row[(2+PlasmidCount+ArgNumber)]=="1":
#                            SubDataSet[PlasmidNumber][ArgNumber]+=1
#                PlasmidNumber+=1
#            
#    #generate a list of plasmid types with subtype information removed for analysis
#    PlasmidIDs = []
#    for PlasmidFile in PlasmidData:
#        for PlasmidType in PlasmidFile:
#            PlasmidIDs.append(re.sub(r"^([A-Za-z0-9/]+).*",r"\1",PlasmidType[0]))
#    PlasmidIDs.append("NoType")
#    
#    #generate a list of arg types with subtype information removed for analysis    
#    ArgIDs = []
#    for ArgFile in ArgData:
#        for ArgType in ArgFile: 
#            TempReplace = re.sub(r"([A-Za-z0-9()]+).+",r"\1",ArgType[0])
#            ArgIDs.append(re.sub(r"(.+)\w\d$",r"\1",TempReplace))
#            
#    #export the subdataset for r to a tab delimited file with column and row names included
#    outputFile = PathAddition+"results/SubDataSet.txt"
#    with open(outputFile, 'w') as file:
#        file.writelines('\t'+ Arg for Arg in ArgIDs)
#        file.writelines("\n")
#        PlasmidCounter = 0
#        for PlasmidType in SubDataSet:
#            TempStr=""
#            for Arg in PlasmidType:
#                if TempStr=="":
#                    TempStr=PlasmidIDs[PlasmidCounter] + "\t" + str(Arg)
#                else:
#                    TempStr = TempStr + "\t" + str(Arg)
#            file.writelines(TempStr + "\n")
#            PlasmidCounter+=1
#    print "Exported plasmid versus arg dataset to", outputFile
#    
#    
#    
#    
#    
#    
#    
#    #r part
##    from numpy import *
##    import scipy as sp
##    from pandas import *
##    from rpy2.robjects.packages import importr
##    import rpy2.robjects as ro
##    import pandas.rpy.common as com
##
##    rdf = com.convert_to_r_dataframe(df)    
##        
#    import pandas as pd
#    PandaDataSet = pd.DataFrame(columns=(ArgIDs)) #sets column names
#    for PlasmidNumber in range(len(SubDataSet)):
#        PandaDataSet.loc[PlasmidNumber] = SubDataSet[PlasmidNumber]
#    PandaDataSet.index = PlasmidIDs #sets row names
#
#    
#   
#    
#    
#    return SubDataSet