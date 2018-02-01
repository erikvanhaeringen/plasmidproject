#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
----------------------------------------
       PLASMID RESISTANCE MATCHER
----------------------------------------
by Erik van Haeringen and Goncalo Macedo
31th Jan 2018

[OPTIONS]
-o <FILE_NAME>      Generates a tab delimited text file of the output to /alignment/results/ 
                    with a name supplied by the user, or as 'Output_DataSet.txt' by default
                    
"""

import setup.master as setup
import alignment.master as alignment
import statistics.master as statistics
import sys

def main(argv):
    print "----------------------------------------\n       PLASMID RESISTANCE MATCHER\n----------------------------------------\nby Erik van Haeringen and Goncalo Macedo\n31th Jan 2018"
    UserOutputFile = None
    for arg in range(len(sys.argv)):
        #checks of option -o is added as an argument
        if sys.argv[arg] == "-o":
            #checks if another argument is given after that and saves this as the outputPath
            if sys.argv[arg+1]:
                UserOutputFile=sys.argv[arg+1] 
            else:
                UserOutputFile="Output_DataSet.txt"
    
    #STEP 1: Setup
    SampleDf,SampleFileNames,SampleCount,PlasmidDf,PlasmidCount,ArgDf,ArgCount = setup.main("setup/")
    
    #STEP 2: Allignment
    ResultDf = alignment.main(SampleDf, SampleFileNames, SampleCount, PlasmidDf, ArgDf, UserOutputFile, "alignment/")
    
    #STEP 3: Statistics
    SubDataSet = statistics.main(ResultDf,PlasmidCount,ArgCount,PlasmidDf,ArgDf,"statistics/")
    
if __name__ == "__main__":
    main(sys.argv)