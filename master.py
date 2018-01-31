#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import setup.master as setup
import alignment.master as alignment
import statistics.master as statistics
import sys

def main(argv):
    for arg in sys.argv:
        print "argv:",str(arg),"\n"
    
    #STEP 1: Setup
    SampleData, SampleFileNames, SampleCount, PlasmidData, ArgData = setup.main("setup/")
    
    #STEP 2: Allignment
    DataSet = alignment.main(SampleData, SampleFileNames, SampleCount, PlasmidData, ArgData, "alignment/")
    
    #STEP 3: Statistics
    result = statistics.main(DataSet)
    print result

if __name__ == "__main__":
    main(sys.argv)