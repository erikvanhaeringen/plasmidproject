# Plasmid Type as factor for spread of ARG

Project developed under the scope of the masters course "Practical Bioinformatics for Biologists", edition of 2018.
The main goal of this work was to develop a working pipeline, in Python, that would simplify the life of researchers analyzing a data set of plasmid sequences.

## Developed by
Erik van Haeringen
Gonçalo Macedo

January 31st 2018

## Getting Started
In the main folder, the user can find this README file, containing information about the general purpose, requirement and options of this program. It also contains a **Master.py** script which controls the three steps of this program: data setup, sequence alignment and statistical analysis. These steps have each their own subfolder, containing a **master.py** script that controls the operations and possible subscripts. For a detailed explanation of the process of each step please see the comments provided in these scripts.

### Prerequisites
This pipeline was developed with Ubuntu v16.04.3. Therefore, it is recommended to run this script in the same version of Ubuntu in which it was developed.
It also requires Python ([picard]https://www.python.org/, min v2.7.14) and R ([picard]https://www.r-project.org/) to be installed. For the required list of modules and packages see [List of modules](#List of Modules)



## Running the tests
* Clone repository from github
'''
git clone git://github.com/erikvanhaeringen/plasmidproject.git
'''
* Go to the cloned folder
'''
cd 'path to repository location'
'''
* Run the script
'''
python Master.py
'''


### Step 1: Setup
This step is the setup for further analysis.
The data consists of plasmid sample sequences (sample), plasmid type sequences (plasmid) and resistance gene sequences (arg).
Plasmid type and resistance gene sequences are first downloaded from the repository and extracted.
Then sample, plasmid, and arg headers and sequences are extracted from the files, reformatted with regular expressions and organized in a dataset.

The function takes a zipfile object, extracts the names of files and removes the first level folder structure to also extract the files in subdirectories to the specified output directory.

### Step 2: Alignment
This step will perform the match of the sequences.
As input, it will take the files produced in "Step 1", and it will generate a tab delimited text file as output
This file will contain the plasmid types in the first column, and the resistance genes in the rest of columns. The values are the number of matches found for each plasmid type.

Note: plasmids with no matches on type are designated as "No type"


### Step 3: Statistics
This step will perform the analysis of the file produced in the Alignment step.
Data will be inserted in R and processed further for development of graphs (ggplot2).


## List of scripts
* /master.py
* /setup/master.py
* /alignment/master.py
* /alignment/matchSequences.py
* /statistics/master.py
* /statistics/scripts/summarize.r


## List of modules
Used in Python
* sys
* urllib
* os
* zipfile
* re
* pandas
* multiprocessing
* functools
* subprocess
* numpy

Used in R
* dplyr
* data.table