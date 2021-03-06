#! /usr/bin/Rscript

#---
#title: "PlasmidProject"
#author: "Goncalo"
#date: "January 30, 2018"
#output: html_document
#editor_options: 
#  chunk_output_type: console
#---

#setwd("~/Documents/plasmidproject/statistics/scripts/")

# Load needed libraries needed (package)
#install.packages("data.table")
#install.packages("dplyr")
library(dplyr)
library(data.table)

# Load dataset (make sure to change all spaces, empty cells and strange characters)
alldata<-read.table("../results/SubDataSet.txt", sep="\t", header=TRUE,  na.strings=c("", "NA"))


# Confirm consistency of data loaded/generated
#str(alldata)
#ncol(alldata)
#nrow(alldata)
#summary(alldata)


# Data manipulations

alldata_treated <- alldata
#colnames(alldata_treated)

#check caracteristics of new dataset
#str(alldata_treated)
#ncol(alldata_treated)
#nrow(alldata_treated)


# make new variables (family of antibiotic) based on sum of different genes of that family

alldata_treated$TET <- apply(alldata_treated[,c(2:105)], 1, sum) #tetracycline

alldata_treated$BLA <- apply(alldata_treated[,c(106:1489)], 1, sum) #Beta-lactams

alldata_treated$MAC <- apply(alldata_treated[,c(1490:1621)], 1, sum) #Macrolide

alldata_treated$NIM <- apply(alldata_treated[,c(1622:1635)], 1, sum) #Nitroimidazole

alldata_treated$COL <- apply(alldata_treated[,c(1636:1655)], 1, sum) #colistin

alldata_treated$FOS <- apply(alldata_treated[,c(1656:1696)], 1, sum) #fosfomycin

alldata_treated$AMN <- apply(alldata_treated[,c(1698:1865)], 1, sum) #Aminoglycoside

alldata_treated$QNL <- apply(alldata_treated[,c(1866:1968)], 1, sum) #Quinolone

alldata_treated$PHN <- apply(alldata_treated[,c(1969:2002)], 1, sum) #Phenicol

alldata_treated$RIF <- apply(alldata_treated[,c(2007:2015)], 1, sum) #rifampicin

alldata_treated$VAN <- apply(alldata_treated[,c(2016:2038)], 1, sum) #vancomycin

alldata_treated$SUL <- apply(alldata_treated[,c(2162:2210)], 1, sum) #sulphonamide

alldata_treated$FUS <- apply(alldata_treated[,c(2211:2212)], 1, sum) #fusidic acid

alldata_treated$TRM <- apply(alldata_treated[,c(2213:2265)], 1, sum) #trimethprim


# New table with only sums of ARG, by antibiotic class
alldata_final <- alldata_treated[,c(1,2266:2279)]
#str(alldata_final)

# remove rows that are all zeros
alldata_final <- alldata_final[!!rowSums(abs(alldata_final[-c(1:2)])),]
#str(alldata_final)

print(alldata_final) #show output
write.table(alldata_final, file = '../results/summaryTable.csv')
return(alldata_final)
