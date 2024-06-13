#####QT Analysis#####
#This file plots QT intervals over time for
#Steller sea lions performing stationary dives
#Aligns the HR profiles from the start of the dive
#and plots them altogether
#then calculates average change in HR over the dive

#####Load libraries#####
library(tidyr)
library(dplyr)
library(ggplot2)

#####Read in the data#####

#create a list of files
mydir <- '/Users/rhea/Library/CloudStorage/OneDrive-UBC/PhD/QT/Analysis/SSL-Diving-QT/01. Raw data'
myfiles <- list.files(path=mydir, pattern="*.csv", full.names=TRUE)
datalist = list()

#read in files
for (i in 1:length(myfiles)) {
  #read in file
  df <- read.csv(file = myfiles[1],
                 sep = ",",
                 stringsAsFactors = FALSE,
                 strip.white = TRUE,
                 na.strings = c("NA",""))
  

#####Plot HR over time#####
#####Summary stats#####
#####Statistics######