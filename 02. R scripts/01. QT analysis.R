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
mydir <- '/Users/rstorlund/Library/CloudStorage/OneDrive-UBC/QT/Analysis/SSL-Diving-QT/01. Raw data'
# '/Users/rhea/Library/CloudStorage/OneDrive-UBC/PhD/QT/Analysis/SSL-Diving-QT/01. Raw data'
myfiles <- list.files(path=mydir, pattern="*.csv", full.names=TRUE)
datalist = list()

# # #read in files
# for (i in 1:length(myfiles)) {
  #read in file
  df <- read.csv(file = myfiles[1],
                 sep = ",",
                 stringsAsFactors = FALSE,
                 strip.white = TRUE,
                 na.strings = c("NA",""))

  df <- df %>% 
    mutate(QTav = rowMeans(across(c(X1, X2, X3)), na.rm = TRUE)) %>% #calculate mean of 3 repeat measures
    mutate(same = QT == QTav) #confirm that mean calculations are correct (double check)
# move ahead assuming QT column is correct
  
  #need time since dive start (x) and QT (y)
  #remove all QT = NA, but only after timing variable is correct
  #then plot QT interval over dive time
  #start time = Behaviour == Underwater, Status == Start --> time 0
  #QT before time 0 can be at negative values (or start recording at T=0?)
  
  

#####Plot HR over time#####
#####Summary stats#####
#####Statistics######