# ----------------------------------------------------------------- #
# name              Percentage_responses.R                          #
# created by        Shauni Van Herck                                #
# description       used to calculate percentage significant ASSRs  #
# version.string    R version 3.5.1 (2018-07-02)                    #
# platform          x86_64-w64-mingw32                              #
# date created      17/01/2020                                      #
# ----------------------------------------------------------------- #

# load libraries -------------------
# reading
library(readxl)
# reshaping
library(stringr)      # split data into numbers/letters
library(reshape)
library(dplyr)



# prepare data ---------------------

setwd("C:/Users/u0125029/Documents/3. Onderzoek/3. Post-test 2019/5. Analyse/ASSR/Data_469")

ASSRdata              = read.csv("ASSRdata_469.csv", header=TRUE, sep=",")

# only keep rows for electrodes 'All' 
ASSRdata              = subset(ASSRdata, Channel%in%c('All')) 

#strsplit only works with character arguments
ASSRdata$Recording    = as.character(ASSRdata$Recording)

# create column subject and group
ASSRdata$subject      = sapply(strsplit(ASSRdata$Recording, split='_', fixed=TRUE), function(x) (x[1]))
  #x[1]: keep the part before '_', 2: keep the part after '_'
groups                = read_excel("c1groups.xlsx")
ASSRdata$group        = groups$Groep[match(unlist(ASSRdata$subject), groups$I_code)]       
  #use groups dataframe to find correct group for each i-code

# create column stimulus type (AM/PULS)
ASSRdata$stimtype     = (str_extract(sapply(strsplit(ASSRdata$Recording, split='_', fixed=TRUE), function(x) (x[2])), "[aA-zZ]+"))
  #strsplit to keep AM4, ... + str_extract to only keep letter part

# create column frequency
ASSRdata$frequency    = ASSRdata$StimulationFrequency
ASSRdata              = subset(ASSRdata, frequency%in%c('20')) # otherwise you will also keep 19 & 21

# create column condition (stimtype + freq)
ASSRdata$condition    = sapply(strsplit(ASSRdata$Recording, split='_', fixed=TRUE), function(x) (x[2]))
#ASSRdata              = subset(ASSRdata, condition%in%c('AM20')) # keep the condition for which you want to check percentage responses

# create column testing phase (pre/post)
ASSRdata$prepost      = sapply(strsplit(ASSRdata$Recording, split='_', fixed=TRUE), function(x) (x[3]))

# only keep relevant columns
ASSRdata              = ASSRdata[c("subject", "condition", "stimtype", "frequency", "AssrHt2Statistic")]



# subset data
colnames(ASSRdata)= c("subject", "cond", "stimtype", "freq", "ht2")
head(ASSRdata)


# analysis percentage -----------------


TF <- ASSRdata$ht2 > 0.05
TF <- as.data.frame(TF)
percentage <- table(TF)
cbind(percentage,prop.table(percentage))

