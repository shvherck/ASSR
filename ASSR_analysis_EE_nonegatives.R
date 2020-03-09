# --------------------------------------------------------- #
# name              ASSR_analysis_EE_nonegatives.R          #
# created by        Shauni Van Herck                        #
# description       used to analyse ASSR data (EE effect)   #
# version.string    R version 3.5.1 (2018-07-02)            #
# platform          x86_64-w64-mingw32                      #
# date created      09/03/2020                              #
# --------------------------------------------------------- #

# load libraries -------------------
# reading
library(readxl)
# reshaping
library(stringr)      # split data into numbers/letters
library(reshape)
library(dplyr)
# plotting
library(car) #qqplot
library(ggplot2)
library(gridExtra)
# correlation matrix
library(corrplot)
source("http://www.sthda.com/upload/rquery_cormat.r") # necessary to use rquery.cormat
library(PerformanceAnalytics)
library(corrr)


# prepare data ---------------------

setwd("C:/Users/u0125029/Documents/3. Onderzoek/3. Post-test 2019/5. Analyse/ASSR/R")

d              = read.csv("ASSRdata.csv", header=TRUE, sep=",")

d$AssrHt2BiasedRecordingSnrDb[d$AssrHt2BiasedRecordingSnrDb < 0] <- 0

# only keep rows for electrodes 'Right', 'Left' and 'All'
# only keep rows for frequencies 4 & 20
d              = subset(d, Channel%in%c('Right', 'Left', 'All')) 
d              = subset(d, StimulationFrequency%in%c('4', '20')) 

#strsplit only works with character arguments
d$Recording    = as.character(d$Recording)

# create column subject and group
d$subject      = sapply(strsplit(d$Recording, split='_', fixed=TRUE), function(x) (x[1]))
#x[1]: keep the part before '_', 2: keep the part after '_'

groups         = read_excel("c1groups.xlsx")
d$group        = groups$Groep[match(unlist(d$subject), groups$I_code)]       
#use groups dataframe to find correct group for each i-code

# create column stimulus type (AM/PULS)
d$stimtype     = (str_extract(sapply(strsplit(d$Recording, split='_', fixed=TRUE), function(x) (x[2])), "[aA-zZ]+"))
#strsplit to keep AM4, ... + str_extract to only keep letter part

# create column frequency
d$frequency    = d$StimulationFrequency

# create column condition (stimtype + freq)
d$condition    = sapply(strsplit(d$Recording, split='_', fixed=TRUE), function(x) (x[2]))

# create column testing phase (pre/post)
d$prepost      = sapply(strsplit(d$Recording, split='_', fixed=TRUE), function(x) (x[3]))

# only keep relevant columns
d              = d[c("subject", "group", "condition", "stimtype", "frequency", "prepost", "Channel", "AssrHt2BiasedRecordingSnrDb")]

EE             = read_excel("EE_exposure.xlsx")
d$EE_hours     = EE$ExposureHours[match(unlist(d$subject), EE$iCode)]   
d$EE_accuracy  = EE$Accuracy[match(unlist(d$subject), EE$iCode)]
#use EE dataframe to find correct GG factors for each i-code

# subset data
d_all          = subset(d, Channel%in%c('All'))
d_all          = subset(d, group%in%c('GG_EE'))
colnames(d_all)= c("subject", "group", "cond", "stimtype", "freq", "prepost", "channel", "snr", "hours", "accuracy")
head(d_all)

d_all_melt           = melt(d_all, id.var=c("subject", "group", "cond", "stimtype", "freq", "prepost", "hours", "accuracy"), measure.var=c("snr"))
d_all_wide           = cast(subject + group + cond + stimtype + freq + hours + accuracy ~ prepost , data=d_all_melt, mean)
head(d_all_wide)
d_all_wide$diff      = d_all_wide$post - d_all_wide$pre
colnames(d_all_wide) = c("subject", "group", "cond", "stimtype", "freq", "hours", "accuracy", "post", "pre", "diff") 

d_4                  = subset(d_all_wide, freq%in%c('4'))
d_20                 = subset(d_all_wide, freq%in%c('20'))

corr4                = as.data.frame(d_4[c("hours", "accuracy", "diff")]) # needs to be data frame so corr recognizes column names to use
corr20               = as.data.frame(d_20[c("hours", "accuracy", "diff")])

corr4$accuracy       = as.numeric(corr4$accuracy)   # r defines these columns as character but they need to be numeric for calculations
corr20$accuracy      = as.numeric(corr20$accuracy)

# plots ------------------

setwd("C:/Users/u0125029/Documents/3. Onderzoek/3. Post-test 2019/5. Analyse/ASSR/R/Plots/corr_EE")

## 4 Hz

tiff("corr4_hoursEE.tiff",width=400,height=400)                 #save figure in wd
plot(corr4$hours, corr4$diff,  main="posttest - pretest (4 Hz) vs. hours played EE", xlab="Hours played EE", ylab="Posttest - pretest (4 Hz) (snr dB)")
abline(lm(diff ~ hours, data=corr4))
dev.off()

tiff("corr4_accuracyEE.tiff",width=400,height=400)
plot(corr4$accuracy, corr4$diff,  main="posttest - pretest (4 Hz) vs. accuracy EE", xlab="Accuracy EE", ylab="Posttest - pretest (4 Hz) (snr dB")
abline(lm(diff ~ accuracy, data=corr4))
dev.off()

## 20 Hz

tiff("corr20_hoursEE.tiff",width=400,height=400)                 #save figure in wd
plot(corr20$hours, corr20$diff,  main="posttest - pretest (20 Hz) vs. hours played EE", xlab="Hours played EE", ylab="Posttest - pretest (20 Hz) (snr dB)")
abline(lm(diff ~ hours, data=corr20))
dev.off()

tiff("corr20_accuracyEE.tiff",width=400,height=400)
plot(corr20$accuracy, corr20$diff,  main="posttest - pretest (20 Hz) vs. accuracy EE", xlab="Accuracy EE", ylab="Posttest - pretest (20 Hz) (snr dB")
abline(lm(diff ~ accuracy, data=corr20))
dev.off()


# analyses ----------------

## 4 Hz

cor.test(corr4$diff, corr4$hours, method="pearson")
cor.test(corr4$diff, corr4$accuracy, method="pearson") 

round(cor(corr4, method = c("pearson"), use="complete.obs"), 2)

rquery.cormat(corr4)

chart.Correlation(corr4, histogram=TRUE, header=TRUE)

## 20 Hz

cor.test(corr20$diff, corr20$hours, method="pearson")
cor.test(corr20$diff, corr20$accuracy, method="pearson") 

round(cor(corr20, method = c("pearson"), use="complete.obs"), 2)

rquery.cormat(corr20)

chart.Correlation(corr20, histogram=TRUE, header=TRUE)
