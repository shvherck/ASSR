# --------------------------------------------------------- #
# name              Append_datasets.R                       #
# description       used to combine ASSR csv files in 1 csv #
# version.string    R version 3.5.1 (2018-07-02)            #
# platform          x86_64-w64-mingw32                      #
# date created      07/10/2019                              #
# --------------------------------------------------------- #

library(plyr)
library(readr)

setwd(dir = "C:/Users/u0125029/Documents/3. Onderzoek/3. Post-test 2019/5. Analyse/ASSR/Data")

list = list.files(path = "C:/Users/u0125029/Documents/3. Onderzoek/3. Post-test 2019/5. Analyse/ASSR/Data")

dat_csv = ldply(list, read_csv)               # for each element of a list, apply a function and then combine the results into a data frame
dat_csv

write.csv(dat_csv, file="ASSRdata.csv")


