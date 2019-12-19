# -----------------------------------------------------------#
# name              epochcount.R                             #
# description       used to check epochcount < 448           #
# version.string    R version 3.5.1 (2018-07-02)             #
# platform          x86_64-w64-mingw32                       #
# date created      08/10/2019                               #
# -----------------------------------------------------------#


setwd(dir = "C:/Users/u0125029/Documents/3. Onderzoek/3. Post-test 2019/5. Analyse/ASSR/Data")

list = list.files(path = "C:/Users/u0125029/Documents/3. Onderzoek/3. Post-test 2019/5. Analyse/ASSR/Data")


for(i in 1:length(list)){
  csv = read.csv(list[i])
  csv$notenoughepochs = csv$EpochCount < 448
  if (csv$notenoughepochs[1] == "TRUE"){
    print(list[i])
  }
}