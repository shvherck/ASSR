# ------------------------------------------------------------------------------ #
# name              mismatch_filename_metadata.R                                 #
# description       used to check mismatch between filename and metadata in file #
# version.string    R version 3.5.1 (2018-07-02)                                 #
# platform          x86_64-w64-mingw32                                           #
# date created      08/10/2019                                                   #
# ------------------------------------------------------------------------------ #

setwd(dir = "C:/Users/u0125029/Documents/3. Onderzoek/3. Post-test 2019/5. Analyse/ASSR/Data")

list = list.files(path = "C:/Users/u0125029/Documents/3. Onderzoek/3. Post-test 2019/5. Analyse/ASSR/Data")


for(i in 1:length(list)){
  csv = read.csv(list[i])
  csv$filename = sapply(strsplit(list[i], split='.', fixed=TRUE), function(x) (x[1]))
  csv$match = csv$Recording == csv$filename
  if (csv$match[1] == "FALSE"){
    print(list[i])
  }
}

