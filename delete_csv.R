# ---------------------------------------------------------- #
# name              delete_CSV.R                             #
# description       used to delete .csv based on name        #
# version.string    R version 3.6.0 (2019-04-26)             #
# platform          x86_64-w64-mingw32                       #
# date created      04/02/2019                               #
# ---------------------------------------------------------- #

# folders
ROOT = 'D:/CohortI_c1project/4_Raw_data/ASSR'

# get folder structure
all_folders     = list.files(ROOT)                              
subject_folders = all_folders[startsWith(all_folders, 'i')]
conditions      = c('pre', 'post')

for (sub in subject_folders){
  for (cond in conditions){
    # glue path together
    CURRENT = paste(ROOT, sub, cond, 'eeg', sep = '/')  
    # delete files
    # list files ending with .csv, and containing AM or PULS
    csv     = list.files(path    = CURRENT,
                         pattern = "(.*)AM4(.*)469(.*)csv$|(.*)PULS(.*)469(.*)csv$")
    # copy these files to DEST
    file.remove(file.path(CURRENT, csv))
  }
}