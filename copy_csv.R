# ---------------------------------------------------------- #
# name              copy_CSV.R                               #
# description       used to copy .csv based on name and cond #
# version.string    R version 3.6.0 (2019-04-26)             #
# platform          x86_64-w64-mingw32                       #
# date created      07/10/2019                               #
# ---------------------------------------------------------- #

# folders
ROOT = 'D:/CohortI_c1project/4_Raw_data/ASSR'
DEST = 'C:/Users/u0125029/Documents/3. Onderzoek/3. Post-test 2019/5. Analyse/ASSR/Data'

# get folder structure
all_folders     = list.files(ROOT)                              
subject_folders = all_folders[startsWith(all_folders, 'i')]
conditions      = c('pre', 'post')

for (sub in subject_folders){
  for (cond in conditions){
    # glue path together
    CURRENT = paste(ROOT, sub, cond, 'eeg', sep = '/')  
    # list files ending with .csv, and containing AM or PULS
    csv     = list.files(path    = CURRENT,
                         pattern = "(.*)AM(.*)csv$|(.*)PULS(.*)csv$")
    # copy these files to DEST
    file.copy(file.path(CURRENT, csv), 
              DEST)
  }
}
