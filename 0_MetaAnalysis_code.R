#############################
# Data
#############################
files <- dir()
files <- files[grep(pattern = "C:\\Users\\Luis\\Desktop\\Análise_V2_microglial markers\\R_CODE_microglial_markers_invitro.txt", x = files)]
data_original <- read.delim(file = files, header = TRUE, sep = "\t", dec = ",")
source(file = paste(dirname(getwd()), "C:\\Users\\Luis\\Desktop\\MetaA_GregorAdapt_202403.R", sep = "/"))