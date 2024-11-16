#############################
# Data
#############################
files <- dir()
files <- files[grep(pattern = "R_CODE_microglial markers_invivo_adjust.txt", x = files)]
data_original <- read.delim(file = files, header = TRUE, sep = "\t", dec = ",")
source(file = paste(dirname(getwd()), "MetaA_GregorAdapt_202403.R", sep = "/"))

