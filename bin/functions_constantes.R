
####################
#### Constantes ####
####################

# Inspired of OkabeIto color-blind friendly palette
color_code = c(
  "CTRL" = "#009E73",
  "CTRL2" = "#009E73",
  "DON" = "#0072B2",
  "2DG" = "#E69F00",
  "AOA" = "#D30707",
  "CTRLak" = "#66FF33",
  "DONaK" = "#56B4E9",
  "2DGaK" = "#F0E442",
  "AOAaK" = "#FF8989",
  "VPA" = "#FF59B5"
)


###################
#### Functions ####
###################

# Finds last subfolder generated in the parent folder
pic_last_dir = function(parent_folder){
  dir = list.dirs(path = parent_folder, recursive = F, full.names = F)
  dir = dir[length(dir)]
  return(paste0(parent_folder,dir))
}

# Loads an RData file and returns it
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}
