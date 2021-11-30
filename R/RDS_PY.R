## ------------------------------------------------------------------------
library(raster)
library(reticulate)

np <- import("numpy")

trial <- readRDS('data/Trial_Design.rds')
wwm = aperm(as.array(trial), c(3,1,2))

wfile = './data/Trial_Design.npy'
np$save(wfile, wwm)

wwnames = as.array(names(trial))

names_file = './data/Trial_Design_names.npy'
np$save(names_file, wwnames)


wwcoords = as.array(coordinates(trial))

coords_file = './data/Trial_Design_coords.npy'
np$save(coords_file, wwcoords)



