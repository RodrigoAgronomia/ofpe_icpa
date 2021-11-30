## ------------------------------------------------------------------------
library(sf)
library(sp)
library(raster)
library(tibble)
library(dplyr)
library(gstat)
library(GWmodel)
library(reticulate)
np <- import("numpy")

set.seed(12345)

trial <- readRDS('data/Trial_Design.rds')
rst <- raster(trial)

wwr = stack(wwf)
wwm = as.array(wwr)

wfile = './data/weigths.npy'
np$save(wfile, wwm)


## ------------------------------------------------------------------------
trial_grd <- rasterToPoints(trial, spatial = TRUE)
gridded(trial_grd) <- TRUE
pts <- coordinates(trial_grd)


bw = 100
wwf = list()
i = 1
for(i in 1:nrow(pts)){
  print(i)
  dist.vi <- gw.dist(pts, pts[i,])
  rst$W.i <- gw.weight(dist.vi, bw, kernel = "gaussian", adaptive = TRUE)
  wwf[[paste0('W', i)]] = round(255 * rst$W.i)
}
wwr = stack(wwf)
wwm = as.array(wwr)

wfile = './data/weigths.npy'
np$save(wfile, wwm)


i = 1
for(bw in bws){
  rst$W.i <- gw.weight(dist.vi, bw, kernel = "gaussian", adaptive = TRUE)
  wwf[[paste0('W', bw)]] = round(255 * rst$W.i)
}
wwr = stack(wwf)
wwm = as.array(wwr)

wfile = './data/weigths_zoom.npy'
np$save(wfile, wwm)







