## ------------------------------------------------------------------------
library(sf)
library(raster)
library(gstat)
library(spmoran)
library(reticulate)
np <- import("numpy")

## ------------------------------------------------------------------------

set.seed(12345)

trial_rst <- readRDS('data/Trial_Design.rds')

## ------------------------------------------------------------------------

coords <- coordinates(trial_rst)
meig <- meigen_f(coords=coords, model = 'gau', enum = 100)
mm = matrix(meig$sf)
dim(mm) = c(300, 200, 100)
sfile = './data/Trial_meig.npy'
np$save(sfile, mm)


## ------------------------------------------------------------------------


trial_grd <- rasterToPoints(trial_rst, spatial = TRUE)
gridded(trial_grd) <- TRUE

m <- vgm(psill = 1, model = "Gau", range = 100, nugget = 0.1)
g.dummy <- gstat(formula = z ~ 1, dummy = TRUE, beta = 0, model = m, nmax = 10)

rst_sim <- predict(g.dummy, trial_grd, nsim = 200)

rst_sim <- scale(stack(rst_sim))
saveRDS(rst_sim, 'data/GSTAT_sim.rds')

wwm = aperm(as.array(rst_sim), c(3,1,2))

sfile = './data/Trial_sim.npy'
np$save(sfile, wwm)



