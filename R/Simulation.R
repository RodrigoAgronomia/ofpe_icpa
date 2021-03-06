## ------------------------------------------------------------------------
library(sf)
library(raster)
library(tibble)
library(dplyr)
library(GWmodel)

set.seed(12345)

rst <- stack('data/Trial_Design.tif')
trial <- aggregate(raster(rst), 10)


## ------------------------------------------------------------------------
trial_grd <- rasterToPoints(trial, spatial = TRUE)
gridded(trial_grd) <- TRUE

m <- gstat::vgm(psill = 1, model = "Gau",
                range = 100,
                nugget = 0.5)
g.dummy <- gstat::gstat(
  formula = z ~ 1,
  dummy = TRUE, beta = 0,
  model = m, nmax = 10
)
# rst_sim <- predict(g.dummy, trial_grd, nsim = 200)
# rst_sim <- scale(stack(rst_sim))
# rst_sim <- disaggregate(rst_sim, 10)
# rst_b0 <- rst_sim[[1:100]]
# rst_b1 <- rst_sim[[101:200]]

rst_sim <- predict(g.dummy, trial_grd, nsim = 2)
rst_sim <- scale(stack(rst_sim))
rst_sim <- disaggregate(rst_sim, 10)
rst_b0 <- 5000 + 500 * rst_sim[[1]]
rst_b1g <- 20
rst_b2g <- -0.1
rst_b1 <- rst_b1g + 3 * rst_sim[[2]]
rst_optr = -0.5 * (rst_b1 + rst_b1g) / rst_b2g
rst_err <- 100 * rnorm(ncell(rst_sim))

inp_rate = 100 + 50 * rst[[12]]
# inp_rate[is.na(inp_rate)] = 100
rst_yield =  rst_b0 + rst_b1 * inp_rate + rst_b2g * inp_rate**2 + rst_err

plot(rst_yield)
plot(inp_rate, rst_yield)
plot(inp_rate)


grdf = aggregate(stack(inp_rate, rst_yield), 10)
names(grdf) = c('inp_rate', 'yield')
trial = st_as_sf(rasterToPolygons(grdf))
pred_coords = rasterToPoints(is.na(grdf), spatial = TRUE)
gridded(pred_coords) = TRUE
trial_sp <- st_set_geometry(trial, st_centroid(st_geometry(trial)))
trial_sp <- as(trial_sp, "Spatial")
pts <- coordinates(trial_sp)
gwr.fml <- as.formula(yield~poly(inp_rate, 2))
bwcv <- 20

gwr.model <- gwr.basic(gwr.fml, trial_sp, bw = bwcv,
                       regression.points = pred_coords,
                       adaptive = TRUE, kernel = "gaussian")


rst_b1_pred = stack(gwr.model$SDF)


gwr_r <- gwr.model$SDF$Intercept
names(gwr_r) <- paste0(c('Intercept', 'inp_rate', 'inp_rate2'), "_gwr")
trial_gwr <- bind_cols(trial, gwr_r)





