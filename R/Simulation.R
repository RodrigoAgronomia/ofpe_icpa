## ------------------------------------------------------------------------
library(sf)
library(raster)
library(tibble)
library(dplyr)
library(mgwrsar)

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

inp_rate = 100 + 50 * rst[[20]]
inp_rate[is.na(inp_rate)] = 0
rst_yield =  rst_b0 + rst_b1 * inp_rate + rst_b2g * inp_rate**2 + rst_err

plot(rst_yield)
plot(inp_rate, rst_yield)


grdf = aggregate(stack(inp_rate, rst_yield), 10)
trial_grd = rasterToPoints(grdf, spatial = TRUE)
gridded(trial_grd) = TRUE
pts <- sp::coordinates(trial_grd)
dMat <- GWmodel::gw.dist(pts, pts)
gwr.fml <- as.formula(Yield_Obs ~ TreatD)

bwcv = 20


gwr.model <- GWmodel::gwr.basic(gwr.fml, trial_grd,
                                bw = bwcv, dMat = dMat,
                                kernel = "gaussian"
)
print(gwr.model)


gwr_r <- gwr.model$SDF
gridded(gwr_r) = TRUE
gwr_rst = stack(gwr_r)
gwr_rst$Treat = gwr_rst$TreatD
plot(gwr_rst)

plot(gwr_rst$Treat[], rst$NR_gwr[], asp = 1)
cor(gwr_rst$Treat[], rst$NR_gwr[])
plot(rst$NR_gwr)
plot(gwr_rst$Treat)

plot(gwr_rst$Intercept)
plot(rst$Yield_Ref)

plot(gwr_rst$Intercept[], rst$Yield_Ref[], asp = 1)
cor(gwr_rst$Intercept[], rst$Yield_Ref[])




