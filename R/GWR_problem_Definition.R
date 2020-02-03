## ------------------------------------------------------------------------
library(sf)
library(raster)
library(fasterize)
library(DIFMR)

set.seed(12345)

f_coords <- base::cbind(c(0, 400), c(0, 300))
pol <- sf::st_as_sfc(sf::st_bbox(sf::st_multipoint(f_coords)))
field <- sf::st_sf(pol, crs = 32616)
rst <- raster(field, res = 10)
trial <- raster::aggregate(rst, c(4, 1))
trial$id <- 1:ncell(trial)
trial$pcol = 1:ncol(trial)
trial$prow = rep(1:nrow(trial), each = ncol(trial))

pts = rasterToPoints(rst, spatial = TRUE)
pts = st_as_sf(pts)
cnt = st_centroid(st_geometry(field))

pols <- rasterToPolygons(trial)
pols <- st_as_sf(pols)
pols <- get_block_ids(pols)

Treat_rates <- c(0, 50, 100, 150, 200)
pols_treat <- design_treat_graeco(pols, length(Treat_rates))
pols_treat$Treat <- Treat_rates[pols_treat$F1]

## ------------------------------------------------------------------------
trial$Treat <- pols_treat$Treat[pols_treat$id]
plot(trial$Treat)

rst$Treat <- raster::resample(trial$Treat, rst, method = "ngb")
plot(rst$Treat)

## ------------------------------------------------------------------------
trial_grd <- rasterToPoints(rst, spatial = TRUE)
gridded(trial_grd) <- TRUE

m <- gstat::vgm(psill = 1, model = "Gau",
                range = 50,
                nugget = 0.1)
g.dummy <- gstat::gstat(
  formula = z ~ 1,
  dummy = TRUE, beta = 0,
  model = m, nmax = 10
)
rst_sim <- predict(g.dummy, trial_grd, nsim = 1)
rst$Yield_Ref <- raster::scale(raster::stack(rst_sim))
rst$Yield_Ref = 1e3 * rst$Yield_Ref
plot(rst$Yield_Ref)

m <- gstat::vgm(psill = 1, model = "Gau",
                range = 100,
                nugget = 0.1)
g.dummy <- gstat::gstat(
  formula = z ~ 1,
  dummy = TRUE, beta = 0,
  model = m, nmax = 10
)
rst_sim <- predict(g.dummy, trial_grd, nsim = 1)
rst$Error = 2e2 * rnorm(ncell(rst_sim))
rst$Treat_Resp <- raster::scale(raster::stack(rst_sim))
rst$NR_gwr = 10 + 2.5 * rst$Treat_Resp
rst$NR2_gwr = -0.05
rst$NR_opt = -0.5 * rst$NR_gwr / rst$NR2_gwr

hist(rst$NR_opt[])

rst$TreatD <- rst$Treat - 100
rst$TreatD2 <- rst$TreatD**2
rst$NRD_opt = rst$TreatD - rst$NR_opt
hist(rst$NRD_opt[])

rst$Treat_Yield = rst$NR_gwr * rst$TreatD + rst$NR2_gwr * rst$TreatD
rst$Yield_Obs = 1e4 + rst$Yield_Ref + rst$Error + rst$Treat_Yield

lm0 = lm(Yield_Obs ~ poly(TreatD, 2, raw = TRUE), data.frame(rst[]))
summary(lm0)
rst$Yield_Res <- residuals(lm0)


trial_grd = rasterToPoints(rst, spatial = TRUE)
gridded(trial_grd) = TRUE
pts <- sp::coordinates(trial_grd)
dMat <- GWmodel::gw.dist(pts, pts)
gwr.fml <- as.formula(Yield_Obs ~ TreatD)

bwcv = 20
# bwcv <- GWmodel::bw.gwr(gwr.fml, trial_grd,
#                         approach = "AIC", dMat = dMat,
#                         kernel = "gaussian"
# )


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




