## ------------------------------------------------------------------------
library(sf)
library(raster)
library(fasterize)
library(agricolae)
library(tibble)
library(dplyr)

set.seed(12345)

f_coords <- cbind(c(0, 900), c(0, 600))
pol <- st_as_sfc(st_bbox(st_multipoint(f_coords)))
field <- st_sf(pol, crs = 32616)
rst <- raster(field, res = 1)

trial <- aggregate(rst, 3)

n_treats = 5
treat <- seq(1, n_treats)
n_blocks = round(ncell(trial)/(n_treats ** 2))
n_rows = nrow(trial)/n_treats
n_cols = ncol(trial)/n_treats

tidx = rep(1:n_treats, 5) + rep(1:n_treats, each = 5)
tidx = tidx %% n_treats + 1
tidx = matrix(tidx, ncol = n_treats)
new_idx = lapply(1:n_blocks, function(x) tidx[sample(treat), sample(treat)])
new_idx = unlist(new_idx)
dim(new_idx) = c(n_treats, n_treats, n_rows, n_cols)
new_idx = aperm(new_idx, c(1, 3, 2, 4))
dim(new_idx) = c(n_treats * n_rows, n_treats * n_cols)
trial$treat = new_idx


# trial$id <- 1:ncell(trial)
# trial$pcol = 1:ncol(trial)
# trial$prow = rep(1:nrow(trial), each = ncol(trial))
# trial$block = 1e6 * ceiling(trial$pcol/n_treats) + ceiling(trial$prow/n_treats)
# trial$block = as.numeric(as.factor(trial$block[]))
# trial$col = (trial$pcol - 1) %% n_treats + 1
# trial$row = (trial$prow - 1) %% n_treats + 1
# 
# 
# rdf <- as_tibble(cbind(row = trial$row[], col = trial$col[], block = trial$block[]))
# df <- as_tibble(design.lsd(treat)$book) %>% mutate_all(as.integer)
# row <- lapply(1:n_blocks, function(x) sample(treat)[df$row]) %>% unlist()
# col <- lapply(1:n_blocks, function(x) sample(treat)[df$col]) %>% unlist()
# block <- rep(1:n_blocks, each = n_treats ** 2)
# tdf <- as_tibble(cbind(row, col, block)) %>% 
#   left_join(df, by = c("row", "col")) %>% 
#   right_join(rdf, by = c("row", "col", "block"))
#   
# 
# trial$treat <- tdf$treat
# plot(trial$treat)


rst$treat <- disaggregate(trial$treat, 3)
plot(rst$treat)

writeRaster(rst, 'data/Trial_Design.tif', overwrite = TRUE, datatype = 'INT1U')

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




