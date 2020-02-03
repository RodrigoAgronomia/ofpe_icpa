## ------------------------------------------------------------------------
library(sf)
library(raster)
library(fasterize)
library(agricolae)
library(tibble)
library(dplyr)

set.seed(12345)

tidx = rep(1:n_treats, 5) + rep(1:n_treats, each = 5)
tidx = tidx %% n_treats + 1
tidx = matrix(tidx, ncol = n_treats)


f_coords <- cbind(c(0, 900), c(0, 600))
pol <- st_as_sfc(st_bbox(st_multipoint(f_coords)))
field <- st_sf(pol, crs = 32616)
rst <- raster(field, res = 1)

plot_sizes = c(3, 10, 30, 60)
n_treats = 5
treat <- seq(1, n_treats)

plot_size = 10
for( plot_size in plot_sizes){
  trial <- aggregate(raster(rst), plot_size)
  trial$id <- 1:ncell(trial)
  trial$pcol = 1:ncol(trial)
  trial$rep_col = ceiling(trial$pcol/n_treats)
  trial$prow = rep(1:nrow(trial), each = ncol(trial))
  trial$rep_row = ceiling(trial$prow/n_treats)
  trial$rep = 1e6 * trial$prow + trial$rep_col
  trial$rep = as.numeric(as.factor(trial$rep[]))
  trial$block = 1e6 * trial$rep_row + trial$rep_col
  trial$block = as.numeric(as.factor(trial$block[]))
  
  n_rows = max(trial$rep_row[])
  n_cols = max(trial$rep_col[])
  n_blocks = max(trial$block[])
  
  new_idx = lapply(1:n_blocks, function(x) tidx[sample(treat), sample(treat)])
  new_idx = unlist(new_idx)
  dim(new_idx) = c(n_treats, n_treats, n_rows, n_cols)
  new_idx = aperm(new_idx, c(1, 3, 2, 4))
  dim(new_idx) = c(n_treats * n_rows, n_treats * n_cols)
  trial$treat = new_idx
  
  s_pct = c(1, 5, 10, 50, 100)
  
  
  t100 = trial$pcol > 0
  
  t50 = xor(trial$rep_col %% 2, trial$prow %% 2)
  plot(t50)
  
  t10 = (trial$rep_col + 1) %% 2 == 0 & (trial$prow + 2) %% 5 == 0
  plot(t10)
  sum(t10[])/ncell(t10)
  
  t05 = (trial$rep_col + 2) %% 4 == 0 & (trial$prow + 2) %% 5 == 0
  plot(t05)
  sum(t05[])/ncell(t05)
  
  t01 = (trial$rep_col + 5) %% 10 == 0 & (trial$prow + 5) %% 10 == 0
  if(!sum(t01[])){
    t01 = trial$rep_col == 2 & (trial$prow + 5) %% 10 == 0
  }
  plot(t01)
  sum(t01[])/ncell(t01)
  
  
  
  n_reps = max(trial$rep[])
  n_samples = ceiling(s_pct/100 * n_reps)
  
  block_sample = sample(1:n_blocks, n_samples[3])
  crit = trial$block[] %in% block_sample
  
  rep_sample = sample(1:n_reps, n_samples[3])
  crit = trial$rep[] %in% rep_sample
  trial$streat = trial$treat
  trial$streat[!crit] = 0
  
  plot(trial$streat)
  
  
  new_trial = disaggregate(trial$treat, plot_size)
  names(new_trial) = paste0('treat',plot_size)
  rst = addLayer(rst, new_trial)
}
plot(rst)




# for(s in n_samples){
#   rep_sample = sample(1:n_reps, s)
# }

plot(trial$treat)


writeRaster(rst, 'data/Trial_Design.tif', overwrite = TRUE, datatype = 'INT1U')

## ------------------------------------------------------------------------
# trial_grd <- rasterToPoints(rst, spatial = TRUE)
# gridded(trial_grd) <- TRUE
# 
# m <- gstat::vgm(psill = 1, model = "Gau",
#                 range = 50,
#                 nugget = 0.1)
# g.dummy <- gstat::gstat(
#   formula = z ~ 1,
#   dummy = TRUE, beta = 0,
#   model = m, nmax = 10
# )
# rst_sim <- predict(g.dummy, trial_grd, nsim = 1)
# rst$Yield_Ref <- raster::scale(raster::stack(rst_sim))
# rst$Yield_Ref = 1e3 * rst$Yield_Ref
# plot(rst$Yield_Ref)
# 
# m <- gstat::vgm(psill = 1, model = "Gau",
#                 range = 100,
#                 nugget = 0.1)
# g.dummy <- gstat::gstat(
#   formula = z ~ 1,
#   dummy = TRUE, beta = 0,
#   model = m, nmax = 10
# )
# rst_sim <- predict(g.dummy, trial_grd, nsim = 1)
# rst$Error = 2e2 * rnorm(ncell(rst_sim))
# rst$Treat_Resp <- raster::scale(raster::stack(rst_sim))
# rst$NR_gwr = 10 + 2.5 * rst$Treat_Resp
# rst$NR2_gwr = -0.05
# rst$NR_opt = -0.5 * rst$NR_gwr / rst$NR2_gwr
# 
# hist(rst$NR_opt[])
# 
# rst$TreatD <- rst$Treat - 100
# rst$TreatD2 <- rst$TreatD**2
# rst$NRD_opt = rst$TreatD - rst$NR_opt
# hist(rst$NRD_opt[])
# 
# rst$Treat_Yield = rst$NR_gwr * rst$TreatD + rst$NR2_gwr * rst$TreatD
# rst$Yield_Obs = 1e4 + rst$Yield_Ref + rst$Error + rst$Treat_Yield
# 
# lm0 = lm(Yield_Obs ~ poly(TreatD, 2, raw = TRUE), data.frame(rst[]))
# summary(lm0)
# rst$Yield_Res <- residuals(lm0)
# 
# 
# trial_grd = rasterToPoints(rst, spatial = TRUE)
# gridded(trial_grd) = TRUE
# pts <- sp::coordinates(trial_grd)
# dMat <- GWmodel::gw.dist(pts, pts)
# gwr.fml <- as.formula(Yield_Obs ~ TreatD)
# 
# bwcv = 20
# # bwcv <- GWmodel::bw.gwr(gwr.fml, trial_grd,
# #                         approach = "AIC", dMat = dMat,
# #                         kernel = "gaussian"
# # )
# 
# 
# gwr.model <- GWmodel::gwr.basic(gwr.fml, trial_grd,
#                                 bw = bwcv, dMat = dMat,
#                                 kernel = "gaussian"
# )
# print(gwr.model)
# 
# 
# gwr_r <- gwr.model$SDF
# gridded(gwr_r) = TRUE
# gwr_rst = stack(gwr_r)
# gwr_rst$Treat = gwr_rst$TreatD
# plot(gwr_rst)
# 
# plot(gwr_rst$Treat[], rst$NR_gwr[], asp = 1)
# cor(gwr_rst$Treat[], rst$NR_gwr[])
# plot(rst$NR_gwr)
# plot(gwr_rst$Treat)
# 
# plot(gwr_rst$Intercept)
# plot(rst$Yield_Ref)
# 
# plot(gwr_rst$Intercept[], rst$Yield_Ref[], asp = 1)
# cor(gwr_rst$Intercept[], rst$Yield_Ref[])




