## ------------------------------------------------------------------------
library(sf)
library(raster)
library(gstat)
library(spmoran)
library(tibble)
library(dplyr)
library(tidyr)
library(tmap)

## ------------------------------------------------------------------------
trial_rst <- readRDS("data/Trial_Design.rds")
rst <- raster(trial_rst)

set.seed(12345)

# Fixed set of parameters to control the response function:
cotton_price = 0.4
nitrogen_cost = 1.0
nitrogen_ratio = nitrogen_cost/cotton_price
sq_rate = 100
eonr = 100
sq_yield = 4500
yield_sd = 500
eonr_sd = 25
rst_b2g = -0.02
test_rates = sq_rate + 50 * c(-1,-0.5,0,0.5,1)

n_rate <- sq_rate + 50 * trial_rst$size15_p050
n_rate[is.na(n_rate[])] = sq_rate
plot(n_rate)

brks = seq(0, 200, length.out = 21)
lbls = rep('', length(brks))
crit = seq(1, length(lbls), 5)
lbls[crit] = as.character(brks)[crit]
tm = tm_shape(n_rate) +
  tm_raster(style='cont', labels = lbls,
            breaks = brks,
            legend.is.portrait = TRUE, title = 'N rate (kg/ha)') +
  tm_layout(legend.outside = TRUE,
            legend.text.size = 1,
            legend.title.size = 2,
            legend.outside.size = c(0.15, 0.15),
            legend.outside.position = 'right')

tmap_save(tm, 'figures/Nrate_map.png', width = 9, height = 6)


# This creates the spatial distribution of two random variables:
trial_grd <- rasterToPoints(rst, spatial = TRUE)
gridded(trial_grd) <- TRUE
m <- vgm(psill = 1, model = "Gau", range = 100, nugget = 0.1)
g.dummy <- gstat(formula = z ~ 1, dummy = TRUE, beta = 0, model = m, nmax = 10)
rst_sim <- stack(predict(g.dummy, trial_grd, nsim = 2))

# This adjust the simulated values to the desirade scales:
sq_yield_rst <- sq_yield + yield_sd * rst_sim$sim1
eonr_rst <- eonr + eonr_sd * rst_sim$sim2
eonr_rst <- clamp(eonr_rst, 0, 200)


brks = seq(2000, 7000, length.out = 21)
lbls = rep('', length(brks))
crit = seq(1, length(lbls), 5)
lbls[crit] = as.character(brks)[crit]
tm = tm_shape(sq_yield_rst) +
  tm_raster(style='cont', labels = lbls,
            breaks = brks, palette = rev(terrain.colors(50)),
            legend.is.portrait = TRUE, title = 'Yield (kg/ha)') +
  tm_layout(legend.outside = TRUE,
            legend.text.size = 1,
            legend.title.size = 2,
            legend.outside.size = c(0.15, 0.15),
            legend.outside.position = 'right')

tmap_save(tm, 'figures/Yield_map.png', width = 9, height = 6)


brks = seq(0, 200, 10)
lbls = rep('', length(brks))
crit = seq(1, length(lbls), 5)
lbls[crit] = as.character(brks)[crit]
tm = tm_shape(eonr_rst) +
  tm_raster(style='cont', labels = lbls, breaks = brks,
            legend.is.portrait = TRUE, title = 'EONR (kg/ha)') +
  tm_layout(legend.outside = TRUE,
            legend.text.size = 1,
            legend.title.size = 2,
            legend.outside.size = c(0.15, 0.15),
            legend.outside.position = 'right')

tmap_save(tm, 'figures/EONR_map.png', width = 9, height = 6)

# This is how I use the yield reponse function to get beta0 and beta1,
#  given the yield at the status quo rate and the EONR:
rst_b1 = -2 * rst_b2g * eonr_rst + nitrogen_ratio
rst_b0 = sq_yield_rst - (rst_b1 * sq_rate + rst_b2g * sq_rate **2) 
yield_obs = rst_b0 + rst_b1 * n_rate + rst_b2g * n_rate **2



brks = seq(2000, 7000, length.out = 21)
lbls = rep('', length(brks))
crit = seq(1, length(lbls), 5)
lbls[crit] = as.character(brks)[crit]
tm = tm_shape(yield_obs) +
  tm_raster(style='cont', labels = lbls,
            breaks = brks, palette = rev(terrain.colors(50)),
            legend.is.portrait = TRUE, title = 'Yield (kg/ha)') +
  tm_layout(legend.outside = TRUE,
            legend.text.size = 1,
            legend.title.size = 2,
            legend.outside.size = c(0.15, 0.15),
            legend.outside.position = 'right')

tmap_save(tm, 'figures/Yield_obs_map.png', width = 9, height = 6)


