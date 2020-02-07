## ------------------------------------------------------------------------
library(sf)
library(tmap)
library(raster)
library(tibble)
library(dplyr)
library(broom)
library(tidyr)
library(Metrics)

## ------------------------------------------------------------------------

st_interpolate = function(x){
  x[] = as.numeric(as.factor(x[]))
  w <- matrix(1, 3, 3)
  while (sum(is.na(x[])) > 0) {
    x <- focal(x, w = w, fun = modal, NAonly = TRUE, na.rm = TRUE,  pad = TRUE)
  }
  return(x)
}


set.seed(12345)

trial <- readRDS('data/Trial_Design.rds')

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
rst_b0 <- 5000 + 500 * rst_sim[[1]]
rst_b1g <- 20
rst_b2g <- -0.1
rst_b1 <- rst_b1g + 3 * rst_sim[[2]]
rst_optr = -0.5 * rst_b1 / rst_b2g
rst_err <- 100 * rnorm(ncell(rst_sim))

rst_optrd <- disaggregate(aggregate(rst_optr, 2), 2)

plot(rst_optr)
plot(rst_optrd)

plot(rst_optrd, rst_optr, asp = 1)
cor(rst_optrd[], rst_optr[])
rmse(rst_optrd[], rst_optr[])



inp_rate = 100 + 50 * rst$size15_p100
# inp_rate[is.na(inp_rate)] = 100
rst_yield =  rst_b0 + rst_b1 * inp_rate + rst_b2g * inp_rate**2 + rst_err

plot(rst_yield)
plot(inp_rate, rst_yield)
plot(inp_rate)


grdf = stack(inp_rate, rst_yield, rst$size15_rep)
names(grdf) = c('inp_rate', 'yield', 'block')
grdf = mask(grdf, inp_rate)
grdf$block = st_interpolate(grdf$block)
plot(grdf)


trial = st_as_sf(rasterToPolygons(mask(grdf, inp_rate)))
inp_rate_1_2 = poly(trial$inp_rate, 2)
inp_rate_1_2 = poly(trial$inp_rate, 2, raw = TRUE)
trial$inp_rate_1 = inp_rate_1_2[,1]
trial$inp_rate_2 = inp_rate_1_2[,2]

trial$block = as.factor(trial$block)
lm.fml = yield ~ inp_rate_1 + inp_rate_2 + inp_rate_1 * block
lm0 = lm(lm.fml, trial)
summary(lm0)

trial$yield_pred = predict(lm0)
plot(yield ~ yield_pred, trial, asp = 1)

# dm = model.matrix(lm.fml, trial)
# cm = matrix(coef(lm0))

#Organize predicted coefficients for each block:
lmdf = tidy(lm0) %>% 
  separate(term, into = c('term', 'block'), sep = ':?block', convert = TRUE)

lmdf$block[1:3] = 1
lmdf$term[1] = ''
lmdf$term[lmdf$term == ''] = 'intercept'

lmdf = lmdf %>% select(1:3) %>%
  spread(key = 'term', value = 'estimate') %>%
  fill(inp_rate_2) %>% mutate(block = as.factor(block)) %>% 
  rename(block = block, b0 = intercept, b1 = inp_rate_1, b2 = inp_rate_2)
lmdf$b0[-1] = lmdf$b0[1] + lmdf$b0[-1]
lmdf$b1[-1] = lmdf$b1[1] + lmdf$b1[-1]

trial_pred = left_join(trial, lmdf, by = 'block') %>%
  mutate(yield_calc = b0 + inp_rate_1 * b1 + inp_rate_2 * b2)

plot(yield_calc ~ yield_pred, trial_pred, asp = 1)

trial_pred = trial_pred %>% mutate(optr = -0.5 * b1 / b2)
hist(trial_pred$optr)

trial_block = trial_pred %>% group_by(block) %>% summarise_all(mean)

brks = seq(3000,6000,500)
tm_shape(rst_b0) + tm_raster(breaks = brks, legend.show = FALSE) +
  tm_shape(trial_block) + tm_polygons('b0', breaks = brks) +
  tm_layout(legend.outside=TRUE)

brks = seq(5,30,5)
tm_shape(rst_b1) + tm_raster(breaks = brks, legend.show = FALSE) +
  tm_shape(trial_block) + tm_polygons('b1', breaks = brks) +
  tm_layout(legend.outside=TRUE)


brks = seq(50,150,10)
tm_shape(rst_optr) + tm_raster(breaks = brks, legend.show = FALSE) +
  tm_shape(trial_block) + tm_polygons('optr', breaks = brks) +
  tm_layout(legend.outside=TRUE)

grd_pts = as_tibble(rasterToPoints(grdf$block))
grd_pts$block = as.factor(grd_pts$block)
grd_pts = left_join(grd_pts, trial_block, by = 'block')
rst_optrp = rst_optr
rst_optrp[] = grd_pts$optr

brks = seq(50,150,10)
tm_shape(rst_optrp) + tm_raster(breaks = brks, legend.show = FALSE) +
  tm_shape(trial_block) + tm_polygons('optr', breaks = brks) +
  tm_layout(legend.outside=TRUE)
# Evaluate the betas or the maximum rate? Second, plus economical aspects.

plot(rst_optrp, rst_optr, asp = 1)
cor(rst_optrp[], rst_optr[])
rmse(rst_optrp[], rst_optr[])

