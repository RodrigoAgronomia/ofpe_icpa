## ------------------------------------------------------------------------
library(sf)
library(tmap)
library(raster)
library(tibble)
library(dplyr)
library(broom)
library(tidyr)
library(gstat)
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

cotton_price = 1.25
nitrogen_cost = 1.0
nitrogen_ratio = nitrogen_cost/cotton_price

set.seed(12345)

trial <- readRDS('data/Trial_Design.rds')

## ------------------------------------------------------------------------
trial_grd <- rasterToPoints(trial, spatial = TRUE)
gridded(trial_grd) <- TRUE

m <- vgm(psill = 1, model = "Gau", range = 100, nugget = 0.1)
g.dummy <- gstat(formula = z ~ 1, dummy = TRUE, beta = 0, model = m, nmax = 10)

sq_rate <- 100
sq_yield <- 1000

rst_sim <- predict(g.dummy, trial_grd, nsim = 2)
rst_sim <- scale(stack(rst_sim))
rst_b1g <- nitrogen_ratio
rst_b2g <- -0.05
rst_b1 <- rst_b1g + 2 * rst_sim[[2]]
rst_optr = sq_rate -0.5 * (rst_b1 - nitrogen_ratio) / rst_b2g
hist(rst_optr[])

rst_b0 <- 1000 + 100 * rst_sim[[1]]
rst_b0 = 1000 - 100 * rst_sim[[1]] - rst_b1 * sq_rate - rst_b2g * sq_rate**2

nrate = seq(0,200)
ynr = rst_b1g * (nrate - sq_rate) + rst_b2g * (nrate - sq_rate)**2
plot(ynr ~ nrate, ylim = c(-1000, 1000))

ynr0 =  min(rst_b1[]) * (nrate - sq_rate) + rst_b2g * (nrate - sq_rate)**2
lines(ynr0 ~ nrate)

ynr0 =  max(rst_b1[]) * (nrate - sq_rate) + rst_b2g * (nrate - sq_rate)**2
lines(ynr0 ~ nrate)


inp_rate = sq_rate + 25 * trial$size15_p100

rst_b0 = rst_yield_ref - rst_b1 * sq_rate - rst_b2g * sq_rate**2

rst_yield_ref =  rst_b0 + rst_b1 * sq_rate + rst_b2g * sq_rate**2
rst_yield_obs =  rst_b0 + rst_b1 * inp_rate + rst_b2g * inp_rate**2
rst_yield_opt =  rst_b0 + rst_b1 * rst_optr + rst_b2g * rst_optr**2


rst_net_ref = rst_yield_ref * cotton_price - sq_rate * nitrogen_cost
rst_net_obs = rst_yield_obs * cotton_price - inp_rate * nitrogen_cost 
rst_net_opt = rst_yield_opt * cotton_price - rst_optr * nitrogen_cost 


hist(rst_optr)
hist(rst_yield_ref)
plot(rst_yield_ref)

max_net_profit = (rst_net_opt - rst_net_ref)
rst_trial_loss = (rst_net_obs - rst_net_ref)

hist(max_net_profit[])
hist(rst_trial_loss[])

mean(max_net_profit[])
mean(rst_trial_loss[])


plot(rst_net_opt - rst_net_ref)
plot(rst_net_opt - rst_net_obs)
plot(rst_net_obs - rst_net_ref)

  




grdf = stack(inp_rate, rst_yield, trial$size15_rep)
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

trial_pred = trial_pred %>% mutate(optr = -0.5 * (b1 - nitrogen_ratio) / b2)
hist(trial_pred$optr)

trial_block = trial_pred %>% group_by(block) %>% summarise_all(mean)

brks = seq(5000,15000,1000)
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
# Evaluate the betas or the optimal rate? Second, plus economical aspects.

plot(rst_optrp, rst_optr, asp = 1)
cor(rst_optrp[], rst_optr[])
rmse(rst_optrp[], rst_optr[])


# rst_opt_resd = (rst_optrp - rst_optr)
# 
# rst_yield_opt =  rst_b0 + rst_b1 * rst_optr + rst_b2g * rst_optr**2
# rst_yield_ref =  rst_b0 + rst_b1 * 100 + rst_b2g * 100**2
# 
# rst_yield_cost = rst_yield_ref - rst_yield
# 
# rst_yield_pred = rst_b0 + rst_b1 * rst_optrp + rst_b2g * rst_optrp**2
# 
# 
# rst_yield_gain = rst_yield_pred - rst_yield_ref
# plot(rst_yield_gain)


rst_yield_ref =  rst_b0 + rst_b1 * 100 + rst_b2g * 100**2
rst_yield_pred = rst_b0 + rst_b1 * rst_optrp + rst_b2g * rst_optrp**2
rst_yield_opt =  rst_b0 + rst_b1 * rst_optr + rst_b2g * rst_optr**2

rst_net_ref = rst_yield_ref * cotton_price - 100 * nitrogen_cost
rst_net_pred = rst_yield_pred * cotton_price - rst_optrp * nitrogen_cost 
rst_net_opt = rst_yield_opt * cotton_price - rst_optr * nitrogen_cost 

rst_net_pred_dif = rst_net_pred - rst_net_ref
rst_net_opt_dif = rst_net_opt - rst_net_ref

mean(rst_net_pred_dif[])
mean(rst_net_opt_dif[])


rst_net_rev <- rst_yield_gain * cotton_price

trial$Yield_tdif <- trial$Yield_pred - trial$Yield_ref
trial$Net_dif <- trial$Yield_tdif * cotton_price

trial <- mutate(trial, Yield_gain_NR = Yield_gain_NR - nitrogen_ratio * NR_opt)




#TODO: add economics:
# How many $/ha less we got using this trial compared to the perfect optimum: 
rst_yield_diff = rst_yield_opt - rst_yield_pred
mean(rst_yield_diff[])

# How many $/ha more we got using this trial compared to the status quo: 
rst_yield_diffp = rst_yield_pred - rst_yield_ref
mean(rst_yield_diffp[])

# How many $/ha more we lost using this trial compared to the status quo: 
mean(rst_yield_cost[])






# create the initial x variable
x1 <- rnorm(100, 15, 5)

# x2, x3, and x4 in a matrix, these will be modified to meet the criteria
x234 <- scale(matrix( rnorm(300), ncol=3 ))

# put all into 1 matrix for simplicity
x1234 <- cbind(scale(x1),x234)

# find the current correlation matrix
c1 <- var(x1234)

# cholesky decomposition to get independence
chol1 <- solve(chol(c1))

newx <-  x1234 %*% chol1 

# check that we have independence and x1 unchanged
zapsmall(cor(newx))
all.equal( x1234[,1], newx[,1] )

# create new correlation structure (zeros can be replaced with other r vals)
newc <- matrix( 
  c(1  , 0.9, 0.5, 0.1, 
    0.9, 1  , 0.5  , 0  ,
    0.5, 0.5  , 1  , 0 ,
    0.1, 0  , 0  , 1  ), ncol=4 )

# check that it is positive definite
eigen(newc)

chol2 <- chol(newc)

finalx <- newx %*% chol2 * sd(x1) + mean(x1)

# verify success
mean(x1)
colMeans(finalx)

sd(x1)
apply(finalx, 2, sd)

zapsmall(cor(finalx))
pairs(finalx)

all.equal(x1, finalx[,1])


# Instead of factorial consider individual scenarios.

# Generate the errors as spatial simulations.
# Aply the errors to the coeffecients after GWR, to speed up process.
# 6 types * 5 replicates * 3 rates * 3 optimum * 10 measurements = 2700


# At least 10 replications of the spatial simulation,
# without nugget and the same range for both betas, using a Gaussian variogram.


# How can we consider the quadratic term in this scenario?
# Using always MGWR could work

# Ranges 0, 50, 100, 150, 200 % of planned rate.
# Ranges 50, 75, 100,  125, 150 % of planned rate.
# Ranges 75, 90, 100,  110, 125 % of planned rate.

# Simulate the status quo rate to be at 75, 100, and 125 % of the true optimum uniform rate
# Keep all trial rates relative to the status quo rate

# Estimate yield loss based on true parameters
# Estimate net revenue based on estimated optimal rates

#Use NDVI as an example, but make it clear that any other tool could be used,
# focus on the correlation between yield response and aux variable response.
# Do not use the correlation with yield itself, rather with the response.
# Create a proxy response function with true slope + error to estimate the measuared variable.
# Consider 0.1, 0.5 and 0.9 correlation

# Remember there is a temporal aspect, by time you make the decision you dont know the future,
# the low correlation shoud also be considered to capture that.
# Y = (true_B0 * alpha + rnd * (1 - alpha) ) + (true_B1 * alpha + rnd * (1 - alpha) ) * NR
# Y = true_B0 + true_B1 * NR + 0.05 * rnd
# Ten different scenarios of measured values.

#Use GWR with adaptive bandwith to estimate optimal application maps
#Use only the points in the trial design (Using all could help with intercept
# estimates but would make worse the slope estimates)
#Consider the yield and the cost if those maps were used to get the net revenue.


