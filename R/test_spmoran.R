## ------------------------------------------------------------------------
library(sf)
library(raster)
library(gstat)
library(tibble)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(spmoran)


## ------------------------------------------------------------------------

set.seed(12345)

trial_rst <- readRDS('data/Trial_Design.rds')
rst_sim = readRDS('data/GSTAT_sim.rds')

plot(trial_rst$size03_p010)

bw = 200
cotton_price = 0.4
nitrogen_cost = 1.0
nitrogen_ratio = nitrogen_cost/cotton_price
sq_rate = 100
sq_yield = 4500
rst_b2g = -0.02

b0_true = rst_sim$sim1
b1_true = rst_sim$sim100


opt_pct = 150
rst_optr = opt_pct/100 * sq_rate * (1 + 0.25 * b1_true) 
rst_optr = clamp(rst_optr, 0, 200)

rst_b1 = -2 * rst_b2g * rst_optr + nitrogen_ratio
rst_b2 = rst_b2g * (1 + (0 * rst_b1))
rst_b0 = sq_yield * (1 + 0.05 * b0_true)
rst_b0 = rst_b0 - (rst_b1 * rst_optr + rst_b2g * rst_optr **2) 

plot(rst_b0)
plot(rst_b1)

rep_rst = trial_rst$size60_rep
inp_rate = sq_rate + 50 * trial_rst$size60_p050
inp_rate[is.na(inp_rate)] = sq_rate
rst_yield_obs =  rst_b0 + rst_b1 * inp_rate + rst_b2 * inp_rate**2
plot(rst_yield_obs)

plot(trial_rst$size60_rep %% 5)


# lm_rst = stack(rst_yield_obs, inp_rate, rep_rst) %>%
#   rasterToPoints() %>% as_tibble() %>% 
#   setNames(c('x', 'y','yield', 'n_rate', 'block')) %>% 
#   filter(!is.na(n_rate)) %>% 
#   mutate(block = as.factor(block))
#   
# table(lm_rst$block)
# 
# ggplot(lm_rst, aes(x = n_rate, y = yield, group = n_rate)) + 
#   geom_boxplot() + facet_wrap(~block)
# 
# lm0 = lm(yield ~ poly(n_rate, 2) + poly(n_rate, 1) * block, lm_rst)

lm_rst = stack(rst_yield_obs, inp_rate) %>%
  rasterToPoints() %>% as_tibble() %>% 
  setNames(c('x', 'y','yield', 'n_rate')) %>% 
  filter(!is.na(n_rate))

rst_obs = stack(rst_b0, rst_b1, rst_b2, rst_optr) %>%
  rasterToPoints() %>% as_tibble() %>% 
  setNames(c('x', 'y','b0_obs', 'b1_obs', 'b2_obs', 'eonr_obs')) %>%
  mutate(eonr_pred = (b1_obs - nitrogen_ratio)/(-2 * b2_obs))


n_rate = poly(lm_rst$n_rate, 2)
lm_rst$n_rate1 = n_rate[,1]
lm_rst$n_rate2 = n_rate[,2]

lm_rst$n_rate1 = lm_rst$n_rate
lm_rst$n_rate2 = lm_rst$n_rate ** 2

y	      <- lm_rst$yield
x       <- lm_rst$n_rate1
xconst  <- lm_rst$n_rate2

coords  <- lm_rst[,c("y","x")]
meig 	  <- meigen_f(coords=coords, enum = 500)

res	    <- resf_vc(y=y, x=x, xconst=xconst, meig=meig, x_nvc = FALSE)

lm_rst$beta0 = res$b_vc$`(Intercept)`
lm_rst$beta1 = res$b_vc$V1
lm_rst$beta2 = res$c$Estimate
lm_pred = lm_rst %>%
  mutate(yield_pred = beta0 + beta1 * n_rate1 + beta2 * n_rate2,
         optr_pred = (beta1 - nitrogen_ratio)/(-2 * beta2)
         )

hist(lm_pred$optr_pred)
lm_pred$eonr_obs = rst_obs$eonr_obs
ggplot(lm_pred, aes(x = optr_pred, y = eonr_obs)) + 
         geom_point(pch = '.') + geom_abline()

rst_yield_pred = inp_rate
rst_b0_pred = inp_rate
rst_b1_pred = inp_rate

rst_yield_pred[!is.na(inp_rate)] = lm_pred$yield_pred
rst_b0_pred[!is.na(inp_rate)] = lm_pred$beta0
rst_b1_pred[!is.na(inp_rate)] = lm_pred$beta1

plot(rst_b0)
plot(rst_b0_pred)
plot(rst_b1)
plot(rst_b1_pred)

plot(rst_yield_obs)
plot(rst_yield_pred)

plot(yield ~ yield_pred, lm_pred, pch = '.')

