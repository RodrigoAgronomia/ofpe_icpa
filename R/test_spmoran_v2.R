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

b0_true = rst_sim$sim1
b1_true = rst_sim$sim100

# coords  <- coordinates(trial_rst)
# meig 	  <- meigen_f(coords=coords, enum = 500)

cotton_price = 0.4
nitrogen_cost = 1.0
nitrogen_ratio = nitrogen_cost/cotton_price
sq_rate = 100
sq_yield = 4500
rst_b2g = -0.02

rate_rngs = c(12.5, 25, 50)
opt_pcts = c(50, 75, 100)
tnames =  c("size03_p001", "size03_p005", "size03_p010", "size03_p050", "size03_p100", 
            "size06_p001", "size06_p005", "size06_p010", "size06_p050", "size06_p100", 
            "size15_p001", "size15_p005", "size15_p010", "size15_p050", "size15_p100", 
            "size30_p001", "size30_p005", "size30_p010", "size30_p050", "size30_p100", 
            "size60_p001", "size60_p005", "size60_p010", "size60_p050", "size60_p100", 
            "size01_p001", "size01_p005", "size01_p010", "size01_p050", "size01_p100"
)

for(rate_rng in rate_rngs){
  for(opt_pct in opt_pcts){
    
    rst_optr = opt_pct/100 * sq_rate * (1 + 0.25 * b1_true) 
    rst_optr = clamp(rst_optr, 0, 200)
    
    rst_b1 = -2 * rst_b2g * rst_optr + nitrogen_ratio
    rst_b2 = rst_b2g * (1 + (0 * rst_b1))
    rst_b0 = sq_yield * (1 + 0.05 * b0_true)
    rst_b0 = rst_b0 - (rst_b1 * rst_optr + rst_b2g * rst_optr **2) 
    
    for(tname in tnames){
      
      inp_rate = sq_rate + rate_rng * trial_rst[[tname]]
      inp_rate[is.na(inp_rate)] = sq_rate
      
      df = stack(rst_b0, rst_b1, rst_b2, inp_rate)[] %>%
        as_tibble() %>% setNames(c('b0_obs', 'b1_obs', 'b2_obs', 'n_rate'))
      
      s_name = paste0('data/data_trial_',tname, '_rng_', rate_rng, '_opt_', opt_pct, '_1.rds')
      saveRDS(df, s_name)
      
    }
  }
}



# y	      <- df$yield
# x       <- df$n_rate
# xconst  <- df$n_rate ** 2
# 
# res	    <- resf_vc(y=y, x=x, xconst=xconst, meig=meig, x_nvc = FALSE)
# 
# saveRDS(res, paste0('data/res_rng_', rate_rng, '_opt_', opt_pct, '_1.rds'))
# 
# 
# df <- df %>%
#   mutate(b0_pred = res$b_vc$`(Intercept)`,
#          b1_pred = res$b_vc$V1,
#          b2_pred = res$c$Estimate,
#          yield_pred = b0_pred + b1_pred * n_rate + b2_pred * n_rate2,
#          eonr_obs = (b1_obs - nitrogen_ratio)/(-2 * b2_obs),
#          eonr_pred = (b1_pred - nitrogen_ratio)/(-2 * b2_pred)
#   )
# 
# saveRDS(df, paste0('data/res_rng_', rate_rng, '_opt_', opt_pct, '_1.rds'))

# glimpse(df)
# 
# 
# ggplot(df, aes(x = eonr_pred, y = eonr_obs)) + 
#   geom_point(pch = '.') + geom_abline()
# 
# 
# ggplot(df, aes(x = yield_pred, y = yield)) + 
#   geom_point(pch = '.') + geom_abline()
# 
# rst_yield_pred = inp_rate
# rst_b0_pred = inp_rate
# rst_b1_pred = inp_rate
# rst_eonr_pred = inp_rate
# 
# rst_yield_pred[!is.na(inp_rate)] = df$yield_pred
# rst_b0_pred[!is.na(inp_rate)] = df$b0_pred
# rst_b1_pred[!is.na(inp_rate)] = df$b1_pred
# rst_eonr_pred[!is.na(inp_rate)] = df$eonr_pred
# 
# plot(rst_b0)
# plot(rst_b0_pred)
# plot(rst_b1)
# plot(rst_b1_pred)
# 
# plot(rst_optr)
# plot(rst_eonr_pred)
# 
# plot(rst_yield_obs)
# plot(rst_yield_pred)



