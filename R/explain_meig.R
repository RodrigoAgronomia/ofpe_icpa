## ------------------------------------------------------------------------
library(sf)
library(raster)
library(gstat)
library(tibble)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)


## ------------------------------------------------------------------------
trial_rst <- readRDS("data/Trial_Design.rds")
rst_sim <- readRDS("data/GSTAT_sim.rds")
meig <- readRDS("data/data_meigen.rds")

mm = matrix(meig$sf)
dim(mm) = c(300, 200, 600)
dim(mm)
rst = raster(trial_rst)
rst[] = t(mm[,,101])
plot(rst)

data_file = 'data/data_trial_size03_p100_rng_50_opt_100_1.rds'
ddf <- readRDS(data_file) %>%
  mutate(yield_obs = b0_obs + b1_obs * n_rate + b2_obs * n_rate**2)

res_file = 'data/res_trial_size03_p100_rng_50_opt_100_1.rds'
sdf <- readRDS(res_file)


df = as_tibble(meig$sf) %>% setNames(paste0('X', 1:ncol(meig$sf)))
fml = as.formula(paste0('y~ n_rate * (', paste(names(df)[1:100], collapse = '+'), ')'))
df$y = ddf$yield_obs
df$n_rate = ddf$n_rate

lm0 = lm(fml, df)
slm = broom::tidy(summary(lm0))

mf <- model.frame(lm0, df)
mt <- attr(x = mf, which = "terms")
X <- as_tibble(model.matrix(object = mt, data = mf))
names(X)[1] = "Intercept"

plot(X$'n_rate:X1', X$n_rate * X$X1, pch = '.')


barplot(abs(slm$statistic[3:102]))
barplot(abs(slm$statistic[103:202]))
barplot(abs(slm$statistic[3:202]))
summary(lm0)



