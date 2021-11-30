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

b0_true <- rst_sim$sim2
b1_true <- rst_sim$sim102

coords <- coordinates(trial_rst)
rst <- raster(trial_rst)

cotton_price <- 0.4
nitrogen_cost <- 1.0
nitrogen_ratio <- nitrogen_cost / cotton_price

sq_rate <- 100
sq_yield <- 4500
rst_b2g <- -0.02

rate_rngs <- c(12.5, 25, 50)
opt_pcts <- c(50, 75, 100)
tnames <- c(
  "size03_p001", "size03_p005", "size03_p010", "size03_p050", "size03_p100",
  "size06_p001", "size06_p005", "size06_p010", "size06_p050", "size06_p100",
  "size15_p001", "size15_p005", "size15_p010", "size15_p050", "size15_p100",
  "size30_p001", "size30_p005", "size30_p010", "size30_p050", "size30_p100",
  "size60_p001", "size60_p005", "size60_p010", "size60_p050", "size60_p100",
  "size01_p001", "size01_p005", "size01_p010", "size01_p050", "size01_p100"
)

reps = 2
rate_rngs <- c(12.5, 25, 50)
opt_pcts <- c(50, 75, 100)
tnames <- c(
  "size03_p001", "size03_p005", "size03_p010", "size03_p050", "size03_p100",
  "size06_p001", "size06_p005", "size06_p010", "size06_p050", "size06_p100",
  "size15_p001", "size15_p005", "size15_p010", "size15_p050", "size15_p100",
  "size30_p001", "size30_p005", "size30_p010", "size30_p050", "size30_p100",
  "size60_p001", "size60_p005", "size60_p010", "size60_p050", "size60_p100",
  "size01_p001", "size01_p005", "size01_p010", "size01_p050", "size01_p100"
)


arg_list <- expand.grid(rate_rng = rate_rngs, opt_pct = opt_pcts, tname = tnames, rep = reps) %>%
  mutate(
    df_names = paste0(
      "data/data_trial_", tname,
      "_rng_", rate_rng,
      "_opt_", opt_pct, "_", rep, ".rds"
    ),
    res_names = paste0(
      "data/res_trial_", tname,
      "_rng_", rate_rng,
      "_opt_", opt_pct, "_", rep, ".rds"
    )
  )

arg_list$Done <- unlist(lapply(arg_list$res_names, file.exists))
arg_list <- filter(arg_list, Done)

details <- arg_list %>% 
  separate(tname, into = c("Size", "Pct"), sep = "_") %>% 
  mutate(i = row_number()) %>%
  mutate(
    Size = gsub("size01", "size09", Size),
    Size = as.integer(gsub("\\D", "", Size)),
    Pct = as.integer(gsub("\\D", "", Pct))
  ) %>%
  tibble() %>% 
  select(i, Size, Pct, rate_rng, opt_pct)


i <- 2
df <- bind_rows(lapply(1:nrow(arg_list), function(i) {
  sdf <- bind_cols(readRDS(arg_list$df_name[i]), readRDS(arg_list$res_name[i]))
  sdf$i <- i
  sdf
}))


dfj <- df %>%
  left_join(details, by = "i") %>%
  mutate(
    yield_obs = b0_obs + b1_obs * n_rate + b2_obs * n_rate**2,
    yield_ref = b0_obs + b1_obs * sq_rate + b2_obs * sq_rate**2,
    yield_opt_f = b0_obs + b1_obs * opt_pct + b2_obs * opt_pct**2,
    eonr_obs = (b1_obs - nitrogen_ratio) / (-2 * b2_obs),
    eonr_pred = clamp((b1_pred - nitrogen_ratio) / (-2 * b2_pred),0, 200),
    yield_opt_max = b0_obs + b1_obs * eonr_obs + b2_obs * eonr_obs**2,
    yield_opt_pred = b0_obs + b1_obs * eonr_pred + b2_obs * eonr_pred**2,
    
    net_obs = yield_obs * cotton_price - n_rate * nitrogen_cost,
    net_ref = yield_ref * cotton_price - sq_rate * nitrogen_cost,
    net_opt_f = yield_opt_f * cotton_price - opt_pct * nitrogen_cost,
    net_opt_max = yield_opt_max * cotton_price - eonr_obs * nitrogen_cost,
    net_opt_pred = yield_opt_pred * cotton_price - eonr_pred * nitrogen_cost
  )

glimpse(dfj)


dfm <- dfj %>%
  group_by(i, Size, Pct, rate_rng, opt_pct) %>%
  summarise(
    b2_pred = mean(b2_pred),
    net_loss = mean(net_obs - net_ref),
    net_fixed = mean(net_opt_f - net_ref),
    net_max = mean(net_opt_max - net_ref),
    net_diff = mean(net_opt_pred - net_ref),
    net_pct = net_diff / mean(net_max),
    avg_nrate = mean(eonr_pred)
  ) %>%
  arrange(Size, Pct, rate_rng, opt_pct)

dfm %>%  filter(net_diff > -25) %>% 
  mutate(VRT = net_diff - net_fixed) %>% 
  ggplot(aes(x = as.factor(Pct), y = net_diff, col = as.factor(Size))) +
  geom_point() +
  theme_classic() +
  # ylim(-1, 25) +
  labs(
    y = "Revenue difference ($/ha)",
    x = "Trial coverage (%)",
    col = "Size (m)"
  ) +
  facet_wrap(~ rate_rng + opt_pct, scales = 'free_y')

glimpse(dfm)

dfm %>% filter(net_diff < 0) %>% arrange(net_diff)

# df_rst = df
# coordinates(df_rst) = coords
# gridded(df_rst) = TRUE
# df_rst = stack(df_rst)
#
# ggplot(df, aes(x = eonr_pred, y = eonr_obs)) +
#   geom_point(pch = '.') + geom_abline()
#
# ggplot(df, aes(x = yield_pred, y = yield_obs)) +
#   geom_point(pch = '.') + geom_abline()
#
#
# plot(df_rst$b0_obs)
# plot(df_rst$b0_pred)
# plot(df_rst$b1_obs)
# plot(df_rst$b1_pred)
#
# plot(df_rst$eonr_obs)
# plot(df_rst$eonr_pred)
# plot(df_rst$yield_obs)
# plot(df_rst$yield_pred)


# 
# df = data.frame(matrix(rnorm(100000), ncol = 100))
# fml = as.formula(paste0('y~', paste(names(df), collapse = '+')))
# df$y = rnorm(nrow(df))
# 
# lm0 = lm(fml, df)
# summary(lm0)
# summary.aov(lm0)

  