## ------------------------------------------------------------------------
library(spmoran)
library(dplyr)
## ------------------------------------------------------------------------
trial_rst <- readRDS("data/Trial_Design.rds")
rst = raster(trial_rst)


plot(b0_true)

rstt = stack('figures/b0_true.png')
rreps = raster(nrow = 10, ncol = 10)
extent(rreps) = extent(rstt)
pols = rasterToPolygons(rreps)

pdf('figures/b0_true_col.pdf', width = 8.2, height = 6)
plot(rstt)
plot(pols, add = TRUE)
dev.off()


# meig <- readRDS("data/data_meigen.rds")
meig <- meigen_f(coords = coordinates(rst), model = 'gau', enum = 100)


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


arg_list = expand.grid(rate_rng = rate_rngs, opt_pct = opt_pcts, tname = tnames) %>% 
  mutate(df_names = paste0("data/data_trial_", tname,
                           "_rng_", rate_rng,
                           "_opt_",opt_pct, "_2.rds"),
         res_names = paste0("data/res_trial_", tname,
                            "_rng_", rate_rng,
                            "_opt_", opt_pct, "_2.rds")
  )

# lapply(arg_list$res_names, unlink)
arg_list$Done = unlist(lapply(arg_list$res_names, file.exists))
# arg_list = filter(arg_list, !Done)
# arg_list %>% distinct(res_names)

tidx = 125
df_name = arg_list$df_names[tidx]
res_name = arg_list$res_names[tidx]
  
df <- readRDS(df_name)

sq_rate = 100
df$n_rate = df$n_rate - sq_rate

yield_obs <- df$b0_obs + df$b1_obs * df$n_rate + df$b2_obs * df$n_rate**2
yield_obs_rst = setValues(rst, yield_obs)

plot(yield_obs_rst)

x <- df$n_rate
xconst <- df$n_rate**2
res <- resf_vc(y = yield_obs, x = x,
               # xconst = xconst, 
               meig = meig, x_sel = FALSE)
  
sdf <- data.frame(cbind(res$b_vc$`(Intercept)`, res$b_vc$V1))
names(sdf) <- c("b0_pred", "b1_pred")
sdf <- bind_cols(df, sdf)

# saveRDS(sdf, res_name)

n_rate = setValues(rst, sdf$n_rate) 
b0_obs = setValues(rst, sdf$b0_obs) 
b1_obs = setValues(rst, sdf$b1_obs) 
b0_pred = setValues(rst, sdf$b0_pred) 
b1_pred = setValues(rst, sdf$b1_pred) 

plot(n_rate)
plot(b0_obs)
plot(b0_pred)
plot(b1_obs)
plot(b1_pred)

plot(b0_obs ~ b0_pred, sdf, pch = '.', asp = 1)
abline(a = 0, b = 1, col = 'red')
plot(b1_obs ~ b1_pred, sdf, pch = '.', asp = 1)
abline(a = 0, b = 1, col = 'red')


