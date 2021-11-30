## ------------------------------------------------------------------------
library(doParallel)  
library(spmoran)
library(dplyr)
## ------------------------------------------------------------------------

calc_resf = function(meig, df_name, res_name){
  
  df <- readRDS(df_name)
  
  yield_obs <- df$b0_obs + df$b1_obs * df$n_rate + df$b2_obs * df$n_rate**2
  
  y <- yield_obs
  x <- df$n_rate
  xconst <- df$n_rate**2
  
  res <- resf_vc(y = y, x = x, xconst = xconst, meig = meig, x_nvc = TRUE)
  
  sdf <- data.frame(cbind(res$b_vc$`(Intercept)`, res$b_vc$V1, res$c$Estimate))
  names(sdf) <- c("b0_pred", "b1_pred", "b2_pred")
  
  saveRDS(sdf, res_name)
  return(TRUE)
}

no_cores <- 15
registerDoParallel(cores=no_cores)  
cl <- makeCluster(no_cores) 


meig <- readRDS("data/data_meigen.rds")
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
arg_list = filter(arg_list, !Done)
# arg_list %>% distinct(res_names)

results <- foreach(i=1:nrow(arg_list), .packages='spmoran') %dopar% {
  calc_resf(meig, arg_list$df_names[i], arg_list$res_names[i])
}

stopCluster(cl) 

