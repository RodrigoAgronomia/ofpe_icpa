## ------------------------------------------------------------------------
library(sf)
library(tmap)
library(raster)
library(tibble)
library(dplyr)
library(stringr)

set.seed(12345)

n_treats = 5
treat <- seq(1, n_treats)
tidx = rep(1:n_treats, 5) + rep(1:n_treats, each = 5)
tidx = tidx %% n_treats + 1
tidx = matrix(tidx, ncol = n_treats) - 3


f_coords <- cbind(c(0, 900), c(0, 600))
pol <- st_as_sfc(st_bbox(st_multipoint(f_coords)))
field <- st_sf(pol, crs = 32616)
rst <- raster(field, res = 3)


plot_sizes = c(3, 6, 15, 30, 60)

plot_size = 60
for(plot_size in plot_sizes){
  trial <- aggregate(raster(rst), plot_size/3)
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
  
  c01 = ceiling(n_cols/5)
  c01 = ceiling(n_cols / c01 * (seq(1, c01) - 0.5))
  c05 = ceiling(n_cols/3)
  c05 = ceiling(n_cols / c05 * (seq(1, c05) - 0.5))
  c10 = ceiling(n_cols/2)
  c10 = ceiling(n_cols / c10 * (seq(1, c10) - 0.5))
  
  t01 = (trial$rep_col %in% c01) & ((trial$prow + 9) %% 20 == 0)
  t05 = (trial$rep_col %in% c05) & ((trial$prow + 4) %% 7 == 0)
  if(plot_size == 60){
    t01 = (trial$rep_col %in% c01) & ((trial$prow + 5) %% 10 == 0)
    t05 = (trial$rep_col %in% c05) & ((trial$prow + 2) %% 5 == 0)
  }
  t10 = (trial$rep_col %in% c10) & ((trial$prow + 2) %% 5 == 0)
  t50 = xor(trial$rep_col %% 2, trial$prow %% 2)
  t100 = trial$rep_col > 0
  
  srep = disaggregate(trial$rep, plot_size/3)
  names(srep) = paste0('size', str_pad(plot_size, 2, pad = 0),'_rep')
  rst = addLayer(rst, srep)
  
  s_pct = list('p001'=t01, 'p005'=t05, 'p010'=t10, 'p050'=t50, 'p100'=t100)
  for(s in names(s_pct)){
    streat = trial$treat
    streat[!s_pct[[s]]] = NA
    streat = disaggregate(streat, plot_size/3)
    names(streat) = paste0('size', str_pad(plot_size, 2, pad = 0), '_',s)
    rst = addLayer(rst, streat)
  }
}

plot_size = 3
trial <- raster(rst)
trial$id <- 1:ncell(trial)
trial$pcol = 1:ncol(trial)
trial$rep_col = ceiling(trial$pcol/n_treats)
trial$prow = rep(1:nrow(trial), each = ncol(trial))
trial$rep_row = ceiling(trial$prow / 4)
trial$rep = 1e6 * trial$rep_row + trial$rep_col
trial$rep = as.numeric(as.factor(trial$rep[]))
trial$block = 1e6 * trial$rep_row + trial$rep_col
trial$block = as.numeric(as.factor(trial$block[]))



n_rows = max(trial$rep_row[])
n_cols = max(trial$rep_col[])
n_blocks = max(trial$block[])

trial$treat = c(-2, 0, 2, 0)[1 + trial$prow[] %% 4]

t01 = trial$rep_row == 25
t05 = trial$rep_row %in% c(5, 25, 45)
t10 = trial$rep_row %in% seq(5, 45, 10)
t50 = trial$rep_row %in% seq(1, 50, 2)
t100 = trial$rep_col > 0

srep = trial$rep
names(srep) = paste0('size01_rep')
rst = addLayer(rst, srep)

s_pct = list('p001'=t01, 'p005'=t05, 'p010'=t10, 'p050'=t50, 'p100'=t100)
for(s in names(s_pct)){
  streat = trial$treat
  streat[!s_pct[[s]]] = NA
  names(streat) = paste0('size01_',s)
  rst = addLayer(rst, streat)
}


rst_treat = rst[[grep('_p', names(rst))]]
rst_treat = rst_treat[[as.numeric(t(matrix(1:30, 5, 6)))]]

saveRDS(rst, 'data/Trial_Design.rds')

writeRaster(rst, 'data/Trial_Design.tif')

tms = tm_shape(rst_treat) + tm_raster(legend.show = FALSE, style = 'cat') + tm_facets(nrow = 5)
tmap_save(tms, filename = 'figures/Trial_Design.png', dpi = 600, width = 9, height = 4.8, units = 'in')

