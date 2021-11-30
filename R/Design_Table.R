## ------------------------------------------------------------------------
library(raster)
library(tidyr)
library(dplyr)

trial <- readRDS('data/Trial_Design.rds')

repsn = list()
trial$size03_rep
for(i in 1:6){
  reps = trial[[1 + (i - 1) * 6]]
  for(j in 1:5){
    
    obs = trial[[j + 1 + (i - 1) * 6]]
    repsn[[names(obs)]] = length(unique(reps[!is.na(obs)]))
  }
}

df = data.frame(repsn)

dft = df %>% gather() %>%
  separate(key, c('Size', 'Coverage'), sep = '_') %>%
  spread(Size, value = value)

# Size, Coverage, Observations
dft
write.csv(dft, './data/repn.csv', row.names = FALSE)
