## ------------------------------------------------------------------------
library(readr)
library(tibble)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

dfl = list()
csv_files = c('trial_eco_lin', 'trial_eco_linf', 'trial_eco_quad', 'trial_eco')
type_names = c('Global Quadratic', 'Global Linear', 'Spatial Quadratic', 'Spatial Linear')

for(csv_file in csv_files){
  df = read_csv(paste0('data/', csv_file, '.csv'))
  df$Type = csv_file
  dfl[[csv_file]] = df
}


df = dfl %>% bind_rows() %>% 
  separate(Trial, c('Size', 'Pct'), sep = '_') %>%
  group_by(opt_pct) %>% 
  mutate(Type = factor(Type, levels = csv_files, labels = type_names),
         Size = gsub('size01', 'size12', Size),
         Size = as.factor(as.integer(gsub('\\D', '', Size))),
         Pct = as.factor(as.integer(gsub('\\D', '', Pct))),
         EONR = paste0('Average EONR: ', opt_pct, ' kg/ha'),
         EONR = factor(EONR, levels = unique(EONR)),
         Range = paste0('Range: ', rate_rng, '%'),
         Range = factor(Range, levels = unique(Range)),
         Maxf = max(VRT)
  )



# ggplot(df, aes(x = Pct, y =  AVG, col = Size)) +
#   geom_boxplot() + theme_classic()+ 
#   labs(y = 'Revenue loss ($/ha)',
#        x = 'Trial coverage (%)',
#        col = 'Size (m)')

# ggplot(df, aes(x = Pct, y =  VRT, col = Size)) +
#   geom_boxplot() + theme_classic() + ylim(0,25)+ 
#   labs(y = 'Revenue diffrence ($/ha)',
#        x = 'Trial coverage (%)',
#        col = 'Size (m)')

# df %>% 
#   ggplot(aes(x = Pct, y = VRT, col = Size)) +
#   geom_boxplot() + #theme_classic() +
#   labs(
#     y = "Revenue difference ($/ha)",
#     x = "Trial coverage (%)",
#     col = "Size (m)"
#   ) +
#   facet_grid(Range ~ EONR)
# 
# ggsave(paste0('figures/Revenue', csv_file, '.png'), width = 9, height = 6)


# df %>% mutate(VRT = pmax(VRT, -3)) %>%
#   filter(opt_pct == 100) %>% 
#   ggplot(aes(x = Pct, y = VRT, col = Size)) +
#   geom_violin(scale = "width") + #theme_classic() +
#   labs(
#     y = "Revenue difference ($/ha)",
#     x = "Trial coverage (%)",
#     col = "Size (m)"
#   ) +
#   facet_grid(Range ~ Type)

# ggsave(paste0('figures/Revenue', csv_file, '_pos.png'), width = 9, height = 6)

# Calc the probability of being wihtin $1 of the optimal rate
df %>% filter(opt_pct == 100) %>% 
  group_by(Type, Range, Size, Pct) %>% 
  summarise(Prob = 100 * mean(abs(AVG - opt_pct) < 5)) %>% 
  ggplot(aes(x = Pct, y = Prob, col = Size, group=Size)) +
  geom_line() +
  labs(
    y = "Success rate (%)",
    x = "Trial coverage (%)",
    col = "Size (m)"
  ) +
  facet_grid(Range ~ Type)
ggsave(paste0('figures/AVG_Success_rate.png'), width = 9, height = 6)


df %>% filter(Type == 'Spatial Linear') %>% 
  group_by(EONR, Range, Size, Pct) %>% 
  summarise(Prob = 100 * mean(abs(AVG - opt_pct) < 5)) %>% 
  ggplot(aes(x = Pct, y = Prob, col = Size, group=Size)) +
  geom_line() +
  labs(
    y = "Success rate (%)",
    x = "Trial coverage (%)",
    col = "Size (m)"
  ) +
  facet_grid(Range ~ EONR)
ggsave(paste0('figures/AVG_Prob_EONR.png'), width = 9, height = 6)


df %>% filter(Type == 'Spatial Linear') %>% 
  group_by(EONR, Range, Size, Pct) %>% 
  summarise(Prob = median(VRT)) %>% 
  ggplot(aes(x = Pct, y = Prob, col = Size, group=Size)) +
  geom_line() +
  labs(
    y = "Success rate (%)",
    x = "Trial coverage (%)",
    col = "Size (m)"
  ) +
  facet_grid(Range ~ EONR)
ggsave(paste0('figures/AVG_Spatial_Return.png'), width = 9, height = 6)


df %>% 
  ggplot(aes(x = Pct, y = AVG, col = Size)) +
  geom_boxplot() + #theme_classic() +
  labs(
    y = "Predicted EONR (kg/ha)",
    x = "Trial coverage (%)",
    col = "Size (m)"
  ) +
  facet_grid(Range ~ EONR)
ggsave(paste0('figures/AVG_', csv_file, '.png'), width = 9, height = 6)


df %>% 
  ggplot(aes(x = Pct, y = Loss, col = Size)) +
  geom_boxplot() + #theme_classic() +
  labs(
    y = "Revenue difference ($/ha)",
    x = "Trial coverage (%)",
    col = "Size (m)"
  ) +
  facet_grid(~Range)

ggsave('figures/Trial_loss.png', width = 9, height = 6)



# df %>% filter(opt_pct == 100, Pct == 100) %>% 
#   arrange(VRT)


# ggplot(df, aes(x = Pct, y =  VRTpct, col = Size)) +
#   geom_boxplot() + theme_classic() + ylim(0,100)+ 
#   labs(y = 'Optimization (%)',
#        x = 'Trial coverage (%)',
#        col = 'Size (m)')
# ggsave('figures/Trial_results_pct.png', width = 9, height = 6)


# ggplot(df, aes(x = Pct, y =  Loss , col = Size)) +
#   geom_boxplot() + theme_classic() +
#   labs(y = 'Revenue loss ($/ha)',
#        x = 'Trial coverage (%)',
#        col = 'Size (m)')
# 
# ggsave('figures/Trial_loss.png', width = 9, height = 6)


