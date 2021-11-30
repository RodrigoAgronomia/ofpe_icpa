## ------------------------------------------------------------------------
library(readr)
library(tibble)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)


df = read_csv('data/trial_spatial_linear.csv') %>% 
  separate(Trial, c('Size', 'Pct'), sep = '_') %>%
  group_by(opt_pct) %>% 
  mutate(Size = gsub('size01', 'size12', Size),
         Size = as.factor(as.integer(gsub('\\D', '', Size))),
         Pct = as.factor(as.integer(gsub('\\D', '', Pct))),
         EONR = paste0('Average EONR: ', opt_pct, ' kg/ha'),
         EONR = factor(EONR, levels = unique(EONR)),
         Range = paste0('Range: ', rate_rng, '%'),
         Range = factor(Range, levels = unique(Range)),
         Maxf = max(VRT)
  )




df %>% group_by(Range, Pct) %>% 
  summarise(Loss = mean(-Loss)) %>% 
  ggplot(aes(x = Pct, y = Loss, col = Range, group=Range)) +
  geom_line() + scale_y_log10() +
  labs(
    y = "Average profit loss ($/ha)",
    x = "Trial coverage (%)"
  ) 

ggsave(paste0('figures/AVG_profit_loss.png'), width = 9, height = 6)


# df %>% group_by(EONR, Range, Size, Pct) %>% 
#   summarise(Prob = 100 * mean(Prob)) %>% 
#   ggplot(aes(x = Pct, y = Prob, col = Size, group=Size)) +
#   geom_line() +
#   labs(
#     y = "Success rate (%)",
#     x = "Trial coverage (%)",
#     col = "Size (m)"
#   ) +
#   facet_grid(Range ~ EONR)


df %>% group_by(Range, Size, Pct) %>%
  summarise(Profit = mean(pmax(VRT + Loss, -15))) %>%
  ggplot(aes(x = Pct, y = Profit, col = Size, group=Size)) +
  geom_line() +
  labs(
    y = "Average net profit difference ($/ha)",
    x = "Trial coverage (%)",
    col = "Size (m)"
  ) +
  facet_wrap(Range ~ .)
# ggsave(paste0('figures/AVG_Spatial_Net_Profit.png'), width = 9, height = 3)

df %>% group_by(EONR, Range, Size, Pct) %>% 
  summarise(Profit = mean(pmax(Total + Loss, -100))) %>% 
  ggplot(aes(x = Pct, y = Profit, col = Size, group=Size)) +
  geom_line() +
  labs(
    y = "Average net profit difference ($/ha)",
    x = "Trial coverage (%)",
    col = "Size (m)"
  ) +
  facet_grid(Range ~ EONR)

ggsave(paste0('figures/AVG_Spatial_Net_Profit.png'), width = 9, height = 6)



df %>% group_by(EONR, Range, Size, Pct) %>% 
  summarise(Profit = mean(pmax(VRT, 0) + Loss)) %>% 
  ggplot(aes(x = Pct, y = Profit, col = Size, group=Size)) +
  geom_line() +
  labs(
    y = "Average profit difference ($/ha)",
    x = "Trial coverage (%)",
    col = "Size (m)"
  ) +
  facet_grid(Range ~ EONR)

