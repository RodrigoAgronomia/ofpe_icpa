## ------------------------------------------------------------------------
library(sf)
library(raster)
library(tibble)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

cotton_price = 1.5
nitrogen_cost = 0.9
nitrogen_ratio = nitrogen_cost/cotton_price
sq_rate = 100
sq_yield = 1500

opt_pct = c(50, 75, 100)
opt_rate = opt_pct/100 * sq_rate

b2 = -0.02
b1 = -2 * b2 * opt_rate + nitrogen_ratio
b0 = sq_yield - (b1 * opt_rate + b2 * opt_rate **2)

# * (1 + 0.25 * b1) 
# rst_optr = np.clip(rst_optr, 0, 200)
# rst_b0 = sq_yield * (1 + 0.05 * b0_true)
# rst_b0 = rst_b0 - (rst_b1 * rst_optr + rst_b2g * rst_optr **2) 
df = crossing(N = seq(0,200),
              b0 = c(4000, 4500, 5000)/3,
              b1 = b1) %>% 
  group_by(b0, b1) %>% 
  mutate(Trial = cur_group_id()) %>% ungroup() %>% 
  mutate(Trial = as.factor(Trial))

df0 = df %>% 
  mutate(b0 = b0 - (mean(b1) * sq_rate + b2 * sq_rate **2),
         Yield = b0 + b1 * N + b2 * N**2,
         NetYield = Yield - nitrogen_ratio * N
  ) %>% group_by(Trial)


df1 = df %>% 
  mutate(b0 = b0 - (b1 * sq_rate + b2 * sq_rate **2),
         Yield = b0 + b1 * N + b2 * N**2,
         NetYield = Yield - nitrogen_ratio * N
  ) %>% group_by(Trial)


df2 = df %>% 
  mutate(eonr = (b1 - nitrogen_ratio) / (-2 * b2),
         b0 = b0 - (b1 * eonr + b2 * eonr **2),
         Yield = b0 + b1 * N + b2 * N**2,
         NetYield = Yield - nitrogen_ratio * N
  ) %>% group_by(Trial)


ggplot(df0, aes(x = N, y = Yield, group = Trial, col = Trial)) +
  geom_line() + theme_classic() + ylim(1000, 2000) +
  geom_point(data = df0 %>% slice(which.max(NetYield)), size = 3) +
  labs(x = 'Nitrogen rate (kg/ha)', y = 'Yield (kg/ha)')
ggsave('figures/Yield_response_zerocenter.png', width = 9, height = 6)

ggplot(df1, aes(x = N, y = Yield, group = Trial, col = Trial)) +
  geom_line() + theme_classic() + ylim(1000, 2000) +
  geom_point(data = df1 %>% slice(which.max(NetYield)), size = 3) +
  labs(x = 'Nitrogen rate (kg/ha)', y = 'Yield (kg/ha)')
ggsave('figures/Yield_response_sqcenter.png', width = 9, height = 6)

ggplot(df2, aes(x = N, y = Yield, group = Trial, col = Trial)) +
  geom_line() + theme_classic() + ylim(1000, 2000) +
  geom_point(data = df2 %>% slice(which.max(NetYield)), size = 3) +
  labs(x = 'Nitrogen rate (kg/ha)', y = 'Yield (kg/ha)')
ggsave('figures/Yield_response_eorncenter.png', width = 9, height = 6)



rst_sim = readRDS('data/GSTAT_sim.rds')

b0_true = rst_sim$sim2
b1_true = rst_sim$sim102

sq_yield = 1500
yield_sd = 100
eonr_sd = 25
rst_b2g = -0.02
sq_rate = 100

opt_pcts = c(50, 75, 100)
opt_pct = opt_pcts[3]
rst_optr = opt_pct + eonr_sd * b1_true
rst_optr = clamp(rst_optr, 0, 200)

plot(rst_optr)

dfh = crossing(opt_pct = opt_pcts, b1_true = rnorm(10000)) %>% 
  mutate(eonr = clamp(sq_rate * (opt_pct/100 + (0.25 * b1_true)),0, 200))

dft = crossing(treat = seq(-2,2), rate_rng = 2*c(12.5, 25, 50)) %>% 
  mutate(rate = sq_rate * (1 + rate_rng/100 * 0.5*treat),
         y = 0.015 + 0.001 * as.integer(as.factor(rate_rng)))

ggplot(dfh, aes(x = eonr, fill = as.factor(opt_pct))) +
  geom_density(alpha = 0.2) +
  geom_point(data = dft, aes(x = rate, y = y, col = as.factor(rate_rng)))+
  labs(x = 'EONR (kg/ha)', y = 'Density', fill = 'Opt (%)', col = 'Range (%)')
ggsave('figures/EONR_hist.png', width = 9, height = 6)

q = c(0.10,0.25,0.5,0.75,0.90)
dfr = dfh %>% group_by(opt_pct) %>% 
  summarise(eonr = quantile(eonr, q), q = 100 * q) %>% 
  crossing(N = seq(0,200),
           b0 = c(4000, 4500, 5000)
  ) %>% 
  group_by(opt_pct, q, b0) %>% 
  mutate(Trial = cur_group_id()) %>% ungroup() %>% 
  mutate(Trial = as.factor(Trial)) %>% 
  mutate(
    b1 = -2 * b2 * eonr + nitrogen_ratio,
    b0n = b0 - (b1 * eonr + b2 * eonr **2),
    Yield = b0n + b1 * N + b2 * N**2,
    NetYield = Yield - nitrogen_ratio * N,
    EONR = paste0('Average EONR: ', opt_pct, ' kg/ha'),
    EONR = factor(EONR, levels = unique(EONR))
  ) %>% group_by(Trial)

ggplot(dfr, aes(x = N, y = Yield, group = Trial, col = q)) +
  geom_line() + theme_classic() + ylim(3300, 5100) +
  geom_point(data = dfr %>%
               slice(which.max(NetYield)),
             size = 3) +
  labs(x = 'Nitrogen rate (kg/ha)',
       y = 'Yield (kg/ha)',
       col = 'EONR quantile (%)'
       ) +
  theme(plot.subtitle = element_text(hjust = 0.5),
        legend.position = c(0.5, 0.15),
        legend.direction = 'horizontal'
        ) +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5)) +
  facet_wrap(~EONR)

ggsave('figures/Yield_response_eonr.png', width = 9, height = 6)



dfr %>% filter(b0 == 4500) %>% 
  inner_join(dft, by = c('N' = 'rate')) %>% 
  mutate(Range = paste0('Range: ', rate_rng, '%'),
         Range = factor(Range, levels = unique(Range))
  ) %>% 
  ggplot(aes(x = N, y = Yield, group = Trial, col = q)) +
  geom_line() + theme_classic() + 
  geom_point(data = dfr %>%
               filter(b0 == 4500) %>%
               slice(which.max(NetYield)),
             size = 3) +
  labs(x = 'Nitrogen rate (kg/ha)',
       y = 'Yield (kg/ha)',
       col = 'EONR quantile (%)'
  ) +
  theme(plot.subtitle = element_text(hjust = 0.5),
        legend.position = c(0.5, 0.1),
        legend.direction = 'horizontal'
  ) +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5)) +
  xlim(0, 200) + ylim(3800, 4800) +
  facet_grid(Range ~ EONR)

ggsave('figures/Yield_response_eonr_range.png', width = 9, height = 6)

