## ------------------------------------------------------------------------
library(sf)
library(raster)
library(gstat)
library(spmoran)
library(tibble)
library(dplyr)
library(tidyr)

## ------------------------------------------------------------------------

set.seed(12345)

# Fixed set of parameters to control the response function:
cotton_price = 1.5 #$/kg fiber
nitrogen_cost = 0.9
nitrogen_ratio = nitrogen_cost/cotton_price
sq_rate = 100
eonr = 50
sq_yield = 1500
yield_sd = 100
eonr_sd = 25
rst_b2g = -0.02
test_rates = sq_rate + 50 * c(-1,-0.5,0,0.5,1)

# This creates a simple RCBD trial with 15X15 plots:
f_coords <- cbind(c(0, 900), c(0, 600))
pol <- st_as_sfc(st_bbox(st_multipoint(f_coords)))
field <- st_sf(pol, crs = 32616)
rst <- raster(field, res = 3)
trial <- aggregate(rst, c(5, 5))
block <- aggregate(trial, c(5, 1))
block <- disaggregate(setValues(block, 1:ncell(block)), c(5,1))
trial$n_rate <- tibble(block = block[]) %>%
  group_by(block) %>%
  mutate(n_rate = sample(test_rates, 5)) %>% 
  pull(n_rate)

n_rate <- disaggregate(trial, c(5, 5))
plot(n_rate)

# This creates the spatial distribution of two random variables:
trial_grd <- rasterToPoints(rst, spatial = TRUE)
gridded(trial_grd) <- TRUE
m <- vgm(psill = 1, model = "Gau", range = 100, nugget = 0.1)
g.dummy <- gstat(formula = z ~ 1, dummy = TRUE, beta = 0, model = m, nmax = 10)
rst_sim <- stack(predict(g.dummy, trial_grd, nsim = 2))

# This adjust the simulated values to the desired scales:
sq_yield_rst <- sq_yield + yield_sd * rst_sim$sim1
eonr_rst <- eonr + eonr_sd * rst_sim$sim2
eonr_rst <- clamp(eonr_rst, 0, 200)


# This is how I use the yield response function to get beta0 and beta1,
#  given the yield at the status quo rate and the EONR:
rst_b1 = -2 * rst_b2g * eonr_rst + nitrogen_ratio
rst_b0 = sq_yield_rst - (rst_b1 * sq_rate + rst_b2g * sq_rate **2) 
yield_obs = rst_b0 + rst_b1 * n_rate + rst_b2g * n_rate **2

# This is from the spmoran function and is
#  a spatial decomposition of the coordinates into eigenvectors:
# meig <- meigen_f(coords = coordinates(rst), model = 'gau', enum = 100)


# This fits the model, with yield as a function of
# spatially varying coefficients for the intercept and N rates
# plus a fixed coefficient for N squared:
y <- yield_obs[]
x <- n_rate[]
xconst <- x**2
res <- resf_vc(y = y, x = x, xconst = xconst, meig = meig, x_sel = FALSE)

# Gets the fitted coeffcients from the model:
sdf <- data.frame(cbind(res$b_vc$`(Intercept)`, res$b_vc$V1, res$c$Estimate))
names(sdf) <- c("b0_pred", "b1_pred", "b2_pred")

# Fill the raster with the fitted values
b0_pred = setValues(rst, sdf$b0_pred) 
b1_pred = setValues(rst, sdf$b1_pred) 
plot(rst_b1)
plot(b1_pred)

# This is the economical analysis:
as_tibble(stack(rst_b0, rst_b1, n_rate)[]) %>%
  setNames(c("b0_obs", "b1_obs", "n_rate")) %>% 
  bind_cols(sdf) %>% 
  mutate(
    b2_obs = rst_b2g,
    yield_obs = b0_obs + b1_obs * n_rate + b2_obs * n_rate**2,
    yield_ref = b0_obs + b1_obs * sq_rate + b2_obs * sq_rate**2,
    yield_opt_f = b0_obs + b1_obs * eonr + b2_obs * eonr**2,
    eonr_obs = (b1_obs - nitrogen_ratio) / (-2 * b2_obs),
    eonr_pred = clamp((b1_pred - nitrogen_ratio) / (-2 * b2_pred),0, 200),
    yield_opt_max = b0_obs + b1_obs * eonr_obs + b2_obs * eonr_obs**2,
    yield_opt_pred = b0_obs + b1_obs * eonr_pred + b2_obs * eonr_pred**2,
    
    net_obs = yield_obs * cotton_price - n_rate * nitrogen_cost,
    net_ref = yield_ref * cotton_price - sq_rate * nitrogen_cost,
    net_opt_f = yield_opt_f * cotton_price - eonr * nitrogen_cost,
    net_opt_max = yield_opt_max * cotton_price - eonr_obs * nitrogen_cost,
    net_opt_pred = yield_opt_pred * cotton_price - eonr_pred * nitrogen_cost
  ) %>%
  summarise(
    net_loss = mean(net_obs - net_ref),
    net_fixed = mean(net_opt_f - net_ref),
    net_max = mean(net_opt_max - net_ref),
    net_diff = mean(net_opt_pred - net_ref),
    net_pct = net_diff / mean(net_max),
    avg_nrate = mean(eonr_pred)
  )


