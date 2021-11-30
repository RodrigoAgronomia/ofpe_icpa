## ------------------------------------------------------------------------
library(sf)
library(raster)
library(tibble)
library(dplyr)
library(mgwrsar)

set.seed(12345)

rst <- stack('data/Trial_Design.tif')
trial <- aggregate(raster(rst), 10)


## ------------------------------------------------------------------------
trial_grd <- rasterToPoints(trial, spatial = TRUE)
gridded(trial_grd) <- TRUE

m <- gstat::vgm(psill = 1, model = "Gau",
                range = 100,
                nugget = 0.5)
g.dummy <- gstat::gstat(
  formula = z ~ 1,
  dummy = TRUE, beta = 0,
  model = m, nmax = 10
)
# rst_sim <- predict(g.dummy, trial_grd, nsim = 200)
# rst_sim <- scale(stack(rst_sim))
# rst_sim <- disaggregate(rst_sim, 10)
# rst_b0 <- rst_sim[[1:100]]
# rst_b1 <- rst_sim[[101:200]]

rst_sim <- predict(g.dummy, trial_grd, nsim = 2)
rst_sim <- scale(stack(rst_sim))
rst_sim <- disaggregate(rst_sim, 10)
rst_b0 <- 5000 + 500 * rst_sim[[1]]
rst_b1g <- 20
rst_b2g <- -0.1
rst_b1 <- rst_b1g + 3 * rst_sim[[2]]
rst_optr = -0.5 * (rst_b1 + rst_b1g) / rst_b2g
rst_err <- 100 * rnorm(ncell(rst_sim))

inp_rate = 100 + 50 * rst[[12]]
# inp_rate[is.na(inp_rate)] = 100
rst_yield =  rst_b0 + rst_b1 * inp_rate + rst_b2g * inp_rate**2 + rst_err

plot(rst_yield)
plot(inp_rate, rst_yield)
plot(inp_rate)


grdf = aggregate(stack(inp_rate, rst_yield), 10)
names(grdf) = c('inp_rate', 'yield')
trial = st_as_sf(rasterToPolygons(grdf))
coords <- st_coordinates(st_centroid(st_geometry(trial)))

mgwr.fml <- as.formula(yield~poly(inp_rate, 2))
bwcv <- 20
fx_var = c("poly(inp_rate, 2)2")

mgwr.model <- MGWRSAR(mgwr.fml, trial, 
                      coord = coords, fixed_vars = fx_var, 
                      kernels = c("gauss_adapt"),
                      H = bwcv, Model = "MGWR", control = list(SE = TRUE)
)


print(mgwr.model)

gwr_r <- as_tibble(mgwr.model$Betav)
names(gwr_r) <- paste0(c('Intercept', 'inp_rate'), "_gwr")
trial_gwr <- bind_cols(trial, gwr_r)


mf <- model.frame(mgwr.fml, trial)
mt <- attr(x = mf, which = "terms")
X <- as_tibble(model.matrix(object = mt, data = mf))
names(X)[1] = "Intercept"

Betav <- as_tibble(mgwr.model$Betav)
Betac = rep(mgwr.model$Betac, each = nrow(mgwr.model$Betav))
dim(Betac) = c(nrow(mgwr.model$Betav), length(mgwr.model$Betac))
Betac = as_tibble(Betac)
names(Betac) = names(mgwr.model$Betac)
Betas = cbind(Betav, Betac)
Betas = Betas[names(X)]


iro = seq(min(trial$inp_rate), max(trial$inp_rate), length.out = 100)
dfp = data.frame(inp_rate = iro)

inp_ratep = poly(trial$inp_rate, degree = 2)
trial$beta1 = (inp_ratep[,1] * Betas$`poly(inp_rate, 2)1`)






