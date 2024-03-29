## ------------------------------------------------------------------------
library(sf)
library(raster)
library(tibble)
library(dplyr)
library(gstat)
library(mgwrsar)
library(reticulate)
np <- import("numpy")

set.seed(12345)

trial <- readRDS('data/Trial_Design.rds')


## ------------------------------------------------------------------------
trial_grd <- rasterToPoints(trial, spatial = TRUE)
gridded(trial_grd) <- TRUE


dMat <- gw.dist(trial_grd, trial_grd)

wwf = list()

i = 620
dist.vi <- dMat[, i]

i = 1
for(bw in bws){
  rst$W.i <- gw.weight(dist.vi, bw, kernel = "gaussian", adaptive = TRUE)
  wwf[[paste0('W', bw)]] = round(255 * rst$W.i)
}
wwr = stack(wwf)
wwm = as.array(wwr)

wfile = './data/weigths_zoom.npy'
np$save(wfile, wwm)



m <- vgm(psill = 1, model = "Gau", range = 100, nugget = 0.1)
g.dummy <- gstat(formula = z ~ 1, dummy = TRUE, beta = 0, model = m, nmax = 10)

rst_sim <- predict(g.dummy, trial_grd, nsim = 2)
rst_sim <- scale(stack(rst_sim))
rst_b0 <- 10000 + 1000 * rst_sim[[1]]
rst_b1g <- 20
rst_b2g <- -0.1
rst_b1 <- rst_b1g + 3 * rst_sim[[2]]
rst_optr = -0.5 * (rst_b1 + rst_b1g) / rst_b2g
rst_err <- 100 * rnorm(ncell(rst_sim))

inp_rate = 100 + 50 * rst[[12]]
# inp_rate[is.na(inp_rate)] = 100
rst_yield =  rst_b0 + rst_b1 * inp_rate + rst_b2g * inp_rate**2

plot(rst_yield)
plot(inp_rate, rst_yield)
plot(inp_rate)

# Fit a global model first to speed up.

grdf = stack(inp_rate, rst_yield)
names(grdf) = c('inp_rate', 'yield')
ptsdf <- rasterToPoints(grdf, spatial = TRUE)@data
coords <- coordinates(ptsdf)

mgwr.fml <- as.formula(yield~poly(inp_rate, 2))
bwcv <- 20
fx_var = c("poly(inp_rate, 2)2")

mgwr.model <- MGWRSAR(mgwr.fml, ptsdf, 
                      coord = ptsdf, fixed_vars = fx_var, 
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






