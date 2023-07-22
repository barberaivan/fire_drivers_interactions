# In this code I will try to correct the difference between pairs according
# to the distance between them. i.e., remove the distance effect from the
# difference between pixels in each covariate.

# I will compute the absolute difference for each covariate.
# Then, I will estimate which would that difference be if the distance
# would have been the mean distance (using the cumulative probability).
# Then, each pair of predictors values will be put symmetrically closer or
# further apart according to the expected difference based on the minimum
# distance between pairs.
# To avoid moving the predictor values out of range, the corrected difference
# will be computed for the smallest distance; hence, absolute differences will
# only be shrunk.

options(scipen = 999)

# Packages ----------------------------------------------------------------

library(tidyverse); theme_set(theme_bw())
library(viridis)

library(grid)
library(egg)      # has its own ggarrange! much better than ggpubr
library(ggh4x)    # varying strip theme for veg_types and all together

library(mgcv)     # bam()
library(parallel) # for bam in parallel
library(DHARMa)
library(bayestestR)      # highest density intervals
library(rstan)    # hierarchical model
library(survival) # clogit model

library(trialr) # rlkjcorr

# functions ---------------------------------------------------------------

# Normalize
normalize <- function(x) x / sum(x)

softmax <- function(x) exp(x) / sum(exp(x))

# Functions to summaryze
hdi_lower <- function(x, ci = 0.8) {
  samples <- as.numeric(x)
  return(hdi(samples, ci = ci)[["CI_low"]])
}

hdi_upper <- function(x, ci = 0.8) {
  samples <- as.numeric(x)
  return(hdi(samples, ci = ci)[["CI_high"]])
}

hdmean <- function(x, ci = 0.8, name = "mu") {
  ci <- hdi(x, ci = ci)
  result <- c(ci$CI_low, mean(x), ci$CI_high)
  names(result) <- paste(rep(name, 3), c("lower", "mean", "upper"), sep = "_")
  result
}

hdmedian <- function(x, ci = 0.8, name = "mu") {
  ci <- hdi(x, ci = ci)
  result <- c(ci$CI_low, median(x), ci$CI_high)
  names(result) <- paste(rep(name, 3), c("lower", "median", "upper"), sep = "_")
  result
}

quantiles <- function(x, ci = 0.8, name = "mu") {
  p_low <- (1 - ci) / 2
  p_upp <- 1 - p_low
  result <- quantile(x, probs = c(p_low, 0.5, p_upp), method = 8)
  names(result) <- paste(rep(name, 3), c("lower", "median", "upper"), sep = "_")
  result
}

mean_ci <- function(x, ci = 0.8, name = "mu") {
  p_low <- (1 - ci) / 2
  p_upp <- 1 - p_low
  qq <- quantile(x, probs = c(p_low, p_upp), method = 8)
  result <- c(mean(x), qq)
  names(result) <- paste(rep(name, 3), c("mean", "lower", "upper"), sep = "_")
  return(result)
}

# To share legend (not used I think)
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# custom ggplot theme -----------------------------------------------------

# from https://rpubs.com/mclaire19/ggplot2-custom-themes

theme_mine <- function() {
  font <- "Arial"   #assign font family up front
  marg <- 2 # figure margin in mm

  theme_bw() %+replace%    #replace elements we want to change

    theme(

      #grid elements
      #panel.grid.major = element_blank(),    #strip major gridlines
      panel.grid.minor = element_blank(),    #strip minor gridlines
      #axis.ticks = element_blank(),          #strip axis ticks

      #text elements
      plot.title = element_text(             #title
        family = font,            #set font family
        size = 16,                #set font size
        #face = 'bold',            #bold typeface
        hjust = -0.1,                #left align
        vjust = 1),

      # plot.subtitle = element_text(          #subtitle
      #   family = font,            #font family
      #   size = 14),               #font size

      axis.title = element_text(             #axis titles
        family = font,            #font family
        size = 12),

      # para separar el eje y de los nros
      axis.title.y = element_text(
        margin = margin(t = 0, r = 2, b = 0, l = 0, "mm"),
        angle = 90),

      axis.text = element_text(              #axis text
        family = font,            #axis family
        size = 9),                #font size

      legend.title = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 9, family = font),

      strip.text = element_text(size = 12, family = font, color = "white"),
      strip.text.x = element_text(margin = margin(1.2,0,1.2,0, "mm")), # tamaño de la cajita
      strip.text.y = element_text(margin = margin(0,1.2,0,1.2, "mm")),
      strip.background = element_rect(fill = "gray10", color = "gray10"),

      plot.margin = unit(c(marg, marg, marg, marg), "mm")
    )
}

theme_set(theme_mine())


# Load data ----------------------------------------------------------------

d <- read.csv("data_edge-core_filtered.csv")
burned_veg <- read.csv("data_burned_and_available_area_by_vegetation_dryforest2.csv")
# to recode vegetation


# Prepare data ------------------------------------------------------------

names(d)
d$vegetation_code <- d$vegetation

# code dry forest A as subalpine forest
d$vegetation_code[d$vegetation_code == 4] <- 2
# plantation as dry forest B
d$vegetation_code[d$vegetation_code == 9] <- 5
# anthrop as shrubland
d$vegetation_code[d$vegetation_code == 8] <- 6

d <- left_join(d, burned_veg[, c("vegetation_code", "vegetation_class")],
               by = "vegetation_code")

# order by fire and pair_id
d <- d[order(d$fire, d$pair_id), ]
# View(d)

# aggregate data at fire level
data_fires <- aggregate(cbind(area_ha, fwi, dstnc_p) ~ fire_id + year, d, mean)

# factorize fire_id
data_fires$fire_id <- factor(data_fires$fire_id)
d$fire_id <- factor(d$fire_id, levels = levels(data_fires$fire_id))

# numeric response
d$class_bin <- as.numeric(factor(d$class, levels = c("edge", "core"))) - 1

# northing
d$northing <- cos(d$aspect * (pi / 180))

# covariates and metadata names
cov_all <- c(
  "vegetation", "ndvi",
  # "unburnable", # not used
  # "pp", "temp",

  "elevation", "slope", "northing", #"aspect",
  "TPI300", "TPI1k", "TPI2k",
  "HLI_south", "ERR450", "Rough450",

  "dist_human", "dist_roads"
)

# without vegetation class
cov_cont <- cov_all[-1]

meta_names <- c("fire_id", "year", "area_ha", "fwi", "pair_id", "dstnc_p")

# pair_id does not have to be repeated
pair_temp <- paste(d$fire_id, d$pair_id, sep = "//")
pair_temp <- factor(pair_temp)
d$pair_id_all <- as.numeric(pair_temp)

# remove the non burnable pairs
d_burnable <- d[d$vegetation_class != "Non burnable" &
                  d$dstnc_p >= 90, ]
pair_complete <- aggregate(pair_id ~ pair_id_all, d_burnable, length)
pair_complete <- pair_complete[pair_complete$pair_id > 1, ]
d_burnable <- d_burnable[d_burnable$pair_id_all %in% pair_complete$pair_id_all, ]

# keep only fires with enough pixels
enes <- aggregate(ndvi ~ fire_id, d_burnable, length)
quantile(enes$ndvi, seq(0, 1, by = 0.05))
# subset fires with less than 10 pixels
fires_10 <- enes$fire[enes$ndvi >= 10]
d_burnable <- d_burnable[d_burnable$fire_id %in% fires_10, ]


# Compute differences
d_diff <- cbind(
  d_burnable[d_burnable$class == "edge", meta_names],
  d_burnable[d_burnable$class == "core", cov_cont] -
    d_burnable[d_burnable$class == "edge", cov_cont]
)

# Correct differences
diff_correct <- matrix(NA, nrow(d_diff), length(cov_cont))
colnames(diff_correct) <- cov_cont
par(mfrow = c(1, 2))
for(v in cov_cont) {
  # v = "dist_human"
  print(v)
  dd <- d_diff[, c("dstnc_p", v)]
  names(dd) <- c("dist", "y")
  dd$y <- abs(dd$y)

  # replace zeroes and scale
  dd$y[dd$y == 0] <- 1e-9
  dd$dist <- (dd$dist - min(dd$dist)) / sd(dd$dist)

  ggplot(dd, aes(dist, y)) +
    geom_smooth() +
    geom_point()

  # fit model to difference as a function of distance
  mm <- gam(list(y ~ dist, ~ dist), data = dd,
            family = gammals(link = c("identity", "identity")),
            method = "REML")

  # shape = 1 / phi,
  # rate = 1 / (phi * mu)
  # fitted returns the mu and the log_phi
  # (see ?gammals example code)
  fitted_mu <- fitted(mm)[, 1]
  fitted_phi <- fitted(mm)[, 2] %>% exp

  # compute cumulative proabibilities
  pobs <- pgamma(dd$y, shape = 1 / fitted_phi, rate = 1 / (fitted_phi * fitted_mu))
  hist(pobs, breaks = 20, main = v)

  # get counterfactual quantiles
  phi_0 <- exp(coef(mm)[3])
  mu_0 <- exp(coef(mm)[1])

  q_cf <- qgamma(pobs, shape = 1 / phi_0, rate = 1 / (phi_0 * mu_0))

  # check whether corrected is usually smaller
  plot(q_cf ~ dd$y, main = v, col = rgb(0, 0, 0, 0.1), pch = 19,
       ylab = "corrected difference",
       xlab = "observed difference")
  abline(0, 1, col = "red", lwd = 1.5)

  # use the smallest difference
  diff_min <- apply(cbind(dd$y, q_cf), 1, min)
  diff_correct[, v] <- diff_min
}
par(mfrow = c(1, 1))

# shrink predictors so that the absolute difference is corrected
d_mean <- (d_burnable[d_burnable$class == "core", cov_cont] +
           d_burnable[d_burnable$class == "edge", cov_cont]) / 2

signs <- sign(d_diff[, cov_cont])
diff_correct_signed <- diff_correct * signs

# dc is data corrected (diff was computed as core - edge)
dc <- d_burnable
dc[d_burnable$class == "core", cov_cont] <- d_mean + (diff_correct_signed / 2)
dc[d_burnable$class == "edge", cov_cont] <- d_mean - (diff_correct_signed / 2)

# # check
# d_diff_corr <- dc[dc$class == "core", cov_cont] - dc[dc$class == "edge", cov_cont]
# for(v in cov_cont) {
#   plot(d_diff_corr[, v] ~ diff_correct_signed[, v]); abline(0, 1, col = "red")
# }

# # Plot differences as a function of distance
# for(i in cov_cont) {
#   # i = "ndvi"
#   yy <- scale(abs(d_diff_corr[, i])) %>% as.numeric
#   xx <- scale(dc$dstnc_p[dc$class == "core"]) %>% as.numeric
#   mdif <- lm(yy ~ xx)
#   coef_z <- coef(mdif)[2]
#
#   plot(yy ~ xx,
#        main = paste(i, round(unname(coef_z), 4), sep = "; "),
#        ylab = i, xlab = "distance (z)",
#        col = rgb(0, 0, 0, 0.05), pch = 19)
#   abline(coef(mdif), lwd = 2, col = 2)
# }
# differences were detrended

# Variable selection based on univariate models (R2) ----------------------

ncont <- length(cov_cont)
mlist_uni <- vector("list", ncont)
names(mlist_uni) <- cov_cont
r2_uni <- data.frame(var = cov_cont, r2 = NA)

for(i in 1:ncont) {
  # i = 1
  dtemp <- dc[, c("class_bin", "pair_id_all", cov_cont[i])]
  names(dtemp) <- c("y", "pair", "x")
  modelo <- clogit(y ~ x + I(x ^ 2) + strata(pair), data = dtemp)
  mlist_uni[[i]] <- modelo
  p <- predict(modelo, type = "expected")
  r2_uni[i, "r2"] <- var(p) / (var(p) + mean(p * (1 - p)))

  xx <- quantile(dtemp$x, probs = c(0.05, 0.95))
  # xseq <- seq(min(dtemp$x), max(dtemp$x), length.out = 150)
  xseq <- seq(xx[1], xx[2], length.out = 150)
  nd <- data.frame(x = xseq, pair = dtemp$pair[1])
  pp <- predict(modelo, nd, type = "lp")
  plot(exp(pp) ~ xseq, main = cov_cont[i], type = "l")
}

# pairs(dc[, cov_cont], col = rgb(0, 0, 0, 0.01), pch = 19)
cm <- cor(dc[, cov_cont])
# View(round(cm, 4))
# View(round(abs(cm), 4))

r2_uni[order(r2_uni$r2, decreasing = T), ]
candidates <- c(
  "elevation",
  "TPI2k",
  "ndvi",
  "ERR450",
  "slope",
  "Rough450",
  "northing",
  "HLI_south",
  "dist_roads",
  "dist_human"
)

# View(round(abs(cm[candidates, candidates]), 4))
# high correlations:
# TPI2k, ERR
# slope, rough
# slope, HLI_south
# dist_roads, dist_human

# coeffs correlation? removing only ERR and Rough
model_corr_1 <- clogit(
  class_bin ~
    vegetation_class +
    elevation + I(elevation ^ 2) +
    TPI1k + I(TPI1k ^ 2) +
    ndvi + I(ndvi ^ 2) +
    slope +
    HLI_south +
    dist_human +
    dist_roads + I(dist_roads ^ 2) +
    strata(pair_id_all),
  data = dc
)
# View(cov2cor(vcov(model_corr_1))) # OK, they don't change
# View(cov2cor(vcov(model_corr_1))) # OK, they don't change
# coef(model_corr_1)

# vector with selected predictors
predictors_cont <- c(
  "ndvi",
  "elevation",
  "TPI2k",
  "slope",
  "northing",
  "dist_roads",
  "dist_human"
)

# standardize predictors
dc_z <- dc
for(i in predictors_cont) {
  dc_z[, i] <- as.numeric(scale(dc[, i]))
}

# get means and sd to unstandardize
data_scale <- data.frame(
  mean = apply(dc[, predictors_cont], 2, mean),
  sd = apply(dc[, predictors_cont], 2, sd)
) %>% as.matrix %>% t

# summarize predictors for predictions
summaries <- do.call("rbind", lapply(predictors_cont, function(p) {
  mm <- quantiles(dc_z[, p], ci = 0.95)
  mm2 <- matrix(mm, nrow = 1)
  colnames(mm2) <- names(mm)
  mm2 <- as.data.frame(mm2)
  mm2$predictor <- p
  return(mm2)
}))

# get fwi summary
summary_fwi <- quantiles(dc_z$fwi, ci = 0.95, name = "fwi")


# Full model: FWI ---------------------------------------------------------

m_fwi <- clogit(
  class_bin ~
    vegetation_class +
    vegetation_class : fwi +
    elevation + I(elevation ^ 2) +
    elevation : fwi + I(elevation ^ 2) : fwi +
    TPI2k + I(TPI2k ^ 2) +
    TPI2k : fwi + I(TPI2k ^ 2) : fwi +
    ndvi + I(ndvi ^ 2) +
    ndvi : fwi + I(ndvi ^ 2) : fwi +
    slope +
    slope : fwi +
    northing +
    northing : fwi +
    dist_human +
    dist_human : fwi +
    dist_roads + I(dist_roads ^ 2) +
    dist_roads : fwi + I(dist_roads ^ 2) : fwi +

    strata(pair_id_all),
  data = dc_z
)
summary(m_fwi)

# residuals analysis
pfit <- predict(m_fwi, type = "expected")
ysim <- matrix(rbinom(length(pfit) * 3000, prob = pfit, size = 1),
               nrow = length(pfit))
res <- createDHARMa(ysim, observedResponse = dc_z$class_bin)
# plot(res)
# plotResiduals(res, form = dc_z$fire_id)
# plotResiduals(res, form = dc_z$vegetation_class)
# for(p in predictors_cont) plotResiduals(res, form = dc_z[, p], xlab = p)
# everything is OK

# Predictions

# prediction data for continuous predictors

npred <- 100

pdata <- do.call("rbind", lapply(predictors_cont, function(p) {
  # p = "ndvi"
  values_list <- lapply(predictors_cont, function(pred) {
    # pred = "elevation"
    if(pred != p) {
      res <- 0#summaries$mu_median[summaries$predictor == pred]
      # zero to remove them from the linear predictor
    } else {
      res <- seq(summaries$mu_lower[summaries$predictor == pred],
                 summaries$mu_upper[summaries$predictor == pred],
                 length.out = npred)
    }
    return(res)
  })
  names(values_list) <- predictors_cont
  values_list$fwi <- unname(summary_fwi)

  gg <- expand.grid(values_list)
  gg$vegetation_class <- "Subalpine forest"

  gg$varying_name <- p
  gg$varying_value <- as.numeric(gg[, p])

  return(gg)
}))

# compute lp by hand to remove vegetation intercepts
pdata$class_bin <- 1
pdata$pair_id_all <- dc$pair_id[1]
dm <- model.matrix(m_fwi, data = pdata)
# linear predictor removing vegetation intercepts
out_veg <- grep("vegetation", names(coef(m_fwi)))
lp <- dm[, -out_veg] %*% coef(m_fwi)[-out_veg]
pdata$p_mle <- lp # apply softmax later

# Simulate to get CI (softmax impedes the use of se.fit)
coef_sim <- mgcv::rmvn(10000, coef(m_fwi)[-out_veg],
                       vcov(m_fwi)[-out_veg, -out_veg]) %>% t
lp_sim <- dm[, -out_veg] %*% coef_sim

# apply softmax by group
pdata$p_mle <- 0
pdata$p_lower <- 0
pdata$p_upper <- 0

for(f in unique(pdata$fwi)) {
  for(v in predictors_cont) {
    # f = unique(pdata$fwi)[1]
    # v = "ndvi"
    rows <- which(pdata$fwi == f & pdata$varying_name == v)

    probs_sim <- apply(lp_sim[rows, ], 2, softmax)
    probs_ci <- apply(probs_sim, 1, quantile, probs = c(0.025, 0.975), method = 8) %>% t

    pdata$p_lower[rows] <- probs_ci[, 1]
    pdata$p_upper[rows] <- probs_ci[, 2]

    pdata$p_mle[rows] <- softmax(lp[rows])
  }
}

# unstandardize predictors
pdata2 <- pdata
for(p in predictors_cont) {
  # values not used for plotting
  pdata2[, p] <- pdata[, p] * data_scale["sd", p] + data_scale["mean", p]
  # values used for plotting
  filt <- pdata$varying_name == p
  pdata2[filt, "varying_value"] <- pdata[filt, "varying_value"] *
    data_scale["sd", p] +
    data_scale["mean", p]
}

# factorize
pdata2$varying_name <- factor(pdata2$varying_name, levels = predictors_cont)
pdata2$fwi_factor <- factor(as.character(round(pdata2$fwi, 3)),
                            levels = as.character(round(summary_fwi, 3)))

ggplot(data = pdata2,
       mapping = aes(x = varying_value, y = p_mle, ymin = p_lower, ymax = p_upper,
                     colour = fwi_factor, fill = fwi_factor, group = fwi_factor)) +
  geom_ribbon(color = NA, alpha = 0.3) +
  geom_line() +
  facet_wrap(vars(varying_name), scales = "free") +
  scale_color_viridis(option = "B", end = 0.7, discrete = TRUE) +
  scale_fill_viridis(option = "B", end = 0.7, discrete = TRUE) +
  theme(legend.position = c(0.7, 0.15),
        axis.title.x = element_blank(),
        legend.title = element_text()) +
  ylab("Conditional burn probability") +
  labs(fill = "FWI anomaly", color = "FWI anomaly")


# Hierarchical model with stan (FWI) ---------------------------------------

# choose reference vegetation class
enes_veg <- aggregate(ndvi ~ vegetation_class, dc, length)
enes_veg <- enes_veg[order(enes_veg$ndvi, decreasing = TRUE), ]
# make it factor in dc_z
dc_z$vegetation_class <- factor(dc_z$vegetation_class,
                                levels = enes_veg$vegetation_class)
n_veg <- length(levels(dc_z$vegetation_class))

# Prepare data
formula_fixed <- formula(
  ~ vegetation_class +
    ndvi + I(ndvi ^ 2) +
    elevation + I(elevation ^ 2) +
    TPI2k + I(TPI2k ^ 2) +
    slope +
    northing +
    dist_human +
    dist_roads + I(dist_roads ^ 2)
)

# design matrix for data-level predictors
# remove intercept, so shrubland intercept is fixed at zero
X <- model.matrix(formula_fixed, dc_z)[, -1]

# design matrix for group-level predictors
fires_agg_h <- aggregate(fwi ~ fire_id, dc_z, mean)
fires_enes <- enes[enes$fire_id %in% fires_10, ]
# all(fires_agg_h$fire_id == fires_enes$fire_id) # OK
# all(unique(dc_z$fire_id) == fires_enes$fire_id) # OK
Z <- model.matrix(~ fwi,# + I(fwi ^ 2), fires_agg_h)
                  fires_agg_h)

# Indices
N <- nrow(X)  # pixels
K <- ncol(X)  # pixels predictors dimension

NF <- nrow(Z) # fires
J <- ncol(Z)  # fires predictors dimension

N_pairs <- N / 2

# id for start and end of each fire
ends <- cumsum(fires_enes$ndvi)
starts <- c(1, ends[-NF] + 1)
# length(ends) == NF
# cbind(starts, ends)

# ids for core and edge
core_ids <- which(dc_z$class_bin == 1)
edge_ids <- which(dc_z$class_bin == 0)
# all(edge_ids == (core_ids + 1))

sdata <- list(
  N = N, NF = NF, J = J, K = K, N_pairs = N_pairs,
  ends = ends, starts = starts,
  core_ids = core_ids, edge_ids = edge_ids,
  X = X, Z = Z,

  prior_corr_eta = 1,
  prior_sigma_sd = 3,
  prior_b_sd = 5
)

# init function for Lcorr
init_f_random <- function(chain_id) {
  list(
    L_corr = t(chol(rlkjcorr(1, K, eta = 30)))
    # tight corr matrix (near identity)
  )
}

# compile
smodel <- stan_model("interactions model - edge core_01.stan")
# sample
m_fwi_h <- sampling(
  smodel, sdata, refresh = 10,
  control = list(adapt_delta = 0.95),
  cores = 4,
  chains = 4,
  iter = 2500,
  # cores = 1,
  # chains = 1,
  # iter = 10,
  init = init_f_random # initalizes correlation matrix
)
saveRDS(m_fwi_h, "interactions model - edge core - distance corrected_samples.rds")
# vuela. hermoso. 792.786 / 60 = 13.2131 min
sm <- summary(m_fwi_h)[[1]]
summary(sm[, "n_eff"]) # OK; min = 410
summary(sm[, "Rhat"])  # OK; max = 1.0138

# d parameters

# Predictions from hierarchical model.
# Continuous covariates

# Prepare parameters
V_mat <- as.matrix(m_fwi_h, pars = "V")
n_post <- nrow(V_mat)
V <- array(NA, c(K, K, n_post))
for(i in 1:n_post) V[, , i] <- matrix(V_mat[i, ], K, K)

g_mat <- as.matrix(m_fwi_h, pars = "g")
g <- array(NA, c(K, J, n_post))
for(i in 1:n_post) g[, , i] <- matrix(g_mat[i, ], K, J)

# design matrix for fwi at three values
Zt_pred <- t(cbind(rep(1, 3), unname(summary_fwi)))#, unname(summary_fwi) ^ 2))

b <- array(NA, c(K, 3, n_post))
for(i in 1:n_post) b[, , i] <- g[, , i] %*% Zt_pred

# prediction data
pdata_h <- pdata[order(pdata$fwi), ] # order matters!!
pdata_h$vegetation_class <- levels(dc_z$vegetation_class)
X_pred <- model.matrix(formula_fixed, pdata_h)[, -1] # remove intercept
X_pred[, grep("vegetation_class", colnames(X_pred))] <- 0
# View(X_pred)

# split X_pred by fwi
X3 <- array(NA, c(nrow(pdata_h) / 3, K, 3))
for(i in 1:3) X3[, , i] <- X_pred[pdata_h$fwi == summary_fwi[i], ]

n_ran <- 20 # simulated random effects
# predictions array to fill
phat_arr <- array(NA, c(nrow(X3), n_post, 3))
phat_temp <- matrix(NA, nrow(X3), n_ran)

# indexes indicating over which values to take the softmax
predictor_index <- rep(1:length(predictors_cont), each = npred)

for(s in 1:n_post) {
  print(s)
  # s = 1

  # sample random effects
  errors <- mgcv::rmvn(n_ran, rep(0, K), V[, , s]) %>% t
  params <- b[, , s]
  params_ran <- abind::abind(
    lapply(1:n_ran, function(f) {
      params + errors[, f]      ## careful here
    }),
    along = 3
  )

  # Loop over fwi values, random effects and predictors to compute p_hat
  for(v in 1:3) { # fwi values
    for(e in 1:n_ran) { # random effects
      # v = 1; e = 1
      # compute linear predictor
      linpred <- X3[, , v] %*% params_ran[, v, e]

      # turn into probability by predictor group
      for(p in 1:length(predictors_cont)) {
        # p = 1
        filt <- predictor_index == p
        phat_temp[filt, e] <- softmax(linpred[filt])
      }
    }

    # average over random effects and store value
    phat_arr[, s, v] <- rowMeans(phat_temp)
  }
}

# bind fwi values
phat_mat <- do.call(rbind, lapply(1:3, function(i) phat_arr[, , i]))
phat_summ <- apply(phat_mat, 1, mean_ci, name = "p", ci = 0.95) %>% t %>% as.data.frame
sum(phat_summ$p_mean) # perfect

pdata_h_temp <- pdata_h[, c(predictors_cont, "fwi", "varying_name", "varying_value")]
pdata_h_sub <- cbind(pdata_h_temp, phat_summ)

pdata_h_sub$varying_name <- factor(pdata_h_sub$varying_name, levels = predictors_cont)

pdata_h_sub2 <- pdata_h_sub
for(p in predictors_cont) {
  # values not used for plotting
  pdata_h_sub2[, p] <- pdata_h_sub[, p] * data_scale["sd", p] + data_scale["mean", p]
  # values used for plotting
  filt <- pdata_h_sub$varying_name == p
  pdata_h_sub2[filt, "varying_value"] <- pdata_h_sub[filt, "varying_value"] *
    data_scale["sd", p] +
    data_scale["mean", p]
}

pdata_h_sub2$fwi_factor <- factor(pdata_h_sub2$fwi)

ggplot(data = pdata_h_sub2,
       mapping = aes(x = varying_value, y = p_mean, ymin = p_lower, ymax = p_upper,
                     colour = fwi_factor, fill = fwi_factor, group = fwi_factor)) +
  geom_ribbon(color = NA, alpha = 0.3) +
  geom_line() +
  facet_wrap(vars(varying_name), scales = "free") +
  scale_color_viridis(option = "B", end = 0.7, discrete = TRUE) +
  scale_fill_viridis(option = "B", end = 0.7, discrete = TRUE) +
  theme(legend.position = c(0.7, 0.15),
        axis.title.x = element_blank(),
        legend.title = element_text()) +
  ylab("Conditional burn probability") +
  labs(fill = "FWI anomaly", color = "FWI anomaly")
ggsave("figures/climate-fuel-interactions_clogis_m1_hierarchical-at-mean_dc.png",
       width = 17, height = 15, units = "cm")

# Las medias son muy afectadas por los valores extremos, entonces me parece
# mejor mirar la mediana. (Los intervalos de las medias se van al diablo)


## Vegetation plot

pdata_veg <- aggregate(cbind(ndvi, elevation, TPI2k, slope, northing,
                             dist_roads, dist_human) ~ vegetation_class,
                       dc_z, mean)
dm_veg <- model.matrix(formula_fixed, pdata_veg)[, -1]

# array to fill
phat_arr_veg <- phat_arr[1:n_veg, , ]
for(s in 1:n_post) {
  print(s)
  # s = 1
  params <- b[, , s]
  linpred <- dm_veg %*% params
  phat_arr_veg[, s, ] <- apply(linpred, 2, softmax)
}

pveg <- apply(phat_arr_veg, c(1, 3), mean_ci, name = "p", ci = 0.95)
dimnames(pveg) <- list(summ = dimnames(pveg)[[1]],
                       vegetation_class = levels(dc_z$vegetation_class),
                       fwi_factor = levels(pdata_h_sub2$fwi_factor))
pveg_tab <- as.data.frame.table(pveg)
pveg_wide <- pivot_wider(pveg_tab, names_from = "summ", values_from = "Freq")

ggplot(pveg_wide, aes(x = vegetation_class, y = p_mean,
                      ymin = p_lower, ymax = p_upper,
                      colour = fwi_factor)) +
  geom_point(size = 2) +
  geom_errorbar(alpha = 0.8, linewidth = 0.5, width = 0.2) +
  facet_wrap(vars(fwi_factor)) +
  scale_color_viridis(option = "B", end = 0.7, discrete = TRUE) +
  scale_fill_viridis(option = "B", end = 0.7, discrete = TRUE) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position = "none")
ggsave("figures/climate-fuel-interactions_clogis_m1_hierarchical-at-median_veg_dc.png",
       width = 17, height = 10, units = "cm")

