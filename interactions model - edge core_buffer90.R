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

d <- read.csv("data_edge-core_filtered_buffer90.csv")
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

# covariates and metadata names
cov_all <- c(
  "vegetation", "ndvi",
  # "unburnable", # not used
  # "pp", "temp",

  "elevation", "slope", "aspect",
  "TPI300", "TPI1k", "TPI2k",
  "HLI_south", "ERR450", "Rough450",

  "dist_human", "dist_roads"
)

# without vegetation class
cov_cont <- cov_all[-1]

meta_names <- c("fire_id", "year", "area_ha", "fwi", "pair_id", "dstnc_p")

# Compute differences
d_diff <- cbind(
  d[d$class == "edge", meta_names],
  d[d$class == "core", cov_cont] - d[d$class == "edge", cov_cont]
)


# pair_id does not have to be repeated
pair_temp <- paste(d$fire_id, d$pair_id, sep = "//")
pair_temp <- factor(pair_temp)
d$pair_id_all <- as.numeric(pair_temp)

# remove the non burnable pairs
d_burnable <- d[d$vegetation_class != "Non burnable", ]
pair_complete <- aggregate(pair_id ~ pair_id_all, d_burnable, length)
pair_complete <- pair_complete[pair_complete$pair_id > 1, ]
d_burnable <- d_burnable[d_burnable$pair_id_all %in% pair_complete$pair_id_all, ]

## USE d_burnable

# Plot differences as a function of distance -------------------------------

# Plot differences as a function of distance
for(i in cov_cont) {
  # i = "ndvi"
  yy <- scale(abs(d_diff[, i])) %>% as.numeric
  xx <- scale(d_diff$dstnc_p) %>% as.numeric
  mdif <- lm(yy ~ xx)
  coef_z <- coef(mdif)[2]

  # plot(yy ~ xx,
  #      main = paste(i, round(unname(coef_z), 4), sep = "; "),
  #      ylab = i, xlab = "distance (z)",
  #      col = rgb(0, 0, 0, 0.05), pch = 19)
  # abline(coef(mdif), lwd = 2, col = 2)
}
# In general, differences increase little with distance between points

# # distance as a function of fire area
# plot(dstnc_p ~ log(area_ha), data_fires)
# plot(dstnc_p ~ area_ha, data_fires)
# even with fixed buffer, distance increases with fire area

# and the median distance?
data_fires_med <- aggregate(cbind(area_ha, fwi, dstnc_p) ~ fire_id + year, d, median)
# plot(dstnc_p ~ log(area_ha), data_fires_med)


# Variable selection based on univariate models (R2) ----------------------

ncont <- length(cov_cont)
mlist_uni <- vector("list", ncont)
names(mlist_uni) <- cov_cont
r2_uni <- data.frame(var = cov_cont, r2 = NA)

for(i in 1:ncont) {
  # i = 1
  dtemp <- d_burnable[, c("class_bin", "pair_id_all", cov_cont[i])]
  names(dtemp) <- c("y", "pair", "x")
  modelo <- clogit(y ~ x + I(x ^ 2) + strata(pair), data = dtemp)
  mlist_uni[[i]] <- modelo
  p <- predict(modelo, type = "expected")
  r2_uni[i, "r2"] <- var(p) / (var(p) + mean(p * (1 - p)))

  xx <- quantile(dtemp$x, probs = c(0.05, 0.95))
  # xseq <- seq(min(dtemp$x), max(dtemp$x), length.out = 150)
  xseq <- seq(xx[1], xx[2], length.out = 150)
  nd <- data.frame(x = xseq, pair = rep(1, 150))
  pp <- predict(modelo, nd, type = "lp")
  # plot(exp(pp) ~ xseq, main = cov_cont[i], type = "l")
}

# pairs(d_burnable[, cov_cont], col = rgb(0, 0, 0, 0.01), pch = 19)
cm <- cor(d_burnable[, cov_cont])
# View(round(cm, 4))
# View(round(abs(cm), 4))

r2_uni[order(r2_uni$r2, decreasing = T), ]
candidates <- c(
  "elevation",
  "TPI1k",
  "ndvi",
  "ERR450",
  "slope",
  "HLI_south",
  "Rough450",
  "dist_roads",
  "dist_human"
)

# View(round(abs(cm[candidates, candidates]), 4))
# high correlations:
# TPI1k, ERR
# slope, rough
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
  data = d_burnable
)
# View(cov2cor(vcov(model_corr_1))) # OK, they don't change
# View(cov2cor(vcov(model_corr_1))) # OK, they don't change
# coef(model_corr_1)

# vector with selected predictors
predictors_cont <- c(
  "ndvi",
  "elevation",
  "TPI1k",
  "slope",
  "HLI_south",
  "dist_roads",
  "dist_human"
)


# standardize predictors
d_burnable_z <- d_burnable
for(i in predictors_cont) {
  d_burnable_z[, i] <- as.numeric(scale(d_burnable[, i]))
}

# get means and sd to unstandardize
data_scale <- data.frame(
  mean = apply(d_burnable[, predictors_cont], 2, mean),
  sd = apply(d_burnable[, predictors_cont], 2, sd)
) %>% as.matrix %>% t

# summarize predictors for predictions
summaries <- do.call("rbind", lapply(predictors_cont, function(p) {
  mm <- quantiles(d_burnable_z[, p], ci = 0.95)
  mm2 <- matrix(mm, nrow = 1)
  colnames(mm2) <- names(mm)
  mm2 <- as.data.frame(mm2)
  mm2$predictor <- p
  return(mm2)
}))

# get fwi summary
summary_fwi <- quantiles(d_burnable_z$fwi, ci = 0.95, name = "fwi")


# Full model: FWI ---------------------------------------------------------

m_fwi <- clogit(
  class_bin ~

    vegetation_class +
      vegetation_class : fwi +
    elevation + I(elevation ^ 2) +
      elevation : fwi + I(elevation ^ 2) : fwi +
    TPI1k + I(TPI1k ^ 2) +
      TPI1k : fwi + I(TPI1k ^ 2) : fwi +
    ndvi + I(ndvi ^ 2) +
      ndvi : fwi + I(ndvi ^ 2) : fwi +
    slope +
      slope : fwi +
    HLI_south +
      HLI_south : fwi +
    dist_human +
      dist_human : fwi +
    dist_roads + I(dist_roads ^ 2) +
      dist_roads : fwi + I(dist_roads ^ 2) : fwi +

    strata(pair_id_all),
  data = d_burnable_z
)
summary(m_fwi)

# residuals analysis
pfit <- predict(m_fwi, type = "expected")
ysim <- matrix(rbinom(length(pfit) * 3000, prob = pfit, size = 1),
               nrow = length(pfit))
res <- createDHARMa(ysim, observedResponse = d_burnable_z$class_bin)
# plot(res)
# plotResiduals(res, form = d_burnable_z$fire_id)
# plotResiduals(res, form = d_burnable_z$vegetation_class)
# for(p in predictors_cont) plotResiduals(res, form = d_burnable_z[, p], xlab = p)
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
pdata$pair_id_all <- 1
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
  labs(fill = "FWI anomaly", color = "FWI anomaly")# guides(fill = "FWI", color = "FWI")pdata2$fwi
ggsave("figures/climate-fuel-interactions_clogis_m1_complete-pooling.png",
       width = 17, height = 15, units = "cm")


############

# Veg plot
pdata_veg <- aggregate(cbind(ndvi, elevation, TPI1k, slope, HLI_south,
                            dist_roads, dist_human) ~ vegetation_class,
                       d_burnable_z, mean)
pdata_veg <- do.call(rbind, lapply(1:3, function(f) {
  pp <- pdata_veg
  pp$fwi <- summary_fwi[f]
  return(pp)
}))

pdata_veg$pair_id_all <- 1
dm_veg <- model.matrix(m_fwi, pdata_veg)

# linpred mle and simulated
lp_veg <- dm_veg %*% coef(m_fwi)
coef_sim_veg <- mgcv::rmvn(10000, coef(m_fwi),
                       vcov(m_fwi)) %>% t
lp_sim_veg <- dm_veg %*% coef_sim_veg


# apply softmax by group
pdata_veg$p_mle <- 0
pdata_veg$p_lower <- 0
pdata_veg$p_upper <- 0

for(f in unique(pdata_veg$fwi)) {
  rows <- which(pdata_veg$fwi == f)

  probs_sim <- apply(lp_sim_veg[rows, ], 2, softmax)
  probs_ci <- apply(probs_sim, 1, quantile, probs = c(0.025, 0.975), method = 8) %>% t

  pdata_veg$p_lower[rows] <- probs_ci[, 1]
  pdata_veg$p_upper[rows] <- probs_ci[, 2]

  pdata_veg$p_mle[rows] <- softmax(lp_veg[rows])
}

ggplot(pdata_veg, aes(x = vegetation_class,
                  y = p_mle, ymin = p_lower, ymax = p_upper,
                  colour = fwi)) +
  geom_errorbar(position = position_dodge()) +
  geom_point() +
  facet_wrap(vars(fwi), ncol = 2) +
  scale_color_viridis(option = "B", end = 0.7) +
  scale_fill_viridis(option = "B", end = 0.7)

# Horribleeeee



# Hierarchical model with stan (FWI) ---------------------------------------

enes <- aggregate(ndvi ~ fire_id, d_burnable, length)
quantile(enes$ndvi, seq(0, 1, by = 0.05))
# subset fires with less than 10 pixels
fires_10 <- enes$fire[enes$ndvi >= 10]
d_burnable_zh <- d_burnable_z[d_burnable_z$fire_id %in% fires_10, ]
nrow(d_burnable_z); nrow(d_burnable_zh) # 38 pixels removed

# choose reference vegetation class
enes_veg <- aggregate(ndvi ~ vegetation_class, d_burnable, length)
enes_veg <- enes_veg[order(enes_veg$ndvi, decreasing = TRUE), ]
# make it factor in d_burnable_z
d_burnable_z$vegetation_class <- factor(d_burnable_z$vegetation_class,
                                        levels = enes_veg$vegetation_class)

# Prepare data
formula_fixed <- formula(
  ~ vegetation_class +
    elevation + I(elevation ^ 2) +
    TPI1k + I(TPI1k ^ 2) +
    ndvi + I(ndvi ^ 2) +
    slope +
    HLI_south +
    dist_human +
    dist_roads + I(dist_roads ^ 2)
  )

# design matrix for data-level predictors
# remove intercept, so shrubland intercept is fixed at zero
X <- model.matrix(formula_fixed, d_burnable_zh)[, -1]

# design matrix for group-level predictors
fires_agg_h <- aggregate(fwi ~ fire_id, d_burnable_zh, mean)
fires_enes <- enes[enes$fire_id %in% fires_10, ]
# all(fires_agg_h$fire_id == fires_enes$fire_id) # OK
# all(unique(d_burnable_zh$fire_id) == fires_enes$fire_id) # OK
Z <- model.matrix(~ fwi, fires_agg_h)

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
core_ids <- which(d_burnable_zh$class_bin == 1)
edge_ids <- which(d_burnable_zh$class_bin == 0)
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
  init = init_f_random # initalizes correlation matrix
)
# vuela. hermoso.
sm <- summary(m_fwi_h)[[1]]
summary(sm[, "n_eff"]) # OK; min = 407
summary(sm[, "Rhat"])  # OK; max = 1.0096

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
g[, , 1]

# design matrix for fwi at three values
Zt_pred <- t(cbind(rep(1, 3), unname(summary_fwi)))

b <- array(NA, c(K, 3, n_post))
for(i in 1:n_post) b[, , i] <- g[, , i] %*% Zt_pred
b[, , 1]

# prediction data
pdata_h <- pdata[order(pdata$fwi), ] # order matters!!
pdata_h$vegetation_class <- levels(d_burnable_zh$vegetation_class)
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
      params #+ errors[, f]      ## careful here
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
phat_summ <- apply(phat_mat, 1, mean_ci, name = "p") %>% t %>% as.data.frame
sum(phat_summ$p_mean) # perfect

pdata_h_temp <- pdata_h[, c(predictors_cont, "fwi", "varying_name", "varying_value")]
pdata_h_sub <- cbind(pdata_h_temp, phat_summ)

pdata_h_sub$fwi_factor <- pdata2$fwi_factor
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
ggsave("figures/climate-fuel-interactions_clogis_m1_hierarchical-at-median.png",
       width = 17, height = 15, units = "cm")

# Results are odd. TPI1k has a linear effect instead of quadratic.
# could it be because of uncertainty?
# make the same plot using the posterior mode
sample_map <- which.max(as.matrix(m_fwi_h, "lp__"))
# are TPI1k and elevation random effects correlated?
# maybe that's the explanation

# también mirar las corr posteriores entre los eff aleatorios
#
# Parece que es por la ineq de Jensen!! Si miro las medianas da muy similar al
# de eff fijos