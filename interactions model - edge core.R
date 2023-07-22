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
data_fires <- aggregate(cbind(area_ha, fwi) ~ fire_id + year, d, mean)

# factorize fire_id
data_fires$fire_id <- factor(data_fires$fire_id)
d$fire_id <- factor(d$fire_id, levels = levels(data_fires$fire_id))


# covariates and metadata names
cov_all <- c(
  "vegetation", "ndvi", "unburnable",
  "pp", "temp",

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


# Plot differences --------------------------------------------------------

for(i in cov_cont) {
  pp <- sum(as.numeric(d_diff[, i] > 0)) / nrow(d_diff)
  hist(d_diff[, i], breaks = 50, main = paste(i, pp, sep = ";"));
  abline(v = 0, col = 2, lwd = 2)
}


# Plot differences as a function of distance
for(i in cov_cont) {
  # i = "ndvi"
  yy <- scale(abs(d_diff[, i])) %>% as.numeric
  xx <- scale(d_diff$dstnc_p) %>% as.numeric
  mdif <- lm(yy ~ xx)
  coef_z <- coef(mdif)[2]

  plot(yy ~ xx,
       main = paste(i, round(unname(coef_z), 4), sep = "; "),
       ylab = i, xlab = "distance (z)",
       col = rgb(0, 0, 0, 0.05), pch = 19)
  abline(coef(mdif), lwd = 2, col = 2)
}
# In general, differences increase little with distance between points

# Plot differences in large fires
ff_use <- data_fires$fire_id[data_fires$area_ha > 400]
ff_use <- data_fires$fire_id[data_fires$area_ha < 1000]
diff_large <- d_diff[d_diff$fire_id %in% ff_use, ]

for(i in cov_cont) {
  # pp <- sum(as.numeric(diff_large[, i] > 0)) / nrow(diff_large)
  pp <- mean(diff_large[, i]) / sd(diff_large[, i])
  hist(diff_large[, i], breaks = 50, main = paste(i, pp, sep = ";"));
  abline(v = 0, col = 2, lwd = 2)
}
# differences seem to be small, but a bit larger in large fires


# overall distributions
for(i in cov_cont) {
  cc <- d[d$class == "core", i]
  ee <- d[d$class == "edge", i]

  dcc <- density(cc)
  dee <- density(ee)

  xx <- range(c(cc, ee))
  yy <- c(0, max(dcc$y, dee$y))

  plot(dcc$x, dcc$y, type = "l", xlim = xx, ylim = yy, col = 2,
       main = i)
  lines(dee$x, dee$y)
}
# hay varios efectos cuadráticos, por eso muchas diferencias no se notan.




# Conditional regression is different -------------------------------------
# and we could code it in stan

# https://github.com/dcmuller/stan_clogit/blob/master/clogit.stan
# https://github.com/paul-buerkner/brms/issues/560

# http://mc-stan.org/rstanarm/reference/stan_clogit.html
# http://mc-stan.org/rstanarm/articles/binomial.html

# https://sci-hub.ru/https://doi.org/10.1093/biomet/68.3.703

# Kruske tiene una buena explicación, pero no sé si la veo aplicable.

# library(rstanarm) # has clogit() or something equivalent
# example(example_model)
# rstan::get_stanmodel(example_model$stanfit) # no es muy legible, pero bue

# por ahora, ajustarlo con el paq survival y luego ver qué onda.
# Luego vemos si hace falta bayesianear.

# All predictors interacting with FWI or log(area_ha)


# Fitting the model with survival (fwi) ------------------------------------

# pair_id does not have to be repeated
pair_temp <- paste(d$fire_id, d$pair_id, sep = "//")
pair_temp <- factor(pair_temp)
d$pair_id_all <- as.numeric(pair_temp)

d_burnable <- d[d$vegetation_class != "Non burnable", ]
pair_complete <- aggregate(pair_id ~ pair_id_all, d_burnable, length)
pair_complete <- pair_complete[pair_complete$pair_id > 1, ]
d_burnable <- d_burnable[d_burnable$pair_id_all %in% pair_complete$pair_id_all, ]
nrow(d_burnable)

# try simple model: ndvi, tpi300, and HLI
m1 <- clogit(
  class_bin ~
    vegetation_class + vegetation_class : fwi +
    ndvi + I(ndvi ^ 2) + ndvi : fwi + I(ndvi ^ 2) : fwi +
    TPI300 + I(TPI300 ^ 2) + TPI300 : fwi + I(TPI300 ^ 2) : fwi +
    HLI_south + HLI_south : fwi +
    elevation + I(elevation ^ 2) + elevation : fwi + I(elevation ^ 2) : fwi +    strata(pair_id_all),
  data = d_burnable
)
summary(m1)



# Predictions
quad <- TRUE
pred <- "elevation"
data <- d_burnable
exts <- c(min(data[, pred]), max(data[, pred]))
fwi_vals <- quantile(data$fwi, prob = c(0.1, 0.5, 0.9))

# get HDI and median by predictor and by vegetation types
predictors_cont <- c("ndvi",  # cont for continuous
                     "TPI300", "elevation", "HLI_south")
# quadratic effect for predictors?
predictors_quad <- c(TRUE, TRUE, TRUE, FALSE)
names(predictors_quad) <- predictors_cont

summaries <- do.call("rbind", lapply(predictors_cont, function(p) {
  mm <- quantiles(data[, p], ci = 0.95)
  mm2 <- matrix(mm, nrow = 1)
  colnames(mm2) <- names(mm)
  mm2 <- as.data.frame(mm2)
  mm2$predictor <- p
  return(mm2)
}))

# get fwi summary
summary_fwi <- quantiles(data$fwi, ci = 0.95, name = "fwi")

# prediction data for continuous predictors
pdata <- do.call("rbind", lapply(predictors_cont, function(p) {
  # p = "ndvi"
  values_list <- lapply(predictors_cont, function(pred) {
    # pred = "elevation"
    if(pred != p) {
      res <- summaries$mu_median[summaries$predictor == pred]
    } else {
      res <- seq(summaries$mu_lower[summaries$predictor == pred],
                 summaries$mu_upper[summaries$predictor == pred],
                 length.out = 200)
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

# add strata to use model.matrix()
pdata$pair_id_all <- 1
# get design matrix
dm <- model.matrix(m1, data = pdata)

# compute linear predictor
linpred <- dm %*% coef(m1)

# for each covariate and fwi value, repeat the linear predictor value
# corresponding to the mean or maximum linear predictor, depending on the
# effect being quadratic or not
linpred_ref <- linpred

# add x value at the pdata for plotting
pdata$ref_value <- NA

for(p in predictors_cont) {
  for(v in unique(pdata$fwi)) {
    # p = "ndvi"
    # v = max(summary_fwi)

    target_rows <- which(pdata$fwi == v & pdata$varying_name == p)
    # target_rows <- which(pdata$varying_name == p)

    # define reference value depending on quadratic effect or not
    if(predictors_quad[p]) {
      row_reference <- target_rows[which.max(linpred[target_rows])]
    } else {
      midpoint <- ceiling(length(target_rows) / 2)
      row_reference <- target_rows[midpoint]
    }

    lp_reference <- linpred[row_reference]
    pred_reference <- pdata$varying_value[row_reference]

    linpred_ref[target_rows] <- lp_reference
    pdata$ref_value[target_rows] <- pred_reference
  }
}

# compute probability
pdata$prob_mle <- apply(cbind(linpred, linpred_ref), 1, softmax) %>% t %>% "["(, 1)

ggplot(pdata, aes(x = varying_value, y = prob_mle, colour = fwi, group = fwi)) +
  geom_line() +
  facet_wrap(vars(varying_name), scales = "free") +
  scale_color_viridis(option = "B", end = 0.7)


# Fitting the model with survival (area_ha) -------------------------------

d_burnable$log_area <- log(d_burnable$area_ha)

m2 <- clogit(
  class_bin ~
    vegetation_class + vegetation_class : log_area +
    ndvi + I(ndvi ^ 2) + ndvi : log_area + I(ndvi ^ 2) : log_area +
    TPI300 + I(TPI300 ^ 2) + TPI300 : log_area + I(TPI300 ^ 2) : log_area +
    HLI_south + HLI_south : log_area +
    elevation + I(elevation ^ 2) + elevation : log_area + I(elevation ^ 2) : log_area +

    # interactions with distance between pixels
    vegetation_class : dstnc_p + vegetation_class : log_area : dstnc_p +
    ndvi : dstnc_p + I(ndvi ^ 2) : dstnc_p + ndvi : log_area : dstnc_p + I(ndvi ^ 2) : log_area : dstnc_p +
    TPI300 : dstnc_p + I(TPI300 ^ 2) : dstnc_p + TPI300 : log_area : dstnc_p + I(TPI300 ^ 2) : log_area : dstnc_p +
    HLI_south : dstnc_p + HLI_south : log_area : dstnc_p +
    elevation : dstnc_p + I(elevation ^ 2) : dstnc_p + elevation : log_area : dstnc_p + I(elevation ^ 2) : log_area : dstnc_p +

    strata(pair_id_all),
  data = d_burnable
)
summary(m2)


log_area_vals <- log(c(1000, 5000, 20000))

# get log_area summary
summary_log_area_nn <- quantiles(data$log_area, ci = 0.95, name = "log_area") %>% names
summary_log_area <- log_area_vals
names(summary_log_area) <- summary_log_area_nn

# prediction data for continuous predictors
pdata <- do.call("rbind", lapply(predictors_cont, function(p) {
  # p = "ndvi"
  values_list <- lapply(predictors_cont, function(pred) {
    # pred = "elevation"
    if(pred != p) {
      res <- summaries$mu_median[summaries$predictor == pred]
    } else {
      res <- seq(summaries$mu_lower[summaries$predictor == pred],
                 summaries$mu_upper[summaries$predictor == pred],
                 length.out = 200)
    }
    return(res)
  })
  names(values_list) <- predictors_cont
  values_list$log_area <- unname(summary_log_area)

  gg <- expand.grid(values_list)
  gg$vegetation_class <- "Subalpine forest"
  gg$dstnc_p <- 60

  gg$varying_name <- p
  gg$varying_value <- as.numeric(gg[, p])

  return(gg)
}))

# add strata to use model.matrix()
pdata$pair_id_all <- 1
# get design matrix
dm <- model.matrix(m2, data = pdata)

# compute linear predictor
out <- which(is.na(coef(m2)))
linpred <- dm[, -out] %*% coef(m2)[-out]

# for each covariate and log_area value, repeat the linear predictor value
# corresponding to the mean or maximum linear predictor, depending on the
# effect being quadratic or not
linpred_ref <- linpred

# add x value at the pdata for plotting
pdata$ref_value <- NA

for(p in predictors_cont) {
  for(v in unique(pdata$log_area)) {
    # p = "ndvi"
    # v = max(summary_log_area)

    target_rows <- which(pdata$log_area == v & pdata$varying_name == p)
    # target_rows <- which(pdata$varying_name == p)

    # define reference value depending on quadratic effect or not
    if(predictors_quad[p]) {
      row_reference <- target_rows[which.max(linpred[target_rows])]
    } else {
      midpoint <- ceiling(length(target_rows) / 2)
      row_reference <- target_rows[midpoint]
    }

    lp_reference <- linpred[row_reference]
    pred_reference <- pdata$varying_value[row_reference]

    linpred_ref[target_rows] <- lp_reference
    pdata$ref_value[target_rows] <- pred_reference
  }
}

# compute probability
pdata$prob_mle <- apply(cbind(linpred, linpred_ref), 1, softmax) %>% t %>% "["(, 1)

ggplot(pdata, aes(x = varying_value, y = prob_mle, colour = log_area, group = log_area)) +
  geom_line() +
  facet_wrap(vars(varying_name), scales = "free_x") +
  scale_color_viridis(option = "B", end = 0.7)

# Including interactions with distance and fixing the distance at 90 for predictions,
# yields very similar results, although they change a little bit.


# Difference importance in clogis -----------------------------------------

n1 <- 100
x1 <- rnorm(n1, 1, 1)
x2 <- rep(0, n1)
factor_diff <- 2

data_small <- data.frame(
  y = rep(c(1, 0), each = n1),
  id = rep(1:n1, 2),
  x = c(x1, x2)
)

data_large <- data.frame(
  y = rep(c(1, 0), each = n1),
  id = rep(1:n1, 2),
  x = c(x1, x2) * factor_diff
)

sum(abs(data_small$x))
sum(abs(data_large$x))

p1 <- clogit(y ~ x + strata(id), data = data_small)
p2 <- clogit(y ~ x + strata(id), data = data_large)

coef(p1)
coef(p2)

# Always, when x has larger differences, its coefficients are smaller.
# It's like saying that to have {0, 1} one can allow very different values.
# We should try a fixed buffer distance.

bd <- aggregate(bufdist ~ fire_id, d, mean)
hist(bd$bufdist)
plot(ecdf(bd$bufdist
         ))


dvecs <- terra::vect("data_edge-core_full.shp")
bd <- aggregate(cbind(bufdist, dstnc_p, fwi, area_ha) ~ fire_id, dvecs, median)
bd <- bd[order(bd$bufdist), ]
sum(bd$bufdist > 60)
sum(bd$bufdist > 90)

bdlarge <- bd[bd$bufdist > 90, ]
plot(area_ha ~ fwi, bdlarge)
ggplot(bdlarge, aes(y = area_ha, x = fwi)) +
  geom_smooth() +
  geom_point()

# use fixed buffer at 90 m