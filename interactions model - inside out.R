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

library(scales)   # log scale
library(circular) # density.circular, for aspect
library(brms)     # not used yet

# functions ---------------------------------------------------------------

# Normalize
normalize <- function(x) x / sum(x)

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

# We will later use prior correction to compute annual fire probabilities,
# following King and Zeng 2001.
# tau is the population proportion, and phi is the sample proportion.
population_intercept <- function(intercept, tau, phi) {
  b0 <- intercept - log(((1 - tau) / tau) * (phi / (1 - phi)))
  return(b0)
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

# Inside-outside analysis -------------------------------------------------

# for now, just use the large raw dataset
d <- read.csv("data_inside-outside.csv")
nrow(d)
d_year <- aggregate(burned ~ year, d, mean)
d_year_sum <- aggregate(burned ~ year, d, sum)
sum(d_year_sum$burned) / nrow(d) # 0.002670833 es la annual fire prob

# Parece que hay muchos más en burned porque los NA en las predictoras en
# el altoandino no permiten que haya la misma cantidad de puntos.

# recode vegetation and remove unuseful categories
burned_veg <- read.csv("data_burned_and_available_area_by_vegetation_dryforest2.csv")
d$vegetation_code <- d$vegetation_type

# code dry forest A as subalpine forest
d$vegetation_code[d$vegetation_code == 4] <- 2
# plantation as dry forest B
d$vegetation_code[d$vegetation_code == 9] <- 5
# anthrop as shrubland
d$vegetation_code[d$vegetation_code == 8] <- 6

d <- left_join(d, burned_veg[, c("vegetation_code", "vegetation_class")],
               by = "vegetation_code")

unique(d$vegetation_class)

veg_levels <- c(
  "Wet forest",
  "Subalpine forest",
  "Dry forest B",
  "Shrubland",
  "Steppe and grassland"
)

veg_labels <- c(
  "Wet forest",
  "Subalpine\nforest",
  "Dry forest",
  "Shrubland",
  "Steppe and\ngrassland"
)

d$vegetation_class <- factor(d$vegetation_class, levels = veg_levels,
                             labels = veg_labels)

# get nn points by class and fire
fires_id <- d$fire_id %>% unique
data <- do.call("rbind", lapply(fires_id, function(f) {
  dd <- d[d$fire_id == f, ]
  ids_in <- which(dd$burned == 1)
  ids_out <- which(dd$burned == 0)
  l_in <- length(ids_in)
  l_out <- length(ids_out)
  l_min <- min(l_in, l_out)
  n <- ifelse(l_min < 50, l_min, 50)
  keep_in <- sample(ids_in, n)
  keep_out <- sample(ids_out, n)
  return(dd[c(keep_in, keep_out), ])
}))

# Fit model
cl <- makeCluster(10, type = "FORK")
m4 <- bam(
  burned ~
    vegetation_class + vegetation_class : fwi +
    s(ndvi, by = vegetation_class, k = 5, bs = "cr") +
      ti(ndvi, fwi, k = 3) +
    s(elevation, by = vegetation_class, k = 5, bs = "cr") +
      ti(elevation, fwi, k = 3) +
    s(slope, by = vegetation_class, k = 5, bs = "cr") +
      ti(slope, fwi, k = 3) +
    s(aspect, by = vegetation_class, k = 5, bs = "cc") +
      ti(aspect, fwi, k = 5, bs = c("cc", "cr")) +
    s(dist_humans, by = vegetation_class, k = 5, bs = "cr") +
      ti(dist_humans, fwi, k = 3, bs = "cr") +
    s(dist_roads, by = vegetation_class, k = 5, bs = "cr") +
      ti(dist_roads, fwi, k = 3, bs = "cr"),
  knots = list(aspect = c(0, 360)),
  data = data, family = binomial(link = "logit"), method = "REML",
  cluster = cl
)
saveRDS(m4, "burn_probability_model_bam_m4.rds")
m4 <- readRDS("burn_probability_model_bam_m4.rds")

# get HDI and median by predictor and by vegetation types
predictors_cont <- c("ndvi",  # cont for continuous
                     "elevation", "slope", "aspect",
                     "dist_humans", "dist_roads")

# use 100000 points from d (marginal balanced) to get summaries
dsub <- d[sample(1:nrow(d), size = 1e5, replace = F), ]

summaries <- do.call("rbind", lapply(predictors_cont, function(p) {
  agg <- aggregate(dsub[, p] ~ vegetation_class, dsub, quantiles) # change hdmedian for quantiles
  agg <- cbind(agg$vegetation_class, as.data.frame(agg$`dsub[, p]`))
  colnames(agg)[1] <- "vegetation_class"
  agg$predictor <- p
  return(agg)
}))

# correct limits for aspect
summaries$mu_lower[summaries$predictor == "aspect"] <- 0
summaries$mu_upper[summaries$predictor == "aspect"] <- 360

# get fwi summary
summary_fwi <- quantiles(d$fwi, name = "fwi")

# prediction data for continuous predictors
pdata <- do.call("rbind", lapply(predictors_cont, function(p) {
  # p = "ndvi"

  dd <- do.call("rbind", lapply(veg_labels, function(v) {
    # v = "Shrubland"

    # filter data
    ss <- summaries[summaries$vegetation_class == v, ]

    values_list <- lapply(predictors_cont, function(pred) {
      # pred = "slope"
      if(pred != p) {
        res <- ss$mu_median[ss$predictor == pred]
      } else {
        res <- seq(ss$mu_lower[ss$predictor == pred],
                   ss$mu_upper[ss$predictor == pred],
                   length.out = 200)
      }
      return(res)
    })
    names(values_list) <- predictors_cont
    values_list$fwi <- unname(summary_fwi)

    gg <- expand.grid(values_list)
    gg$vegetation_class <- v
    return(gg)
  }))

  dd$varying_name <- p
  dd$varying_value <- as.numeric(dd[, p])

  return(dd)
}))

# tidy names
pdata$predictor <- factor(pdata$varying_name,
                          levels = predictors_cont,
                          labels = c(
                            "NDVI", "Elevation", "Slope",
                            "Aspect", "Dist. sett.",
                            "Dist. roads"
                          ))
pdata$vegetation_class <- factor(pdata$vegetation_class,
                                 levels = veg_labels)

# compute predictions
pcont <- predict(m4, pdata, se.fit = T)

# probs
pdata$p_mle <- plogis(pcont$fit)
pdata$p_lower <- plogis(pcont$fit + qnorm(0.025) * pcont$se.fit)
pdata$p_upper <- plogis(pcont$fit + qnorm(0.975) * pcont$se.fit)

# Plots
# plot
plist <- vector("list", length(predictors_cont))
npred <- length(predictors_cont)
names(plist) <- predictors_cont
for(p in predictors_cont) {
  # p = "ndvi"
  ddd <- pdata[pdata$varying_name == p, ]

  plist[[p]] <-
    ggplot(ddd, aes(varying_value, p_mle,
                    ymin = p_lower, ymax = p_upper,
                    color = fwi, fill = fwi, group = fwi)) +
    geom_ribbon(color = NA, alpha = 0.4) +
    geom_line() +
    facet_grid(rows = vars(vegetation_class), cols = vars(predictor)
               #,scales = "free_y") +
               ) +
    scale_fill_viridis(end = 0.8, option = "B") +
    scale_color_viridis(end = 0.8, option = "B") +
    theme(legend.title = element_text(),
          strip.background.y = element_blank(),
          strip.text.y = element_blank(),
          strip.text.x = element_text(size = 8),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_text(size = 6),
          legend.position = "none",
          plot.margin = unit(c(1, 1, 1, 1), "mm")) +
    ylim(0, 1)
  plist[[p]]
}

plist[["ndvi"]] <- plist[["ndvi"]] +
  theme(axis.title.y = element_text()) +
  ylab("Burned class probability")
plist[[npred]] <- plist[[npred]] +
  theme(strip.background.y = element_rect(),
        strip.text.y = element_text(angle = 270, size = 8))

leg <- g_legend(plist[[1]] + theme(legend.position = "bottom",
                                   legend.key.size = unit(4, 'mm'),
                                   legend.spacing.x = unit(2, 'mm'),
                                   legend.text = element_text(size = 8)))

p1 <- egg::ggarrange(plots = plist, ncol = 6)
p2 <- grid.arrange(p1, leg, nrow = 2, heights = c(20, 1))
ggsave("figures/climate-fuel-interactions_m4_y_0-1.png", plot = p2,
       width = 22, height = 15, units = "cm")


# Same model, using log_area as climate predictor -------------------------
# (total failure)

data$area_log <- log(data$area_ha)
d$area_ha %>% unique %>% log %>% hist

# set cluster
cl <- makeCluster(10, type = "FORK")

# Fit model
m4 <- bam(
  burned ~
    vegetation_class + vegetation_class : area_log +
    s(ndvi, by = vegetation_class, k = 5, bs = "cr") +
      ti(ndvi, area_log, k = 3) +
    s(elevation, by = vegetation_class, k = 5, bs = "cr") +
      ti(elevation, area_log, k = 3) +
    s(slope, by = vegetation_class, k = 5, bs = "cr") +
      ti(slope, area_log, k = 3) +
    s(aspect, by = vegetation_class, k = 5, bs = "cc") +
      ti(aspect, area_log, k = 4, bs = c("cc", "cr")) +
    s(dist_humans, by = vegetation_class, k = 5, bs = "cr") +
      ti(dist_humans, area_log, k = 3, bs = "cr") +
    s(dist_roads, by = vegetation_class, k = 5, bs = "cr") +
      ti(dist_roads, area_log, k = 3, bs = "cr"),
  knots = list(aspect = c(0, 360)),
  data = data, family = binomial(link = "logit"), method = "REML",
  cluster = cl
)
# algorithm did not converge

m4 <- brm(
  burned ~
    vegetation_class + vegetation_class : area_log +
    s(ndvi, by = vegetation_class, k = 5, bs = "cr") +
      # t2(ndvi, area_log, k = 3) +
    s(elevation, by = vegetation_class, k = 5, bs = "cr") +
      # t2(elevation, area_log, k = 3) +
    s(slope, by = vegetation_class, k = 5, bs = "cr") +
      # t2(slope, area_log, k = 3) +
    s(aspect, by = vegetation_class, k = 5, bs = "cc") +
      # t2(aspect, area_log, k = 4, bs = c("cc", "cr")) +
    s(dist_humans, by = vegetation_class, k = 5, bs = "cr") +
      # t2(dist_humans, area_log, k = 3, bs = "cr") +
    s(dist_roads, by = vegetation_class, k = 5, bs = "cr"), #+
      # t2(dist_roads, area_log, k = 3, bs = "cr"),
  knots = list(aspect = c(0, 360)),
  data = data, family = bernoulli(link = "logit"),
  cores = 1,
  chains = 1,
  iter = 10
)
# Deja usar s(), pero no t2... dice
# Error in X %*% diag(diagU[indi]) : non-conformable arguments

m4 <- brm(
  burned ~
    vegetation_class + vegetation_class : area_log +
    # without 3-way interactions including veg type
    t2(ndvi, area_log, k = 4) +
    t2(elevation, area_log, k = 4) +
    t2(slope, area_log, k = 4) +
    t2(aspect, area_log, k = 4, bs = c("cc", "cr")) +
    t2(dist_humans, area_log, k = 4, bs = "cr") +
    t2(dist_roads, area_log, k = 4, bs = "cr"),
  knots = list(aspect = c(0, 360)),
  data = data, family = bernoulli(link = "logit"),
  cores = 6,
  chains = 6,
  iter = 3000,
  warmup = 1000
)
# It took so long that I killed it after a few hours not reaching 10 %



# Without veg-type interactions -------------------------------------------

m5 <- bam(
  burned ~
    vegetation_class + vegetation_class : fwi +
    s(ndvi, k = 5, bs = "cr") +
      ti(ndvi, fwi, k = 3) +
    s(elevation, k = 5, bs = "cr") +
      ti(elevation, fwi, k = 3) +
    s(slope, k = 5, bs = "cr") +
      ti(slope, fwi, k = 3) +
    s(aspect, k = 5, bs = "cc") +
      ti(aspect, fwi, k = 4, bs = c("cc", "cr")) +
    s(dist_humans, k = 5, bs = "cr") +
      ti(dist_humans, fwi, k = 3, bs = "cr") +
    s(dist_roads, k = 5, bs = "cr") +
      ti(dist_roads, fwi, k = 3, bs = "cr"),
  knots = list(aspect = c(0, 360)),
  data = data, family = binomial(link = "logit"), method = "REML"
) # sooooo fast
# saveRDS(m5, "burn_probability_model_bam_m5.rds")
# m5 <- readRDS("burn_probability_model_bam_m5.rds")

# get HDI and median by predictor and by vegetation types
predictors_cont <- c("ndvi",  # cont for continuous
                     "elevation", "slope", "aspect",
                     "dist_humans", "dist_roads")

# use 100000 points from d (marginal balanced) to get summaries
dsub <- d[sample(1:nrow(d), size = 1e5, replace = F), ]

summaries <- do.call("rbind", lapply(predictors_cont, function(p) {
  mm <- quantiles(dsub[, p], ci = 0.95)
  mm2 <- matrix(mm, nrow = 1)
  colnames(mm2) <- names(mm)
  mm2 <- as.data.frame(mm2)
  mm2$predictor <- p
  return(mm2)
}))

# correct limits for aspect
summaries$mu_lower[summaries$predictor == "aspect"] <- 0
summaries$mu_upper[summaries$predictor == "aspect"] <- 360

# get fwi summary
summary_fwi <- quantiles(d$fwi, ci = 0.95, name = "fwi")

# prediction data for continuous predictors
pdata <- do.call("rbind", lapply(predictors_cont, function(p) {
  # p = "ndvi"
  values_list <- lapply(predictors_cont, function(pred) {
    # pred = "slope"
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
  gg$vegetation_class <- "Subalpine\nforest"

  gg$varying_name <- p
  gg$varying_value <- as.numeric(gg[, p])

  return(gg)
}))

# tidy names
pdata$predictor <- factor(pdata$varying_name,
                          levels = predictors_cont,
                          labels = c(
                            "NDVI", "Elevation", "Slope",
                            "Aspect", "Dist. sett.",
                            "Dist. roads"
                          ))
pdata$vegetation_class <- factor(pdata$vegetation_class,
                                 levels = veg_labels)

# compute predictions
pcont <- predict(m5, pdata, se.fit = T)

# probs
pdata$p_mle <- plogis(pcont$fit)
pdata$p_lower <- plogis(pcont$fit + qnorm(0.025) * pcont$se.fit)
pdata$p_upper <- plogis(pcont$fit + qnorm(0.975) * pcont$se.fit)

# Plot
ggplot(pdata, aes(varying_value, p_mle,
                  ymin = p_lower, ymax = p_upper,
                  color = fwi, fill = fwi, group = fwi)) +
  geom_ribbon(color = NA, alpha = 0.4) +
  geom_line() +
  facet_wrap(vars(predictor), scales = "free_x") +
  scale_fill_viridis(end = 0.8, option = "B") +
  scale_color_viridis(end = 0.8, option = "B") +
  theme(plot.margin = unit(c(1, 1, 1, 1), "mm")) +
  ylim(0, 1)


plist <- vector("list", length(predictors_cont))
npred <- length(predictors_cont)
names(plist) <- predictors_cont
for(p in predictors_cont) {
  # p = "ndvi"
  ddd <- pdata[pdata$varying_name == p, ]

  plist[[p]] <-
    ggplot(ddd, aes(varying_value, p_mle,
                    ymin = p_lower, ymax = p_upper,
                    color = fwi, fill = fwi, group = fwi)) +
    geom_ribbon(color = NA, alpha = 0.4) +
    geom_line() +
    facet_grid(rows = vars(vegetation_class), cols = vars(predictor)
               #,scales = "free_y") +
    ) +
    scale_fill_viridis(end = 0.8, option = "B") +
    scale_color_viridis(end = 0.8, option = "B") +
    theme(legend.title = element_text(),
          strip.background.y = element_blank(),
          strip.text.y = element_blank(),
          strip.text.x = element_text(size = 8),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_text(size = 6),
          legend.position = "none",
          plot.margin = unit(c(1, 1, 1, 1), "mm")) +
    ylim(0, 1)
  plist[[p]]
}

plist[["ndvi"]] <- plist[["ndvi"]] +
  theme(axis.title.y = element_text()) +
  ylab("Burned class probability")
plist[[npred]] <- plist[[npred]] +
  theme(strip.background.y = element_rect(),
        strip.text.y = element_text(angle = 270, size = 8))

leg <- g_legend(plist[[1]] + theme(legend.position = "bottom",
                                   legend.key.size = unit(4, 'mm'),
                                   legend.spacing.x = unit(2, 'mm'),
                                   legend.text = element_text(size = 8)))

p1 <- egg::ggarrange(plots = plist, ncol = 6)
p2 <- grid.arrange(p1, leg, nrow = 2, heights = c(20, 1))
ggsave("figures/climate-fuel-interactions_m5_y_0-1.png", plot = p2,
       width = 22, height = 15, units = "cm")


# Quizás el problema está en que no estoy mirando la diferencia -----------

# Tendré que parear píxeles para ver si así encuentro un efecto mayor.
# Y quizás mirando la diferencia se note más? No lo sé.

# Sí, tiene todo el sentido! Imaginate que siempre los quemados tengan menor
# NDVI, y te tocan muchos casos del tipo
#  burned, not
#  [10, 20]
#  [100, 110]
#  [50, 60]

# En promedio, el efecto del NDVI va a ser muy pequeño, porque varía en un rango
# re amplio (a escala grande), pero a escala pequeña, mirando pares,
# siempre es menor en el quemado. Entonces, haciendo eso sí podemos ver efectos.

# El problema de las diferencias es que no deja que haya eff cuadráticos
# ni que haya predictoras categóricas.
# Aunque el eff cuadrático podría mirarse poniendo como predictoras a
# diff(x) + diff(x) : mean(x)

# Quizás simplemente este método no permita evaluar el efecto del tipo de
# vegetación.