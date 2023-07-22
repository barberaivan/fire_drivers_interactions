options(scipen = 999)

# Packages ----------------------------------------------------------------

library(tidyverse); theme_set(theme_bw())
library(viridis)

library(grid)
library(egg)      # has its own ggarrange! much better than ggpubr
library(ggh4x)    # varying strip theme for veg_types and all together

library(mgcv)     # bam()
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


# Prepare data ------------------------------------------------------------

# for now, just use the large raw dataset
d <- read.csv("data_dynamic_burn_prob_more_ones.csv")
nrow(d)
d_year <- aggregate(burned ~ year, d, mean)
d_year_sum <- aggregate(burned ~ year, d, sum)
# no hay nada de burned. Hay que hacer muestreo estratificado y trabajar con
# datos desproporcionados en realción a lo real.
sum(d_year_sum$burned) / nrow(d) # 0.002670833 es la annual fire prob

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

# # make less-imbalanced dataset
# N_unburned <- 3000
# id_unburned <- which(d$burned == 0)
# ids_keep_unburned <- sample(id_unburned, size = N_unburned, replace = F)
# data <- d[c(ids_keep_unburned, which(d$burned == 1)), ]
data <- d[d$vegetation_type > 1, ]
nrow(data)

burn_prob_sample <- mean(data$burned)
aggregate(burned ~ vegetation_code, data, length)


# Simple model ------------------------------------------------------------

# We fit a logistic regression model oversampling the ones (8000 pixels),
# with less unburned pixels than the real proportion.
# (The population proportion is 0.002670833)
burn_prob_pop <- 0.002670833


m1 <- bam(
  burned ~
    vegetation_class * fwi +
    s(ndvi, by = vegetation_class, k = 7, bs = "cr") +
      ti(ndvi, fwi, k = 4) +
    s(elevation, by = vegetation_class, k = 7, bs = "cr") +
      ti(elevation, fwi, k = 4) +
    s(slope, by = vegetation_class, k = 7, bs = "cr") +
      ti(slope, fwi, k = 4) +
    s(aspect, by = vegetation_class, k = 7, bs = "cc") +
      ti(aspect, fwi, k = 4, bs = c("cc", "cr")) +
    s(dist_humans, by = vegetation_class, k = 7, bs = "cr") +
      ti(dist_humans, fwi, k = 4, bs = "cr") +
    s(dist_roads, by = vegetation_class, k = 7, bs = "cr") +
      ti(dist_roads, fwi, k = 4, bs = "cr"),
  knots = list(aspect = c(0, 360)),
  data = data, family = binomial(link = "logit"), method = "REML"
)

saveRDS(m1, "burn_probability_model_bam_m1.rds")
summary(m1)
# plot(m1)
# gam.check(m1)

# res_m1 <- simulateResiduals(m1, n = 3000, integerResponse = TRUE)
# # plot(res_m1)
# plotResiduals(res_m1, form = data$ndvi, rank = F)
# plotResiduals(res_m1, form = data$elevation, rank = F)
# plotResiduals(res_m1, form = data$slope, rank = F)
# plotResiduals(res_m1, form = data$aspect, rank = F)
# plotResiduals(res_m1, form = data$fwi, rank = F)
# plotResiduals(res_m1, form = as.factor(data$year))
# plotResiduals(res_m1, form = data$vegetation_class, rank = F)
# plotResiduals(res_m1, form = data$dist_humans, rank = F)
# plotResiduals(res_m1, form = data$dist_roads, rank = F)
# # Anda hermoso. Plotear predicciones lueguito.



# Predictions -------------------------------------------------------------

# get HDI and median by predictor and by vegetation types
predictors_cont <- c("ndvi",  # cont for continuous
                     "elevation", "slope", "aspect",
                     "dist_humans", "dist_roads")

# use 100000 points from d (marginal balanced) to get summaries
dsub <- d[sample(1:nrow(d), size = 1e5, replace = F), ]

summaries <- do.call("rbind", lapply(predictors_cont, function(p) {
  agg <- aggregate(dsub[, p] ~ vegetation_class, dsub, hdmedian)
  agg <- cbind(agg$vegetation_class, as.data.frame(agg$`dsub[, p]`))
  colnames(agg)[1] <- "vegetation_class"
  agg$predictor <- p
  return(agg)
}))

# correct limits for aspect
summaries$mu_lower[summaries$predictor == "aspect"] <- 0
summaries$mu_upper[summaries$predictor == "aspect"] <- 360

# get fwi summary
summary_fwi <- hdmedian(d$fwi, name = "fwi")



# Effects of static predictors --------------------------------------------

# For every vegetation type separately,
# vary every predictor in its 0.8 HDI, while keeping the others in their medians,
# and varying fwi at 3 values (summary_fwi).

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
pcont <- predict(m1, pdata, se.fit = T)
pcont$fit_correct <- population_intercept(pcont$fit, burn_prob_pop, burn_prob_sample)

# probs
pdata$p_mle <- plogis(pcont$fit_correct)
pdata$p_lower <- plogis(pcont$fit_correct + qnorm(0.025) * pcont$se.fit)
pdata$p_upper <- plogis(pcont$fit_correct + qnorm(0.975) * pcont$se.fit)

# probs balanced
pdata$pbal_mle <- plogis(pcont$fit)
pdata$pbal_lower <- plogis(pcont$fit + qnorm(0.025) * pcont$se.fit)
pdata$pbal_upper <- plogis(pcont$fit + qnorm(0.975) * pcont$se.fit)

# odds
pdata$odds_mle <- exp(pcont$fit_correct)
pdata$odds_lower <- exp(pcont$fit_correct + qnorm(0.025) * pcont$se.fit)
pdata$odds_upper <- exp(pcont$fit_correct + qnorm(0.975) * pcont$se.fit)

# linpred
pdata$linpred_mle <- pcont$fit
pdata$linpred_lower <- pcont$fit + qnorm(0.025) * pcont$se.fit
pdata$linpred_upper <- pcont$fit + qnorm(0.975) * pcont$se.fit

# plot

# # probability
# ggplot(pdata, aes(varying_value, p_mle, ymin = p_lower, ymax = p_upper,
#                   color = fwi, fill = fwi, group = fwi)) +
#   geom_ribbon(color = NA, alpha = 0.4) +
#   geom_line() +
#   facet_grid(rows = vars(vegetation_class), cols = vars(predictor),
#              scales = "free") +
#   scale_fill_viridis(end = 0.8, option = "B") +
#   scale_color_viridis(end = 0.8, option = "B") +
#   theme(legend.title = element_text())
#
# # probability balanced
# ggplot(pdata, aes(varying_value, pbal_mle,
#                   ymin = pbal_lower, ymax = pbal_upper,
#                   color = fwi, fill = fwi, group = fwi)) +
#   geom_ribbon(color = NA, alpha = 0.4) +
#   geom_line() +
#   facet_grid(rows = vars(vegetation_class), cols = vars(predictor),
#              scales = "free") +
#   scale_fill_viridis(end = 0.8, option = "B") +
#   scale_color_viridis(end = 0.8, option = "B") +
#   theme(legend.title = element_text())
#
# # odds
# ggplot(pdata, aes(varying_value, odds_mle,
#                   ymin = odds_lower, ymax = odds_upper,
#                   color = fwi, fill = fwi, group = fwi)) +
#   geom_ribbon(color = NA, alpha = 0.4) +
#   geom_line() +
#   facet_grid(rows = vars(vegetation_class), cols = vars(predictor),
#              scales = "free") +
#   scale_fill_viridis(end = 0.8, option = "B") +
#   scale_color_viridis(end = 0.8, option = "B") +
#   theme(legend.title = element_text())
#
# # linear predictor
# ggplot(pdata, aes(varying_value, linpred_mle,
#                   ymin = linpred_lower, ymax = linpred_upper,
#                   color = fwi, fill = fwi, group = fwi)) +
#   geom_ribbon(color = NA, alpha = 0.4) +
#   geom_line() +
#   facet_grid(rows = vars(vegetation_class), cols = vars(predictor),
#              scales = "free") +
#   scale_fill_viridis(end = 0.8, option = "B") +
#   scale_color_viridis(end = 0.8, option = "B") +
#   theme(legend.title = element_text())


# Better plot

plist <- vector("list", length(predictors_cont))
npred <- length(predictors_cont)
names(plist) <- predictors_cont
for(p in predictors_cont) {
  # p = "ndvi"
  ddd <- pdata[pdata$varying_name == p, ]

  plist[[p]] <-
    ggplot(ddd, aes(varying_value, p_mle, ymin = p_lower, ymax = p_upper,
                    color = fwi, fill = fwi, group = fwi)) +
    geom_ribbon(color = NA, alpha = 0.4) +
    geom_line() +
    facet_grid(rows = vars(vegetation_class), cols = vars(predictor),
               scales = "free") +
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
          plot.margin = unit(c(1, 1, 1, 1), "mm"))
  plist[[p]]
}

plist[["ndvi"]] <- plist[["ndvi"]] +
                    theme(axis.title.y = element_text()) +
                    ylab("Annual burn probability")
plist[[npred]] <- plist[[npred]] +
                  theme(strip.background.y = element_rect(),
                        strip.text.y = element_text(angle = 270, size = 8))

leg <- g_legend(plist[[1]] + theme(legend.position = "bottom",
                                   legend.key.size = unit(4, 'mm'),
                                   legend.spacing.x = unit(2, 'mm'),
                                   legend.text = element_text(size = 8)))

p1 <- egg::ggarrange(plots = plist,#plist[[1]], plist[[2]], plist[[3]], plist[[4]],
                     #plist[[5]], plist[[6]],
                     ncol = 6)
p2 <- grid.arrange(p1, leg, nrow = 2, heights = c(20, 1))
ggsave("figures/climate-fuel-interactions_m1.png", plot = p2,
       width = 22, height = 15, units = "cm")





# NOTAS -------------------------------------------------------------------

# Dará lo mismo hacer el análisis a escala de incendio, fijando la burnprob en
# 0.5? Quizás el efecto sería similar a usar una proporción balanceada... no sé

# Tendríamos menos ruido del paisaje lejano que no se quemó porque no hubo
# igniciones.

# Luego rehacer el plot asumiendo que la prob poblacional es 0.5

pcont$fit_half <- population_intercept(pcont$fit, 0.5, burn_prob_sample)

# probs balanced exactly (0.5)
pdata$phalf_mle <- plogis(pcont$fit_half)
pdata$phalf_lower <- plogis(pcont$fit_half + qnorm(0.025) * pcont$se.fit)
pdata$phalf_upper <- plogis(pcont$fit_half + qnorm(0.975) * pcont$se.fit)

# plot

# Better plot

plist_half <- vector("list", length(predictors_cont))
npred <- length(predictors_cont)
names(plist_half) <- predictors_cont
for(p in predictors_cont) {
  # p = "ndvi"
  ddd <- pdata[pdata$varying_name == p, ]

  plist_half[[p]] <-
    ggplot(ddd, aes(varying_value, phalf_mle,
                    ymin = phalf_lower, ymax = phalf_upper,
                    color = fwi, fill = fwi, group = fwi)) +
    geom_ribbon(color = NA, alpha = 0.4) +
    geom_line() +
    facet_grid(rows = vars(vegetation_class), cols = vars(predictor),
               scales = "free") +
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
          plot.margin = unit(c(1, 1, 1, 1), "mm"))
  plist_half[[p]]
}

plist_half[["ndvi"]] <- plist_half[["ndvi"]] +
  theme(axis.title.y = element_text()) +
  ylab("Annual burn probability")
plist_half[[npred]] <- plist_half[[npred]] +
  theme(strip.background.y = element_rect(),
        strip.text.y = element_text(angle = 270, size = 8))

leg <- g_legend(plist_half[[1]] + theme(legend.position = "bottom",
                                   legend.key.size = unit(4, 'mm'),
                                   legend.spacing.x = unit(2, 'mm'),
                                   legend.text = element_text(size = 8)))

p1 <- egg::ggarrange(plots = plist_half,#plist_half[[1]], plist_half[[2]], plist_half[[3]], plist_half[[4]],
                     #plist_half[[5]], plist_half[[6]],
                     ncol = 6)
p2 <- grid.arrange(p1, leg, nrow = 2, heights = c(20, 1))
ggsave("figures/climate-fuel-interactions_m1_prob-pop-05.png", plot = p2,
       width = 22, height = 15, units = "cm")

# Cambia un poco, pero no tantísimo. El problema es el efecto global del fwi


# IMPORTANTE: -------------------------------------------------------------

# fijé las interacciones (ti) across veg types. Intentar relajar este supuesto.
# Hacer el análisis basado en fuegos, tomando píxeles del interior y de un buffer
# de igual área.
# Además del modelo, cuantificar la proporción de píxeles no quemados que son
# no combustibles. Aunque quizás para eso haga falta ver un borde más chico
# (100 m?)

# Para tener un set balanceado habría que primero tirar pixeles regularmente.
# Luego, ver el N adentro y afuera. Luego, elegir el N min.

# Ahora, si el modelo incluye la categoría non-burnable (que no tendría sentido),
# los datos no van a estar balanceados. Y el modelo es para ver los efectos sobre
# lo no-quemable. O sea, necesitamos los datos balanceados. Podemos re-balancear
# a posteriori.

# Y para ver la prop de no quemadble en el borde, es mejor calcularlo por área.

# Leer paper de Holsinger sobre qué frena a los fuegos.

# Podemos usar el FWI como predictora organizadora o el fire size. Digamos que
# el clima es mejor, pero el fire size puede variar a escala diaria.

