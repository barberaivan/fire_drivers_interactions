# in this code I use a smaller dataset to allow the exploration of
# residual spatial correlation

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

library(terra)    # import and write raster
library(sptotal)  # spatial correlation plot

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

v <- vect("data/data_dynamic_burn_prob_more_ones_smaller_dataset.shp")
vflat <- project(v, "EPSG:5343")
vcoord <- crds(vflat) %>% as.data.frame
range_q <- diff(range(vcoord$y)) / diff(range(vcoord$x))

# for now, just use the large raw dataset
d <- cbind(as.data.frame(v), vcoord)
# remove unburnable
d <- d[d$vegetation > 1, ]
nrow(d) # 11624

d_year <- aggregate(burned ~ year, d, mean)
d_year_sum <- aggregate(burned ~ year, d, sum)
sum(d_year_sum$burned) / nrow(d)
# 0.169 es la prop de pixeles quemados en mi dataset
# 0.002670833 es la annual fire prob real

# recode vegetation and remove unuseful categories
burned_veg <- read.csv("data/data_burned_and_available_area_by_vegetation_dryforest2.csv")
d$vegetation_code <- d$vegetation

# code dry forest A as subalpine forest
d$vegetation_code[d$vegetation_code == 4] <- 2
# plantation as dry forest B
d$vegetation_code[d$vegetation_code == 9] <- 5
# anthrop as shrubland
d$vegetation_code[d$vegetation_code == 8] <- 6

d <- left_join(d, burned_veg[, c("vegetation_code", "vegetation_class")],
                  by = "vegetation_code")

d$vegetation_class[d$vegetation_class == "Dry forest B"] <- "Dry forest"

unique(d$vegetation_class)
veg_levels <- c(
  "Wet forest",
  "Subalpine forest",
  "Dry forest",
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

# remove the burned sampled in the unburned dataset
data <- d[!(d$dataset == "unburned" & d$burned == 1), ]
# remove duplicated coordinates
data <- data[!duplicated(data[, c("x", "y")]), ]
# rename var
names(data)[names(data) == "dist_human"] <- "dist_humans"
# add northing
data$northing <- cos(data$aspect * (pi / 180))


burn_prob_sample <- mean(data$burned)
# aggregate(burned ~ vegetation_code, data, length)

# numeric response

# covariates and metadata names
cov_all <- c(
  "vegetation", "ndvi",
  # "unburnable", # not used
  # "pp", "temp",

  "elevation", "slope", "northing",
  "TPI300", "TPI1k", "TPI2k",
  "HLI_south", "ERR450", "Rough450",

  "dist_humans", "dist_roads"
)

# without vegetation class
cov_cont <- cov_all[-1]

meta_names <- c("fire_id", "year", "area_ha", "fwi")

# Variable selection based on univariate models (R2) ----------------------

ncont <- length(cov_cont)
mlist_uni <- vector("list", ncont)
names(mlist_uni) <- cov_cont
r2_uni <- data.frame(var = cov_cont, r2 = NA)

for(i in 1:ncont) {
  print(i)
  # i = 11
  dtemp <- data[, c("burned", cov_cont[i])]
  names(dtemp) <- c("y", "x")
  modelo <- bam(y ~ s(x, bs = "cr", k = 4), data = dtemp,
                family = "binomial")
  mlist_uni[[i]] <- modelo
  p <- predict(modelo, type = "response")
  r2_uni[i, "r2"] <- var(p) / (var(p) + mean(p * (1 - p)))

  xx <- quantile(dtemp$x, probs = c(0.05, 0.95))
  # xseq <- seq(min(dtemp$x), max(dtemp$x), length.out = 150)
  xseq <- seq(xx[1], xx[2], length.out = 150)
  nd <- data.frame(x = xseq, pair = rep(1, 150))
  pp <- predict(modelo, nd, type = "response")
  plot(pp ~ xseq, main = cov_cont[i], type = "l")
}

cm <- cor(data[, cov_cont])
# View(round(cm, 4))
# View(round(abs(cm), 4))

r2_uni[order(r2_uni$r2, decreasing = T), ]
candidates <- c(
  "elevation",
  "ndvi",
  "northing",
  "Rough450",
  "TPI2k",
  "slope",
  "dist_roads",
  "dist_humans",
  "ERR450"
)

# View(round(abs(cm[candidates, candidates]), 4))
# high correlations:
# TPI2k, ERR
# slope, rough
# dist_roads, dist_human

# coeffs correlation? removing only slope
model_corr_1 <- bam(
  burned ~
    vegetation_class +
    s(elevation, bs = "cr", k = 5) +
    s(ndvi, bs = "cr", k = 5) +
    s(northing, bs = "cr", k = 5) +
    s(Rough450, bs = "cr", k = 5) +
    # s(slope, bs = "cr", k = 5) +
    s(TPI2k, bs = "cr", k = 5) +
    s(ERR450, bs = "cr", k = 5) +
    s(dist_roads, bs = "cr", k = 5) +
    s(dist_humans, bs = "cr", k = 5),
  data = data, family = "binomial"
)
concurvity(model_corr_1)
# slope increases the concurvity. use rough

# FWI model ---------------------------------------------------------------

# We fit a logistic regression model oversampling the ones (8000 pixels),
# with less unburned pixels than the real proportion.
# (The population proportion is 0.002670833)
burn_prob_pop <- 0.002670833

m1 <- bam(
  burned ~
    vegetation_class * fwi +
    s(ndvi, k = 4, by = vegetation_class, bs = "cr") +
      ti(ndvi, fwi, k = 4) +
    s(elevation, k = 4, bs = "cr") +
      ti(elevation, fwi, k = 4) +
    s(northing, k = 4, bs = "cr") +
      ti(northing, fwi, k = 4) +
    s(Rough450, k = 4, bs = "cc") +
      ti(Rough450, fwi, k = 4, bs = c("cc", "cr")) +
    # s(TPI2k, k = 4, bs = "cc") +
    #   ti(TPI2k, fwi, k = 4, bs = c("cc", "cr")) +
    # s(ERR450, k = 4, bs = "cc") +
    #   ti(ERR450, fwi, k = 4, bs = c("cc", "cr")) +
    s(dist_humans, k = 4, bs = "cr") +
      ti(dist_humans, fwi, k = 4, bs = "cr") +
    s(dist_roads, k = 4, bs = "cr") +
      ti(dist_roads, fwi, k = 4, bs = "cr"),
  data = data, family = binomial(link = "logit"), method = "REML"
)
# saveRDS(m1, "burn_probability_model_bam_m1.rds")
# summary(m1)
# plot(m1)
# gam.check(m1)

res_m1 <- simulateResiduals(m1, n = 1000, integerResponse = TRUE)
hist(res_m1$scaledResiduals, breaks = 30)
plot(res_m1)
plotResiduals(res_m1, form = data$ndvi, rank = F)
plotResiduals(res_m1, form = data$elevation, rank = F)
plotResiduals(res_m1, form = data$slope, rank = F)
plotResiduals(res_m1, form = data$aspect, rank = F)
plotResiduals(res_m1, form = data$fwi, rank = F)
plotResiduals(res_m1, form = as.factor(data$year))
plotResiduals(res_m1, form = data$vegetation_class, rank = F)
plotResiduals(res_m1, form = data$dist_humans, rank = F)
plotResiduals(res_m1, form = data$dist_roads, rank = F)
# Perfecto

# correlación espacial?
m1_spat_corr <- testSpatialAutocorrelation(
  res_m1, x = data$x, y = data$y, plot = T
)
# sí, hay mucha


m2 <- bam(
  burned ~
    vegetation_class * fwi +
    s(ndvi, k = 4, by = vegetation_class, bs = "cr") +
    ti(ndvi, fwi, k = 4) +
    s(elevation, k = 4, bs = "cr") +
    ti(elevation, fwi, k = 4) +
    s(northing, k = 4, bs = "cr") +
    ti(northing, fwi, k = 4) +
    s(Rough450, k = 4, bs = "cc") +
    ti(Rough450, fwi, k = 4, bs = c("cc", "cr")) +
    # s(TPI2k, k = 4, bs = "cc") +
    #   ti(TPI2k, fwi, k = 4, bs = c("cc", "cr")) +
    # s(ERR450, k = 4, bs = "cc") +
    #   ti(ERR450, fwi, k = 4, bs = c("cc", "cr")) +
    s(dist_humans, k = 4, bs = "cr") +
    ti(dist_humans, fwi, k = 4, bs = "cr") +
    s(dist_roads, k = 4, bs = "cr") +
    ti(dist_roads, fwi, k = 4, bs = "cr") +

    # spatial effect (to account for spatial correlation)
    te(x, y, bs = "cr", k = c(10, 4)),
  data = data, family = binomial(link = "logit"), method = "REML"
)

res_m2 <- simulateResiduals(m2, n = 1000, integerResponse = TRUE)
(m2_spat_corr <- testSpatialAutocorrelation(
  res_m2, x = data$x, y = data$y, plot = T
))
# Sigue habiendo corr



# visual inspection of residual correlation.
data$res_m1 <- res_m1$scaledResiduals
data$res_m2 <- res_m2$scaledResiduals

sv1 <- sv(data, which(names(data) == "x"), which(names(data) == "y"),
          which(names(data) == "res_m1"), bins = 100, cutoff = 20000)
gc()


sv2 <- sv(data, which(names(data) == "x"), which(names(data) == "y"),
          which(names(data) == "res_m2"), bins = 100, cutoff = 20000)
gc()

par(mfrow = c(1, 2))
plot(gamma ~ dist, sv1, main = "model 1", ylim = c(0, 0.1))
plot(gamma ~ dist, sv2, main = "model 2", ylim = c(0, 0.1))
par(mfrow = c(1, 1))

# la corr espacial de los residuos es re baja; usar m1

# Predictions -------------------------------------------------------------

# get HDI and median by predictor and by vegetation types
predictors_cont <- c(
  "ndvi",
  "elevation",
  "northing",
  "Rough450",
  # "TPI2k",
  # "ERR450",
  "dist_roads",
  "dist_humans"
)

summaries <- do.call("rbind", lapply(predictors_cont, function(p) {
  print(p)
  agg <- aggregate(data[, p] ~ vegetation_class, data, hdmedian)
  agg <- cbind(agg$vegetation_class, as.data.frame(agg$`data[, p]`))
  colnames(agg)[1] <- "vegetation_class"
  agg$predictor <- p
  return(agg)
}))

# get fwi summary
summary_fwi <- hdmedian(d$fwi, name = "fwi", ci = 0.95)

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
                   length.out = 150)
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
                            "NDVI", "Elevation", "Northing",
                            "Rough450", #"TPI 2k", "ERR450",
                            "Dist. roads", "Dist. sett."
                          ))

pdata$vegetation_class <- factor(pdata$vegetation_class,
                          levels = veg_labels)

# add coordinates (in the mid)
pdata$x <- 1540000
pdata$y <- 5540000

# compute predictions
pcont <- predict(m1, pdata, se.fit = T)
pcont$fit_correct <- population_intercept(pcont$fit, burn_prob_pop, burn_prob_sample)

# probs
pdata$p_mle <- plogis(pcont$fit_correct)
pdata$p_lower <- plogis(pcont$fit_correct + qnorm(0.025) * pcont$se.fit)
pdata$p_upper <- plogis(pcont$fit_correct + qnorm(0.975) * pcont$se.fit)


# plot

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

p1 <- egg::ggarrange(plots = plist, ncol = npred)
p2 <- grid.arrange(p1, leg, nrow = 2, heights = c(20, 1))
ggsave("figures/climate-fuel-interactions_m1_img3_smaller_data.png", plot = p2,
       width = 24, height = 15, units = "cm")
ggsave("figures/climate-fuel-interactions_m2_img3_smaller_data.png", plot = p2,
       width = 24, height = 15, units = "cm")


# adding a smooth for coordinates yields high posterior uncertainty, but the
# patterns remain the same.


# Map burn probabilities ----------------------------------------------------
# Using three fwi values: summary_fwi

img <- rast("predictors_image_brc-cholila.tif")
img_vals <- values(img) %>% as.data.frame()

# recode vegetation
summary(img_vals)

img_vals$vegetation_code <- img_vals$vegetation_valdivian

# code dry forest A as subalpine forest
img_vals$vegetation_code[img_vals$vegetation_code == 4] <- 2
# plantation as dry forest B
img_vals$vegetation_code[img_vals$vegetation_code == 9] <- 5
# anthrop as shrubland
img_vals$vegetation_code[img_vals$vegetation_code == 8] <- 6

img_vals <- left_join(img_vals, burned_veg[, c("vegetation_code", "vegetation_class")],
                      by = "vegetation_code")

unique(img_vals$vegetation_class)
img_vals$vegetation_class <- factor(img_vals$vegetation_class,
                                    levels = veg_levels,
                                    labels = veg_labels)

img_vals$northing <- cos(img_vals$aspect * (pi / 180))
head(img_vals)
# design matrix
p_mat <- data.frame(low = rep(NA, nrow(img_vals)),
                    mid = rep(NA, nrow(img_vals)),
                    high = rep(NA, nrow(img_vals)))

not_na <- complete.cases(img_vals)
for(i in 1:ncol(p_mat)) {
  img_vals$fwi <- summary_fwi[i]
  lp_local <- predict(m1, img_vals[not_na, ], type = "link")
  p_mat[not_na, i] <- plogis(population_intercept(lp_local,
                                                  burn_prob_pop,
                                                  burn_prob_sample))
  rm(lp_local)
  gc()
}


# save images with fwi values
img_fwi <- img[[1:3]]
names(img_fwi) <- c("fwi_low", "fwi_mid", "fwi_high")
values(img_fwi) <- p_mat
head(img_fwi)

for(i in 1:3) {
  writeRaster(img_fwi[[i]], paste("predicted_burn_prob_m1_", names(img_fwi)[i],
                                  ".tif", sep = ""))
}



# yearly aggregates
dyear <- aggregate(fwi ~ year, data, mean)
dyear <- dyear[order(dyear$fwi, decreasing = TRUE), ]
summary_fwi
# year         fwi
# 17 2015  2.13690876
# 10 2008  1.28572634
# 11 2009  1.07442468
# 21 2019  1.05802444
# 18 2016  0.92919911
# 23 2021  0.85766246
# 4  2002  0.84591191
# 16 2014  0.74815353
# 1  1999  0.55428037
# 22 2020  0.20878818
# 20 2018 -0.08350853
# 14 2012 -0.13634643
# 6  2004 -0.21088284
# 5  2003 -0.33217322
# 19 2017 -0.39841721
# 24 2022 -0.46953133
# 9  2007 -0.50957321
# 7  2005 -0.59761123
# 8  2006 -0.92682204
# 2  2000 -0.94269612
# 15 2013 -1.05697706
# 13 2011 -1.31273812
# 3  2001 -1.34690247
# 12 2010 -1.37261026


# bring fwi data
fwi_data <- read.csv("data_climate_interannual_fwi.csv")
fwi_data$date <- as.Date(fwi_data$date, format = "%Y-%m-%d")
fwi_data$year <- format(fwi_data$date, format = "%Y") %>% as.numeric
fwi_data$month <- format(fwi_data$date, format = "%m") %>% as.numeric
fwi_data$fseason <- NA
fwi_data$fseason <- fwi_data$year
fwi_data$fseason[fwi_data$month == 12] <- fwi_data$year[fwi_data$month == 12] + 1

fwi_agg <- aggregate(fwi ~ fseason, fwi_data[fwi_data$month %in% c(12, 1:3), ],
                     FUN = mean)
names(fwi_agg) <- c("year", "fwi")

# remaining climatic variables
climate_long <- read.csv("data_climate_interannual.csv")
climate <- pivot_wider(climate_long[, c("variable", "value", "year")],
                       names_from = "variable",
                       values_from = "value")

# merge with fwi
climate <- cbind(climate[climate$year > 1998, ], fwi = fwi_agg$fwi)



# merge with fire data
burned_annual_clim <- left_join(
  burned_annual, climate, by = "year"
)
