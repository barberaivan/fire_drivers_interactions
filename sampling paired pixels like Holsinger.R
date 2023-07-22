library(terra)
library(sf)
library(tidyverse)


# Functions ---------------------------------------------------------------

# Function to get buffer of fire with a given relative area, defined by ratio.
# (ratio = 0.5 returns an inner buffer with half the area of the fire).
# fireg is the sf geometry of a fire (usually, a multipolygon).
area_buffer <- function(fireg, ratio = 0.5, precision = 0.001, maxsteps = 100) {

  # # testing
  # fireg <- fireg
  # fireg <- fires[fires$area_ha == min(fires$area_ha), ]$geometry
  # plot(fireg)
  # ratio = 0.5
  # precision = 0.001
  # maxsteps = 100

  fire_area <- st_area(fireg)
  fire_length <- st_length(st_cast(fireg, "MULTILINESTRING"))

  width <- (ratio - 1) * fire_area / fire_length
  inner_poly <- st_buffer(fireg, width)
  achieved_precision <- as.numeric(st_area(inner_poly) / fire_area) - ratio

  steps <- 1

  while(abs(achieved_precision) > precision & steps < maxsteps) {
    # print(steps)
    inner_perimeter <- st_length(st_cast(inner_poly, "MULTILINESTRING"))
    width <- width - achieved_precision * fire_area / inner_perimeter
    inner_poly <- st_buffer(fireg, width)
    achieved_precision <- as.numeric(st_area(inner_poly) / fire_area) - ratio
    steps <- steps + 1
  }

  if(steps == maxsteps & abs(achieved_precision) > precision) {
    warning("Required precision not reached after ",
            maxsteps, " steps.\nAchieved precision = ", achieved_precision)
  }

  return(list(polygon = inner_poly,
              distance = as.numeric(width)))
}


# Data --------------------------------------------------------------------

# Import fires
fires_wgs <- st_read("../patagonian_fires/patagonian_fires.shp")

# Filter uncertain_year ones
fires_wgs$month <- format(as.Date(fires_wgs$date, "%Y-%m-%d"), "%m") %>% as.numeric
# table(fires_wgs$month) %>% plot
# Later we could filter winter months, but not now.
fires_wgs$obs[is.na(fires_wgs$obs)] <- "it_is_na"
fires_wgs <- fires_wgs[fires_wgs$obs != "uncertain_year", ]
# nrow(fires_wgs)

# Area issues:
# one pixel = 30 * 30 m = 900 m2 = 900 / 10000 = 0.09 ha
# 10 ha = 10 / 0.09 = 111.1111 pixels
fires_wgs$pixels <- fires_wgs$area_ha / (30 * 30 / 10000)
size_table <- data.frame(area = fires_wgs$area_ha,
                         pixels = fires_wgs$pixels,
                         fire_id = fires_wgs$fire_id)
size_table <- size_table[order(size_table$area), ]
size_table$row <- 1:nrow(size_table)
# View(size_table)

# For now, take 200 pixels by fire and class (edge-core)
n_points <- 200

# project to posgar 2007/1
fires <- st_transform(fires_wgs, crs = 5343)

# Get paired points in every fire polygon ----------------------------------

# (loop so if an error occurs, not everything will be lost)
points_data <- vector("list", length(fires$fire_id))

names_from_largest <- fires$fire_id[order(fires$area_ha, decreasing = TRUE)]

for(f in 1:length(fires$fire_id)) {

  # # test
  # fire_id <- "2016_87"#"2015_50"
  fire_id <- names_from_largest[f]
  fire <- fires[fires$fire_id == fire_id, ]
  fireg <- fire$geometry

  # get core area
  print(paste("fire ", fire_id, ": ", "making core area", sep = ""))
  core_list <- area_buffer(fireg)
  core_poly <- core_list$polygon
  buffer_distance <- abs(core_list$distance)

  # make many points along both polygons
  edge_length <- st_length(st_cast(fireg, "MULTILINESTRING")) %>% as.numeric
  edge_points <- st_segmentize(fireg,
                               dfMaxLength = edge_length / (n_points + 1)) %>%
    st_coordinates() %>%
    as.data.frame() %>%
    select(X, Y) %>%
    st_as_sf(coords = c("X", "Y"))

  core_length <- st_length(st_cast(core_poly, "MULTILINESTRING")) %>% as.numeric
  core_points_full <- st_segmentize(core_poly,
                                    dfMaxLength = core_length / (n_points + 1)) %>%
    st_coordinates() %>%
    as.data.frame() %>%
    select(X, Y) %>%
    st_as_sf(coords = c("X", "Y"))

  # keep only n_points for edge
  keep <- seq(1, nrow(edge_points), length.out = n_points) %>% round %>% unique
  edge_points <- edge_points[keep, ]

  # For every edge point, find the closest core point and match it
  edge_points$pair_id <- 1:nrow(edge_points)
  edge_points$point_id <- 1:nrow(edge_points)
  edge_points$distance_pair <- NA
  core_points <- edge_points # just to initialize spdf

  print(paste("fire ", fire_id, ": ", "pairing points", sep = ""))
  for(i in 1:n_points) {
    p <- edge_points[i, ]$geometry

    times <- 2
    ll <- 0
    while(ll < 1) {
      p_buffer <- st_buffer(p, buffer_distance * times)
      in_buffer <- st_intersects(p_buffer, core_points_full)[[1]]
      ll <- length(in_buffer)
      times <- times + 1
    }

    core_try <- core_points_full[in_buffer, ]
    distances <- st_distance(p, core_try) %>% as.numeric
    point_id <- which.min(distances)
    core_point <- core_try[point_id, ]
    core_point$pair_id <- i
    core_point$point_id <- point_id
    core_point$distance_pair <- distances[point_id]
    edge_points$distance_pair[i] <- distances[point_id]

    core_points[i, ] <- core_point
  }

  # # check:
  # par(mfrow = c(2, 1))
  # plot(edge_points[, "distance_pair"], col = "red", reset = F)
  # plot(core_points[, "distance_pair"], col = "blue", add = T)
  # plot(fireg, reset = F)
  # plot(core_poly, add = TRUE, col = "red")
  # par(mfrow = c(1, 1))

  # merge edge and core (long)
  core_points$class <- "core"
  edge_points$class <- "edge"

  points_both <- rbind(edge_points, core_points)
  # plot(points_both)

  points_both$fire_id <- fire_id
  points_both$year <- fire$year
  points_both$area_ha <- fire$area_ha
  points_both$date <- fire$date
  points_both$month <- fire$month
  points_both$bufdist <- buffer_distance

  points_data[[f]] <- points_both
}

# merge list
points_flat <- do.call("rbind", points_data)

# assign crs
st_crs(points_flat) <- st_crs(fires)

# reproject to latlong
points_wgs <- st_transform(points_flat, crs = st_crs(fires_wgs))

# write to disk
# st_write(points_wgs, "paired_points_edge-core.shp")




# Explore distance between paired pixels ----------------------------------
# as a function of fire area

pairs_vec <- vect("paired_points_edge-core.shp")
pairs <- as.data.frame(pairs_vec)
str(pairs)

pairs_agg <- aggregate(cbind(dstnc_p, bufdist, area_ha) ~ fire_id, pairs, median)
pairs_agg_q <- aggregate(cbind(dstnc_p, bufdist, area_ha) ~ fire_id, pairs,
                         quantile, probs = 0.1)

plot(dstnc_p ~ log(area_ha), pairs_agg)
plot(dstnc_p ~ log(area_ha), pairs_agg)
plot(bufdist ~ log(area_ha), pairs_agg[pairs_agg$area_ha >= 400, ])
plot(bufdist ~ area_ha, pairs_agg[pairs_agg$area_ha >= 400, ])

plot(dstnc_p ~ log(area_ha), pairs_agg_q)

nrow(pairs_agg[pairs_agg$dstnc_p > 90, ])
nrow(pairs_agg[pairs_agg$area_ha >= 400, ])


# Exploring FWI variation in large fires ----------------------------------

d <- read.csv("data_inside-outside.csv")
agg2 <- aggregate(fwi ~ fire_id, d, mean)

# merge with distance data
dist_data <- left_join(pairs_agg, agg2, by = "fire_id")

# filter mean distance at 45 m
dist_data_large <- dist_data[dist_data$dstnc_p > 45, ]
nrow(dist_data_large)
hist(dist_data_large$fwi, breaks = 30)
plot(area_ha ~ fwi, dist_data_large)
plot(log(area_ha) ~ fwi, dist_data_large)
min(dist_data_large$area_ha)

# filter mean distance at 60 m
dist_data_large <- dist_data[dist_data$dstnc_p > 60, ]
nrow(dist_data_large) # 97
hist(dist_data_large$fwi, breaks = 30)
plot(area_ha ~ fwi, dist_data_large)
plot(log(area_ha) ~ fwi, dist_data_large)
abline(coef(lm(log(area_ha) ~ fwi, data = dist_data_large)))
min(dist_data_large$area_ha)

# filter mean distance at 90 m
dist_data_large <- dist_data[dist_data$dstnc_p > 90, ]
nrow(dist_data_large) # 53
hist(dist_data_large$fwi, breaks = 30)
plot(area_ha ~ fwi, dist_data_large)
plot(log(area_ha) ~ fwi, dist_data_large)
abline(coef(lm(log(area_ha) ~ fwi, data = dist_data_large)))
min(dist_data_large$area_ha) # 110

# filter 100 ha
dist_data_large <- dist_data[dist_data$area_ha >= 100, ]
nrow(dist_data_large) # 85
hist(dist_data_large$fwi, breaks = 30)
plot(area_ha ~ fwi, dist_data_large)
plot(log(area_ha) ~ fwi, dist_data_large);
abline(coef(lm(log(area_ha) ~ fwi, data = dist_data_large)))
min(dist_data_large$area_ha)
min(dist_data_large$dstnc_p)

# filter 100 ha
dist_data_large <- dist_data[dist_data$area_ha >= 100, ]
nrow(dist_data_large) # 85
hist(dist_data_large$fwi, breaks = 30)
plot(area_ha ~ fwi, dist_data_large)
plot(log(area_ha) ~ fwi, dist_data_large);
abline(coef(lm(log(area_ha) ~ fwi, data = dist_data_large)))
min(dist_data_large$area_ha)
min(dist_data_large$dstnc_p)


# Assign different N according to size:

# filter median distance at 60 m
dist_data_large <- dist_data[dist_data$dstnc_p > 60, ]
dist_data_large <- dist_data_large[order(dist_data_large$area_ha), ]
nrow(dist_data_large) # 97
dist_data_large$area_cat <- cut(dist_data_large$area_ha,
                                breaks = c(0, 100, 1000, 10000, 100000))
dist_data_large$n_pix_max <- factor(dist_data_large$area_cat,
                                    levels = levels(dist_data_large$area_cat),
                                    labels = c("10", "30", "50", "100")) %>%
                              as.character() %>% as.numeric()


sum(dist_data_large$n_pix_max) * 2 # 6620. Sounds good



# Remove near pairs (< 45 m) and duplicated ones --------------------------

pairs_filt_0 <- pairs[pairs$dstnc_p > 45, ]
nrow(pairs_filt_0) / nrow(pairs) # 60 % remains

# remove duplicated:
pairs_filt <- do.call(
  "rbind",
  lapply(unique(pairs_filt_0$fire_id), function(f) {

    # f = "2016_70"
    d0 <- pairs_filt_0[pairs_filt_0$fire_id == f, ]
    # evaluate just one class (otherwise, all are repeated)
    d <- d0[d0$class == "edge", ]

    # duplicated point_ids
    dups <- d$point_d[duplicated(d$point_d)]

    # define which to use between the duplicated ones
    if(length(dups) == 0) return (d0)
    else {
      print(f)
      # define which pairs are not duplicated
      pairs_use <- d[!(d$point_d %in% dups), "pair_id"]

      # within the duplicated ones, choose the maximum distance pair
      pairs_too <- integer(length(dups))
      for(i in 1:length(dups)) {
        # i = 1
        dd <- d[d$point_d == dups[i], c("pair_id", "point_d", "dstnc_p")]
        pair_max_dist <- which.max(dd$dstnc_p)
        pairs_too[i] <- dd$pair_id[pair_max_dist]
      }

      # subset all valid pairs
      pairs_ok <- c(pairs_use, pairs_too)
      return(d0[d0$pair_id %in% pairs_ok, ])
    }
}))

nrow(pairs_filt) / nrow(pairs_filt_0) # there were no duplicated ones
# 56500 points


# N by fire after filtering -----------------------------------------------

pairs_length <- aggregate(pair_id ~ fire_id, pairs_filt, length)
names(pairs_length) <- c("fire_id", "n_pixels")
pairs_area <- aggregate(area_ha ~ fire_id, pairs_filt, median)
pairs_dist <- aggregate(dstnc_p ~ fire_id, pairs_filt, median)

# join data sets
all(pairs_length$fire_id == pairs_area$fire_id)
all(pairs_area$fire_id == pairs_dist$fire_id) # they are ordered

d0 <- cbind(pairs_length,
            area_ha = pairs_area$area_ha,
            distance = pairs_dist$dstnc_p)

# filter median distance > 60 m
d1 <- d0[d0$distance > 60, ]
d1 <- d1[order(d1$area_ha), ]
rownames(d1) <- NULL
nrow(d1) # 142
d1$area_cat <- cut(d1$area_ha,
                   breaks = c(0, 50, 100, 1000, 10000, 100000))
d1$n_pix_max <- factor(d1$area_cat,
                       levels = levels(d1$area_cat),
                       labels = c("0", "20", "40", "100", "200")) %>%
                as.character() %>% as.numeric()
d1$n_pix_max %>% sum %>% "*"(2) # 11600 pixels

n_fires <- aggregate(n_pix_max ~ area_cat, d1, length)
sum(n_fires$n_pix_max[-1]) # 100 fuegos


# Get the Npix most distanced pixels by fire ------------------------------

fires_use <- d1$fire_id[d1$n_pix_max > 0]

pairs_vec_flat <- project(pairs_vec, "EPSG:5343")
crds(pairs_vec)
crds(pairs_vec_flat)

keep <- (pairs_vec_flat$fire_id %in% fires_use) &
        (pairs_vec_flat$dstnc_p > 45) &
        (pairs_vec_flat$month < 5 | pairs_vec_flat$month > 11)
p <- pairs_vec_flat[keep, ]
nrow(p) # OK

# filter points in the edge by distance, keeping the most distanced
pex <- p[p$fire_id == "2015_50" & p$class == "edge", ]
dis <- distance(pex, symmetrical = TRUE)
pex_c <- crds(pex)
dd <- dist(pex_c, upper = TRUE) %>% as.matrix
diag(dd) <- NA
mindist <- apply(dd, 1, min, na.rm = T)
ppp <- pex[order(mindist), ]
ppp <- ppp[1:100, ]

p_distanced <- do.call(
  "rbind",
  lapply(unique(p$fire_id), function(f) {
    # subset points by fire
    # f = "2006_18"
    v1 <- p[p$fire_id == f, ]
    # plot(v1)
    # subset the edge class
    v2 <- v1[v1$class == "edge", ]

    # distance between pixels
    dd <- dist(crds(v2), upper = TRUE) %>% as.matrix
    diag(dd) <- NA
    mindist <- apply(dd, 1, min, na.rm = T)
    v3 <- v2[order(mindist), ]

    # Get N allowed
    NN <- d1[d1$fire_id == f, ]$n_pix_max
    v4 <- v3[1:NN, ]

    v5 <- v1[v1$pair_id %in% v4$pair_id, ]
    # plot(v5)
    return(v5)
  })
)

# p_distanced %>% nrow 10240

# project to latlong
p_dist_wgs <- project(p_distanced, "EPSG:4326")
# write again
# writeVector(p_dist_wgs, "paired_points_edge-core_subset.shp")

# unique(p_distanced$year)[order(unique(p_distanced$year))]
# aggregate(pair_id ~ fire_id, p_distanced, length)


# Filter points in large database -----------------------------------------

# When the data was downloaded previously filtering points, many points were
# NA and GEE removed them. So I apply the filters after trying to
# download all points for each fire (200 by class in all.)

pairs_vec <- vect("data_edge-core_full.shp")
pairs_vec <- project(pairs_vec, "EPSG:5343") # flat projection needed to get distances

# remove winter fires and pairs too close
pairs_vec <- pairs_vec[(pairs_vec$month < 5 | pairs_vec$month > 11) &
                        pairs_vec$dstnc_p >= 45, ]

# Get desired number of pixels by fire
pairs_length <- aggregate(pair_id ~ fire_id, pairs_vec, length)
names(pairs_length) <- c("fire_id", "pixels_available")
pairs_area <- aggregate(area_ha ~ fire_id, pairs_vec, median)
pairs_dist <- aggregate(dstnc_p ~ fire_id, pairs_vec, median)

# join data sets
all(pairs_length$fire_id == pairs_area$fire_id)
all(pairs_area$fire_id == pairs_dist$fire_id) # they are ordered

d0 <- cbind(pairs_length,
            area_ha = pairs_area$area_ha,
            distance = pairs_dist$dstnc_p)
d1 <- d0[order(d0$area_ha), ]
rownames(d1) <- NULL
nrow(d1) # 164
d1$area_cat <- cut(d1$area_ha,
                   breaks = c(0, 50, 100, 1000, 10000, 100000))
d1$n_pix_max <- factor(d1$area_cat,
                       levels = levels(d1$area_cat),
                       labels = c("0", "20", "40", "100", "200")) %>%
  as.character() %>% as.numeric()
d1$n_pix_max %>% sum %>% "*"(2) # 10680 pixels

n_fires <- aggregate(n_pix_max ~ area_cat, d1, length)
sum(n_fires$n_pix_max[-1]) # 92 fuegos


# Get desired number of pixels by fire
d1$pix_get <- apply(cbind(ceiling(d1$pixels_available / 2), d1$n_pix_max), 1, min)

# subset fires where pixels are desired
d2 <- d1[d1$pix_get > 0, ]
nrow(d2) # 92 fires :)

# In those fires, get the most distant edge pixels, associated with their
# core pixels.

p_distanced <- do.call(
  "rbind",
  lapply(unique(d2$fire_id), function(f) {
    # subset points by fire
    # f = "2006_18"
    v1 <- pairs_vec[pairs_vec$fire_id == f, ]
    # plot(v1)
    # subset the edge class
    v2 <- v1[v1$class == "edge", ]
    # plot(v2)

    # distance between pixels
    dd <- dist(crds(v2), upper = TRUE) %>% as.matrix
    diag(dd) <- NA
    mindist <- apply(dd, 1, min, na.rm = T)
    v3 <- v2[order(mindist, decreasing = TRUE), ]

    # Get N allowed
    N_best <- d2[d2$fire_id == f, ]$pix_get
    N_use <- ifelse(nrow(v3) >= N_best, N_best, nrow(v3))
    v4 <- v3[1:N_use, ]
    # plot(v4)

    # check these points have pairs paired (not NA in the core)
    pairs_use <- v4$pair_id[v4$pair_id %in% v1$pair_id[v1$class == "core"]]

    # get original data
    v5 <- v1[v1$pair_id %in% pairs_use]
    return(v5)
  })
)

nrow(p_distanced) # 9610
kk <- aggregate(pair_id ~ fire_id, p_distanced, length)
summary(kk$pair_id)
# min is 16 pixels (total)
summary(as.numeric(as.factor(p_distanced$class)) - 1) # OK, perfectly balanced.

# write data
writeVector(p_distanced, "data_edge-core_filtered.shp")
write.csv(p_distanced, "data_edge-core_filtered.csv")