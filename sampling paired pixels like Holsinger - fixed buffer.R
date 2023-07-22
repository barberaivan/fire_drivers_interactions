library(terra)
library(sf)
library(tidyverse)

# As the difference in predictors values increases with distance between
# pixels, I'm gonna fix the buffer distance at 90 m

# Distance still increases with fire area. Maybe it's because the st_segmentize
# produces more separated points in larger fires, so the distance between the
# edge and core is frequently larger than 90 m.
# Fix maximum distance between points in core segments at 10 m.

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
n_points <- 300 # to have backup

# project to posgar 2007/1
fires <- st_transform(fires_wgs, crs = 5343)

buffer_distance <- 90


# Get paired points in every fire polygon ----------------------------------

# (loop so if an error occurs, not everything will be lost)
points_data <- vector("list", length(fires$fire_id))

names_from_largest <- fires$fire_id[order(fires$area_ha, decreasing = TRUE)]

for(f in 75:length(fires$fire_id)) {

  # test
  # fire_id <- "2000_7" #2000_19"#"2015_50"

  # bug
  # f <- 146

  fire_id <- names_from_largest[f]
  fire <- fires[fires$fire_id == fire_id, ]
  fireg <- fire$geometry

  # get core area
  core_poly <- st_buffer(fireg, -buffer_distance)
  # plot(core_poly)

  # is there any polygon?
  dd <- core_poly %>% st_coordinates() %>% as.data.frame()
  follow <- nrow(dd) > 0

  # only continue if there is a core polygon
  if(follow) {


  # make many points along both polygons
  # edge_length <- st_length(st_cast(fireg, "MULTILINESTRING")) %>% as.numeric
  edge_points <- st_segmentize(fireg,
                               #dfMaxLength = edge_length / (n_points + 1)) %>%
                               dfMaxLength = 90) %>%
    st_coordinates() %>%
    as.data.frame() %>%
    select(X, Y) %>%
    st_as_sf(coords = c("X", "Y"))

  # core_length <- st_length(st_cast(core_poly, "MULTILINESTRING")) %>% as.numeric
  core_points_full <- st_segmentize(core_poly,
                                    dfMaxLength = 10) %>% # small distance to find points at around 90 m
    st_coordinates() %>%
    as.data.frame() %>%
    select(X, Y) %>%
    st_as_sf(coords = c("X", "Y"))

  # keep only n_points for edge (or less if you have less)
  keep <- seq(1, nrow(edge_points), length.out = min(n_points, nrow(edge_points))) %>% round %>% unique
  edge_points <- edge_points[keep, ]

  # For every edge point, find the closest core point and match it
  edge_points$pair_id <- 1:nrow(edge_points)
  edge_points$point_id <- 1:nrow(edge_points)
  edge_points$distance_pair <- NA
  core_points <- edge_points # just to initialize spdf

  print(paste("fire ", fire_id, ": ", "pairing points", sep = ""))

  for(i in 1:nrow(edge_points)) {
    # i = 286
    p <- edge_points[i, ]$geometry

    times <- 2
    ll <- 0
    while(ll < 1) {
      # print(paste("point", i, "times", times))
      p_buffer <- st_buffer(p, (buffer_distance + 20) * times)
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
}

# remove one extra slot
points_data <- points_data[1:length(fires$fire_id)]

# null ones?
nulls <- do.call(c, lapply(points_data, is.null))
points_data <- points_data[!nulls]

# merge list
points_flat <- do.call("rbind", points_data)

# is fire area related to distance between points?
test_dist <- aggregate(cbind(distance_pair, area_ha) ~ fire_id, points_flat, median)
plot(distance_pair ~ log(area_ha), test_dist)
abline(lm(distance_pair ~ log(area_ha), test_dist)) # now there is an inverse relationship
dist_test <- points_flat[points_flat$class == "edge", ]$distance_pair
hist(dist_test[dist_test < 225.53404])
plot(ecdf(dist_test))
quantile(dist_test, probs = seq(0, 1, by = 0.05))
# Later I will have to choose pairs based on distance.

# assign crs
st_crs(points_flat) <- st_crs(fires)

# reproject to latlong
points_wgs <- st_transform(points_flat, crs = st_crs(fires_wgs))

# write to disk
# st_write(points_wgs, "paired_points_edge-core_buffer90.shp", append = FALSE)


# Filter points in large database -----------------------------------------

# When the data was downloaded previously filtering points, many points were
# NA and GEE removed them. So I apply the filters after trying to
# download all points for each fire (200 by class in all.)

pairs_vec <- vect("data_edge-core_full_buffer90.shp")
# pairs_vec <- vect("data_edge-core_full.shp")
pairs_vec <- project(pairs_vec, "EPSG:5343") # flat projection needed to get distances

# remove winter fires and pairs too far apart
quantile(pairs_vec$dstnc_p, probs = seq(0, 1, by = 0.05))
distance_threshold <- 130
# tt <- aggregate(cbind(area_ha, dstnc_p) ~ fire_id,
#                 pairs_vec[pairs_vec$dstnc_p <= distance_threshold, ], median)
# plot(dstnc_p ~ log(area_ha), tt)
# abline(lm(dstnc_p ~ log(area_ha), tt))

pairs_vec <- pairs_vec[(pairs_vec$month < 5 | pairs_vec$month > 11) &
                        pairs_vec$dstnc_p <= distance_threshold, ]
# hist(pairs_vec$dstnc_p)
nrow(pairs_vec)

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
nrow(d2) # 93 fires :)

# In those fires, get the most distant core pixels, associated with their
# edge pixels. For that, make a buffer in each point and count how many core
# pixels intersect. Start choosing the ones that intersect less pixels.
p_distanced <- do.call(
  "rbind",
  lapply(unique(d2$fire_id), function(f) {
    # subset points by fire
    # f = "2015_82"
    print(f)

    points_all <- pairs_vec[pairs_vec$fire_id == f, ]
    points_all <- points_all[order(points_all$class, points_all$pair_id), ]

    # plot(points_all)
    # subset the edge class
    core_points_all <- points_all[points_all$class == "core", ]
    edge_points_all <- points_all[points_all$class == "edge", ]
    # plot(core_points_all)

    # buffers around points
    pbuf <- buffer(core_points_all, 150)

    # evaluate disjoint pairs
    disjoint_mat <- relate(pbuf, relation = "disjoint")

    # touches data
    disjoints <- colSums(disjoint_mat)
    touch <- data.frame(mat_id = 1:length(disjoints),
                        pair_id = core_points_all$pair_id,
                        disjoints = disjoints,
                        inside = NA)
    touch <- touch[order(touch$disjoints, decreasing = TRUE), ]

    # keep the most separated ones
    for(i in 1:nrow(touch)) {
      # i = 1
      if(is.na(touch$inside[i])) {

        # does it have a pair in the edge?
        ee <- nrow(edge_points_all[edge_points_all$pair_id == touch$pair_id[i], ])
        if(ee < 1) {
          touch$inside[i] <- FALSE
        } else {
          # keep it
          touch$inside[i] <- TRUE

          # put out the ones that the this polygon touches
          mat_in <- touch$mat_id[i]
          touched <- which(!disjoint_mat[, mat_in])
          touched <- touched[touched != mat_in]
          # reject the touched ones
          touch$inside[touch$mat_id %in% touched] <- FALSE
        }
      }
    }

    # filter
    candidates <- touch$pair_id[touch$inside]

    # Get N allowed
    N_best <- d2[d2$fire_id == f, ]$pix_get
    N_use <- ifelse(length(candidates) >= N_best, N_best, length(candidates))
    pairs_keep <- candidates[1:N_use]

    points_use <- points_all[points_all$pair_id %in% pairs_keep, ]
    ## check
    # plot(pbuf)
    # plot(pbuf[touch$mat_id[touch$inside], ], add = TRUE, col = 2)
    # plot(points_use, add = TRUE, col = 4)

    return(points_use)
  })
)

nrow(p_distanced) # 5768
kk <- aggregate(pair_id ~ fire_id, p_distanced, length)
summary(kk$pair_id)
# min is 2 pixels (total)
summary(as.numeric(as.factor(p_distanced$class)) - 1) # OK, perfectly balanced.

# write data
writeVector(p_distanced, "data_edge-core_filtered_buffer90.shp",
            overwrite = TRUE)
write.csv(p_distanced, "data_edge-core_filtered_buffer90.csv")

# check data
p1 <- aggregate(pair_id ~ fire_id, p_distanced, length)
names(p1) <- c("fire_id", "pixels_available")
p2 <- aggregate(cbind(fwi, area_ha, dstnc_p) ~ fire_id, p_distanced, median)
p3 <- cbind(p1, p2[, -1])
plot(area_ha ~ fwi, p3)
abline(lm(area_ha ~ fwi, p3))
View(p3)


