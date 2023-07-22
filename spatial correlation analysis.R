library(terra)
library(sptotal)

pll <- vect("data_random_points_spatial_corr.shp")
pll <- pll[sample(1:30000, size = 10000), ]
p <- project(pll, "EPSG:5343")
pd <- values(p)
vars <- colnames(pd)[colnames(pd) != "vegetation"]

pd <- cbind(pd, crds(p))
head(pd)

sv_list <- vector("list", length(vars))
names(sv_list) <- vars

for(i in 1:length(vars)) {
  print(vars[i])
  sv_list[[i]] <- sv(pd, which(names(pd) == "x"), which(names(pd) == "y"),
                     vars[i], bins = 100, cutoff = 20000)
  plot(gamma ~ dist, sv_list[[i]], main = vars[i])
  gc()
}

for(v in vars) plot(gamma ~ dist, sv_list[[v]], main = v)
sv_list$Rough450$np
# very high correlation always