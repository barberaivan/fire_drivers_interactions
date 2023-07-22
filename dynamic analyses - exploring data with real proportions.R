options(scipen = 999)
d <- read.csv("data_dynamic_burn_prob.csv")
d_year <- aggregate(burned ~ year, d, mean)
d_year_sum <- aggregate(burned ~ year, d, sum)
# no hay nada de burned. Hay que hacer muestreo estratificado y trabajar con
# datos desproporcionados en realción a lo real.
sum(d_year_sum$burned) / nrow(d) # 0.002670833 es la annual fire prob
