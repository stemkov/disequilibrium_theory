
source("functions.R")
`%!in%` <- Negate(`%in%`)
library(data.table)
library(fields)

#### Main acclimation and function time series figures

# community climate niches over time

time = c(1:300)
clim_assym = 200
clim_slope = .05
clim_0 = 100
ecosystem_0 = 100
service_slope = 1
service_0 = 1000
lambda = 50
gamma = -0.02

climate <- acclimation.response.func(time, clim_assym = clim_assym, clim_slope = clim_slope, clim_0 = clim_0,
                                     E0 = ecosystem_0, beta = service_slope, alpha = service_0,
                                     lambda = lambda, gamma = gamma,
                               what_return = "climate")

time_start <- 51
climate_sub <- climate[time_start:300]
climate_noisy <- climate_sub +  rnorm(length(climate_sub), 0, 10)

par(mfrow=c(1,2))

par(mar=c(3,3,1,1))

plot(climate_noisy ~ time[time_start:300], type="l", xlab="", ylab="", xaxt="n", yaxt="n")
lines(climate_sub ~ time[time_start:300], lwd=2)
title(ylab="Climate", line=1, cex.lab=1.4)
title(xlab="Time", line=1, cex.lab=1.4)

pal <- colorRampPalette(c("red", "blue"), space="Lab", interpolate = "spline")
cols <- pal(6)

lambdas <- c(10, 30, 50, 100, 200, 500)

for(i in seq_along(lambdas)){
  ecosystem <- acclimation.response.func(time = c(time_start:300),
                                         clim_assym = clim_assym, clim_slope = clim_slope, clim_0 = clim_0,
                                         E0 = ecosystem_0, beta = service_slope, alpha = service_0,
                                         lambda = lambdas[i], gamma = gamma,
                                         what_return = "ecosystem", climate_direct = climate_noisy)
  lines(ecosystem ~ time[time_start:300], col=cols[i], lwd=3)
}

legend("topleft", inset=0.03, legend = lambdas, lty=1, lwd=3, col=cols, title=expression(lambda), cex=1.1,bty="n")

# ecosystem function over time

lambdas <- c(10, 30, 50, 100, 200, 500)
gammas <- c(-.02, -0.018, -0.022)

params <- expand.grid(lambda = lambdas, gamma = gammas)
sims <- mapply(acclimation.response.func, gamma = params$gamma, lambda = params$lambda,
               MoreArgs = list(time = time[time_start:300], clim_assym = 200, clim_slope = .05,
                               clim_0 = 100, E0 = 100,
                               beta = 1, alpha = 1000, plot=F,
                               what_return = "response", climate_direct = climate_sub),
               SIMPLIFY = F)


plot(NA, xlim=c(1,250), ylim=range(sims), xlab="", ylab="", xaxt="n", yaxt="n")
mapply(lines, sims[1:6], col=cols, lwd=3)

draw.poly <- function(time,curve1, curve2, col){
  polygon(c(time, rev(time)),
          c(curve1, rev(curve2)),
          col = adjustcolor(col, 0.3), border=NA)
}
mapply(draw.poly, curve1 = sims[7:12], curve2 = sims[13:18], col=cols, MoreArgs = list(time=c(1:250)))

abline(h=1100, lty=2)
title(ylab="Ecosystem function", line=1, cex.lab=1.4)
title(xlab="Time", line=1, cex.lab=1.4)

legend("bottomleft", 0.1, legend=c(expression(gamma)), lty=1, lwd=3, bty="n")






