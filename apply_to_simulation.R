### Fitting to simulated data and making figure

source("lib/CommunityTempDis.R")
source("lib/SpeciesPoolGen.R")


set.seed(1)

n_sites <- 25

sim.temps <- function(Tmin = 0, Tmax = 15, Tstdev = 2, L_land = n_sites,
                      burnin_yrs = 5000, baseline_yrs = 4000, warming_yrs = 100, final_yrs = 10000,
                      deltaT = 3, stationary_periods = T, identical_temp_noise = F){
  sim_yrs <- burnin_yrs + baseline_yrs + warming_yrs + final_yrs
  Tmean <-seq(Tmin,Tmax,length=L_land)
  temperature<-matrix(NA,sim_yrs,length(Tmean))
  temp_noise <- rnorm(sim_yrs, 0, Tstdev)
  for(iT in 1:length(Tmean)){
    tmp_mean <- generateTemps(meanT=Tmean[iT],changeT=deltaT,stationary=stationary_periods, burnin_yrs, baseline_yrs, warming_yrs, final_yrs)
    
    if(identical_temp_noise){
      temperature[,iT] = tmp_mean + temp_noise
    } else{
      temperature[,iT] = rnorm(sim_yrs,tmp_mean,Tstdev)  
    }
  }
  return(temperature)
}
temp <- sim.temps()

reps <- 100

biomass_sites_all <- array(NA, dim=c(14100, n_sites, reps))
community_niches_all <- array(NA, dim=c(14100, n_sites, reps,2)) # 1 is fundamental, 2 is realized

for(i in 1:reps){
  comm_sim <- sim.chalmandier(burnin_yrs = 5000, baseline_yrs = 4000, warming_yrs = 100, final_yrs = 10000,
                              deltaT = 3, d= .01, Gmax = 0.5, Tstdev = 2,
                              temp_niche="temp_optima", identical_temp_noise = F,
                              Tmin = 0, Tmax = 18,
                              Tmin_sp = 0, Tmax_sp = 20, N = 50, L_land = n_sites,
                              temp_direct = NA) # set temp_direct to NA for independent reps, temp for same climate across reps
  burnin_yrs <- comm_sim$args$burnin_yrs
  sim_yrs <- comm_sim$args$sim_yrs
  baseline_yrs <- comm_sim$args$baseline_yrs
  baseline_period <- seq(burnin_yrs, burnin_yrs+baseline_yrs)
  time <- (comm_sim$args$burnin_yrs+1):comm_sim$args$sim_yrs

  # calculating community niches
  tr <- comm_sim$tr # traits
  E_fund <- apply(comm_sim$spxp[time,,], c(1,2), function(x) weighted.mean(tr[,5], x)) # abundance weighted mean of temp optima (fundamental niche)
  # Temp at site with maximum abundance over time (realized niche)
  ab_by_site <- apply(comm_sim$spxp[baseline_period,,], c(2,3), mean)
  t_sites <- comm_sim$args$T_sites
  site_max <- apply(ab_by_site, 2, function(x) ifelse(all(x < 1.0e-10), NA, which.max(x)))
  #site_max
  tr_realized <- t_sites[site_max]
  E_real <- apply(comm_sim$spxp[time,,], c(1,2), function(x) weighted.mean(tr_realized, x, na.rm=T))
  community_niches_all[,,i,1] <- E_fund
  community_niches_all[,,i,2] <- E_real
  
  biomass_sites <- apply(comm_sim$spxp[time, , ], c(1,2), sum)

  biomass_sites_all[,,i] <- biomass_sites
  print(i)
}
saveRDS(community_niches_all, "community_niches_all_fig4_25sites.RDS")
saveRDS(biomass_sites_all, "biomass_sites_all_fig4_25sites.RDS")

community_niches_all <- readRDS("community_niches_all_fig4_25sites.RDS")
biomass_sites_all <- readRDS("biomass_sites_all_fig4_25sites.RDS")


focal_site <- 16
biomass_all_site <- biomass_sites_all[,focal_site,]
biomass_quants <- apply(biomass_all_site, 1, function(x) quantile(x, probs = c(0.01, 0.5, 0.99)))
time <- 1:ncol(biomass_quants)
temp_smooth <- sapply(seq(0,18,length=n_sites), function(x) generateTemps(meanT=x,changeT=3,stationary=T, 5000, 4000, 100, 10000))[-c(1:5000),]

E_fund_quants <- apply(community_niches_all[,focal_site,,1], 1, function(x) quantile(x, probs = c(0.01, 0.5, 0.99)))
E_real_quants <- apply(community_niches_all[,focal_site,,2], 1, function(x) quantile(x, probs = c(0.01, 0.5, 0.99)))


# back-getting 100 simulated temperature landscapes for plotting
temp_all <- array(NA, dim=c(14100, n_sites, reps))
for(i in 1:reps){
  temp_all[,,i] <- sim.temps(L_land = n_sites, Tmax = 18)[5001:19100,]
}
temp_all_site <- temp_all[,focal_site,]
temp_quants <- apply(temp_all_site, 1, function(x) quantile(x, probs = c(0.01, 0.5, 0.99)))

sim_time <- 1:14100

fit <- fit.acc.model.sim(F_t = biomass_quants[2,], C_t = temp_smooth[,focal_site], t = sim_time,
                         l_guess = 1000, g_guess = .1, a_guess = .1, b_guess = .01, E0_guess = temp_smooth[1,focal_site],
                         method = "BFGS",
                         form = 9, disequil_sqr = T, plot=F)

svg("/home/michael/Documents/Postdoc - Ecosystem Acclimation/Research Projects/acclimation_theory/figures/fitting_data.svg",
    width=3.5, height=7)

  plot_time <- 3000:8000
  time_axis <- seq_along(plot_time)
  layout(matrix(c(1,2,2,2,3,3,3),ncol=1))
  
  par(mar=c(0,4,1,1))
  plot(temp_smooth[plot_time,focal_site] ~ plot_time, type="l", col="white", xlab="", ylab="Temperature (°C)", xaxt="n", ylim=range(temp_quants))
  mtext("   a",3,-2, adj=0)
  polygon(x=c(plot_time, rev(plot_time)),
          y=c(temp_quants[1,plot_time], rev(temp_quants[3,plot_time])),
          col=adjustcolor("black",0.75), border=NA, angle=-60, density=20)
  lines(temp_smooth[plot_time, focal_site] ~ plot_time, col="black", lwd=3)
  
  par(mar=c(4,4,0,1))
  plot(biomass_quants[2,plot_time] ~ time_axis, type="l", ylim=range(c(biomass_quants)), col="white", xlab="Time", ylab="Community biomass")
  polygon(x=c(time_axis, rev(time_axis)),
          y=c(biomass_quants[1,plot_time], rev(biomass_quants[3,plot_time])),
          col=adjustcolor("#ff78b2",0.75), border=NA, angle=-60, density=20)
  lines(biomass_quants[2,plot_time] ~ time_axis, lwd=8, col="#ff78b2")
  lines(fit$F_pred[plot_time] ~ time_axis, lwd=4, col="#2d3bcf", lty="22")
  legend("bottomright", inset=0.05, legend=c("Simulated data", "Predicted"), col=c("#ff78b2", "#2d3bcf"), bty="n", lwd=5, lty=c(1,1))
  mtext("   b",3,-2, adj=0)
  
  par(mar=c(4,4,1,1))
  plot(fit$E_pred[plot_time] ~ E_real_quants[2,plot_time], col="#ffa600", type="l", lwd=5,
       xlim=range(c(E_real_quants, E_fund_quants)),
       xlab="Actual community climate niche (°C)", ylab="Inferred community climate niche (°C)")
  polygon(x=c(E_real_quants[1,plot_time], rev(E_real_quants[3,plot_time])),
          y=c(fit$E_pred[plot_time], rev(fit$E_pred[plot_time])),
          col=adjustcolor("#ffa600",0.75), border=NA, angle=-60, density=20)
  lines(fit$E_pred[plot_time] ~ E_fund_quants[2,plot_time], col="#00f2ff", lwd=5)
  polygon(x=c(E_fund_quants[1,plot_time], rev(E_fund_quants[3,plot_time])),
          y=c(fit$E_pred[plot_time], rev(fit$E_pred[plot_time])),
          col=adjustcolor("#00f2ff",0.75), border=NA, angle=60, density=20)
  legend("bottomright", inset=0.05, legend=c("Realized", "Fundamental"), col=c("#ffa600", "#00f2ff"), bty="n", lwd=5)
  mtext("   c",3,-2, adj=0)

dev.off()

