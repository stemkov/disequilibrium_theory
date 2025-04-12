##### Figure to show disequilibrium dynamics from the metacommunity simulation

source("lib/CommunityTempDis.R")
source("lib/SpeciesPoolGen.R")

set.seed(1)

# to pass same temperature to every replicate (direct)
sim.temps <- function(Tmin = 0, Tmax = 15, Tstdev = 2, L_land = 20,
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
sim_time <- seq_along(temp[,1])

# smooth: climate independent in every run - stochasticity and interannual variability averaged away
# direct: climate same in every run - only stochasticity averaged away
reps <- 100
biomass_sites_all <- array(NA, dim=c(14100, 20, reps))
biomass_by_species_17 <- array(NA, dim=c(14100, 50, reps))
for(i in 1:reps){
  comm_sim <- sim.chalmandier(burnin_yrs = 5000, baseline_yrs = 4000, warming_yrs = 100, final_yrs = 10000,
                              deltaT = 3, d= .01, Gmax = 0.5, Tstdev = 2,
                              temp_niche="temp_optima", identical_temp_noise = F,
                              Tmin = 0, Tmax = 15,
                              Tmin_sp = 0, Tmax_sp = 20, N = 50,
                              temp_direct = temp) # set temp_direct to NA for independent reps
  time <- (comm_sim$args$burnin_yrs+1):comm_sim$args$sim_yrs
  biomass_sites <- apply(comm_sim$spxp[time, , ], c(1,2), sum) # summing across species
  biomass_sites_all[,,i] <- biomass_sites
  biomass_by_species_17[,,i] <- comm_sim$spxp[time, 17, ] # pulling out all species abundances at site 17
  print(i)
}
saveRDS(biomass_sites_all, "biomass_sites_all_temp_independent.RDS")
saveRDS(biomass_sites_all, "biomass_sites_all_temp_direct.RDS")
saveRDS(temp, "temp_direct.RDS")
saveRDS(biomass_by_species_17, "biomass_by_species_17.RDS")
biomass_sites_all <- readRDS("biomass_sites_all_temp.RDS")

biomass_sites_avg <- apply(biomass_sites_all, c(1,2), mean)
site_pal <- colorRampPalette(c("blue","red"))(20)
matplot(biomass_sites_avg[seq(1,dim(biomass_sites_avg)[1], by=2),], type="l", lty=1, col=site_pal)

# simulation time series for fitting
biomass_smooth <- apply(readRDS("biomass_sites_all_temp_independent.RDS"), c(1,2), mean)
biomass_direct <- apply(readRDS("biomass_sites_all_temp_direct.RDS"), c(1,2), mean)
temp_smooth <- sapply(seq(0,15,length=20), function(x) generateTemps(meanT=x,changeT=3,stationary=T, 5000, 4000, 100, 10000))[-c(1:5000),]
temp_direct <- readRDS("temp_direct.RDS")[-c(1:5000),]
biomass_by_species_17 <- readRDS("biomass_by_species_17.RDS")

time_pal_fun <- colorRampPalette(c("#ffc800","pink","purple"))
time_pal_fun <- colorRampPalette(c("#ffdf3d","#ffcb3d","#fc95df","#f873ff","#3d44ff","#3d44ff"))
time_pal <- adjustcolor(time_pal_fun(14100), 0.5)

tr_fund <- comm_sim$tr[,1] # getting species niches
tr_fund_mismatch <- tr_fund - t_sites[17] 
biomass_by_species_17_mean <- apply(biomass_by_species_17, c(1,2), mean) # taking mean for all species across realizations

# figure for box in theory paper - including species abundance changes
svg("/home/michael/Documents/Postdoc - Ecosystem Acclimation/Research Projects/acclimation_theory/figures/simulation_curves4.svg",
    width=9, height=6)
  layout(matrix(c(1,2,2,3,3,4,4,4,4,4), nrow=5, ncol=2), widths = c(1,1.5), heights = c(1,.85,.80,1,1))
  par(mar=c(0,4.5,1.5,1))
  pal <- c("#f79000", "#18c900", "#00bef7")
  plot(temp_direct[,17], type="l", col="gray", xlab="", ylab="", xaxt="n", yaxt="n")
  abline(h=temp_smooth[1,17]+2, col=pal[1], lwd=2)
  abline(h=temp_smooth[1,17], col=pal[2], lwd=2)
  abline(h=temp_smooth[1,17]-2, col=pal[3], lwd=2)
  mtext("Temperature",2,1.5)
  mtext("   a",3,-2, adj=0)
  lines(temp_smooth[,17], type="l", lwd=3)
  
  par(mar=c(0,4.5,0,1))
  tr_fund <- comm_sim$tr[,1] # getting species niches
  tr_fund_mismatch <- tr_fund - t_sites[17] 

  sp_pal <- rep(pal[2], 50); sp_pal[which(tr_fund_mismatch > 2)] <- pal[1]; sp_pal[which(tr_fund_mismatch < -2)] <- pal[3] # color palette by species niches relative to baseline site temp
  select_sp <- seq(1,50, by=3) # species to include in plot
  matplot(biomass_by_species_17_mean[,select_sp],
          type="l", lty=1, lwd=3, col=sp_pal[select_sp],
          ylab="", xlab="", xaxt="n", yaxt="n") # consistent patterns when averaged across realizations
  legend("topleft", inset=0.05, legend=c("warmer","near baseline","colder"), title="Temperature niche", lwd=5, col=pal, bty="n")
  mtext("Species abundance",2,1.5)
  mtext("   b",3,-2, adj=0)
  
  par(mar=c(4,4.5,0,1))
  plot(biomass_direct[,17], col=time_pal, xlab="", ylab="", pch=19, xaxt="n", yaxt="n")
  mtext("Time",1,1.5)
  mtext("Community biomass",2,1.5)
  mtext("   c",3,-2, adj=0)
  
  par(mar=c(4,4,1.5,1))
  plot(diff(biomass_direct[,17]) ~ temp_direct[-1,17], col=time_pal, pch=19,
       xlab="", ylab="", xaxt="n", yaxt="n")
  abline(h=0, lty=2); abline(v=c(head(temp_smooth[,17],1),tail(temp_smooth[,17],1)), col=c(head(time_pal_fun(2),1), tail(time_pal_fun(2),1)), lty=2, lwd=2)
  mtext("Change in biomass per unit time",2,1.5)
  mtext("Temperature",1,1.5)
  mtext("Increase →  ",2,-2, adj=1)
  mtext("  ← Decrease",2,-2, adj=0)
  mtext("d   ",3,-2, adj=1)

dev.off()
