##### Figure to show disequilibrium dynamics from the metacommunity simulation

source("lib/CommunityTempDis.R")
source("lib/SpeciesPoolGen.R")

set.seed(123)

burnin_yrs <- 5000
baseline_yrs <- 4000
warming_yrs <- 100
final_yrs <- 10000
obs_yrs <- baseline_yrs + warming_yrs + final_yrs
n_sp <- 50

temp <- sim.temps(burnin_yrs = burnin_yrs, baseline_yrs = baseline_yrs, warming_yrs = warming_yrs, final_yrs = final_yrs)
sim_time <- seq_along(temp[,1])

comm_sim <- sim.chalmandier(burnin_yrs = burnin_yrs, baseline_yrs = baseline_yrs, warming_yrs = warming_yrs, final_yrs = final_yrs,
                            deltaT = 3, d = .01, Gmax = 0.5, Tstdev = 2,
                            temp_niche="temp_optima", identical_temp_noise = F,
                            Tmin = 0, Tmax = 15,
                            Tmin_sp = 0, Tmax_sp = 20, N = n_sp,
                            temp_direct = temp)
time <- (comm_sim$args$burnin_yrs+1):comm_sim$args$sim_yrs
biomass_sites <- apply(comm_sim$spxp[time, , ], c(1,2), sum) # summing across species
biomass_by_species_17 <- comm_sim$spxp[time, 17,] # pulling out all species abundances at site 17
biomass_17 <- apply(biomass_by_species_17, 1, sum) # total biomass at focal site

site_pal <- colorRampPalette(c("blue","red"))(20)
# matplot(biomass_sites[seq(1,dim(biomass_sites_avg)[1], by=2),], type="l", lty=1, col=site_pal, main="total biomasses at all sites")
# matplot(biomass_by_species_17, type="l", main="species at focal site")

# simulation time series for fitting
temp_smooth <- sapply(seq(0,15,length=20), function(x) generateTemps(meanT=x,changeT=3,stationary=T, burnin_yrs, baseline_yrs, warming_yrs, final_yrs))[-c(1:burnin_yrs),]
temp_direct <- temp[-c(1:burnin_yrs),]
# plot(temp[-c(1:burnin_yrs),17], type="l")
# lines(temp_smooth[,17], col="red")

# fit acclimation model
fit <- fit.acc.model.sim(F_t = biomass_17, C_t = temp_smooth[,17], t = 1:obs_yrs,
                         l_guess = 1000, g_guess = .1, a_guess = .1, b_guess = .01, E0_guess = temp_smooth[1,17],
                         method = "BFGS",
                         form = 1, disequil_abs = T, plot=F)

time_pal_fun <- colorRampPalette(c("#ffdf3d","#ffcb3d","#fc95df","#f873ff","#3d44ff","#3d44ff"))
time_pal <- adjustcolor(time_pal_fun(obs_yrs), 0.5)

# calculating community niches
tr <- comm_sim$tr # traits
E_fund <- apply(comm_sim$spxp[time,,], c(1,2), function(x) weighted.mean(tr[,5], x)) # abundance weighted mean of temp optima (fundamental niche)
# Temp at site with maximum abundance over time (realized niche)
baseline_period <- seq(burnin_yrs, burnin_yrs+baseline_yrs)
ab_by_site <- apply(comm_sim$spxp[baseline_period,,], c(2,3), mean)
t_sites <- comm_sim$args$T_sites
site_max <- apply(ab_by_site, 2, function(x) ifelse(all(x < 1.0e-10), NA, which.max(x)))
tr_realized <- t_sites[site_max]
E_real <- apply(comm_sim$spxp[time,,], c(1,2), function(x) weighted.mean(tr_realized, x, na.rm=T))
# plot(temp[-c(1:burnin_yrs),17], type="l", col="gray", ylim=c(10,20))
# # lines(E_real[,17], ylim=range(c(E_real[,17], E_fund[,17], fit$E_pred)))
# lines(E_real[,17])
# lines(E_fund[,17], col="blue")
# lines(fit$E_pred, col="red")

# Simulation figure - Figure 3
svg("figures/simulation_demo_revision.svg",width=10, height=8)
  # layout(matrix(c(1,2,2,3,3,4,4,4,4,4), nrow=5, ncol=2), widths = c(1,1.5), heights = c(1,.85,.80,1,1))
  layout(matrix(c(1,2,2,3,3,4,4,4,5,5,4,4,4,6,6), nrow=5, ncol=3), widths = c(1,0.65,0.90), heights = c(1,.85,.80,1,1))
  par(mar=c(0,4.5,1.5,1))
  pal <- c("#f79000", "#18c900", "#00bef7")
  plot(temp[-c(1:burnin_yrs),17], type="l", col="gray", xlab="", ylab="", xaxt="n", yaxt="n")
  abline(h=temp_smooth[1,17]+2, col=pal[1], lwd=2)
  abline(h=temp_smooth[1,17], col=pal[2], lwd=2)
  abline(h=temp_smooth[1,17]-2, col=pal[3], lwd=2)
  mtext("Temperature",2,1.5)
  mtext("   a",3,-2, adj=0)
  lines(temp_smooth[,17], type="l", lwd=3)
  abline(v=baseline_yrs, lty=2, col="gray")
  
  par(mar=c(0,4.5,0,1))
  tr_fund <- comm_sim$tr[,1] # getting species niches
  t_sites <- comm_sim$args$T_sites
  tr_fund_mismatch <- tr_fund - t_sites[17] 
  
  sp_pal <- rep(pal[2], n_sp); sp_pal[which(tr_fund_mismatch > 2)] <- pal[1]; sp_pal[which(tr_fund_mismatch < -2)] <- pal[3] # color palette by species niches relative to baseline site temp
  # select_sp <- seq(1,n_sp, by=3) # species to include in plot
  select_sp <- round(seq(1,n_sp, length.out=n_sp)) # species to include in plot
  select_sp <- c(20, 22, 30, 32, 35, 37)
  # select_sp <- c(1:50)[head(order(apply(biomass_by_species_17, 2, mean)),3)] # most abundant species
  # select_sp <- c()
  
  matplot(biomass_by_species_17[,select_sp],
          type="l", lty=1, lwd=3, col=sp_pal[select_sp],
          ylab="", xlab="", xaxt="n", yaxt="n") # consistent patterns when averaged across realizations
  mtext("Species abundance",2,1.5)
  mtext("   b",3,-2, adj=0)
  abline(v=baseline_yrs, lty=2, col="gray")
  
  # text(x = nrow(biomass_by_species_17), y = tail(biomass_by_species_17[,select_sp], 1), labels = select_sp, col = sp_pal[select_sp], pos = 1)
  
  par(mar=c(4,4.5,0,1))
  plot(biomass_17, col=time_pal, xlab="", ylab="", pch=19, xaxt="n", yaxt="n", ylim=range(c(biomass_17, fit$F_pred)), cex=1.2)
  mtext("Time",1,1.5)
  mtext("Community biomass",2,1.5)
  mtext("   c",3,-2, adj=0)
  abline(v=baseline_yrs, lty=2, col="gray")
  lines(fit$F_pred, lwd=8, lty=1, col=adjustcolor("white", 0.25))
  lines(fit$F_pred, lwd=2, lty="33", col=adjustcolor("black", 0.65))
  
  par(mar=c(4,4,1.5,1))
  plot(diff(biomass_17) ~ temp_direct[-1,17], col=time_pal, pch=19,
       xlab="", ylab="", xaxt="n", yaxt="n")
  abline(h=0, lty=2); abline(v=c(head(temp_smooth[,17],1),tail(temp_smooth[,17],1)), col=c(head(time_pal_fun(2),1), tail(time_pal_fun(2),1)), lty=2, lwd=2)
  mtext("Community biomass growth rate",2,1.5)
  mtext("Temperature",1,1.5)
  mtext("Increase →  ",2,-2, adj=1)
  mtext("  ← Decrease",2,-2, adj=0)
  mtext("Warmer →  ",1,-2, adj=1)
  mtext("  ← Colder",1,-2, adj=0)
  mtext("d   ",3,-2, adj=1)
  
  par(mar=c(1,1,1,1))
  plot.new()
  mtext("Legends")
  legend("top", inset=0.05, legend=c("Mean","Variation"), title="Temperature data", lwd=c(5,2), col=c("black", "gray"),
         box.col ="transparent", bg=adjustcolor("white",0.5), cex=1.25)
  legend("center", inset=0.05, legend=c("Warmer","Near baseline","Colder"), title="Species temperature niche", lwd=5, col=pal,
         box.col ="transparent", bg=adjustcolor("white",0.5), cex=1.25)
  legend("topright", legend=NA, title="Time")
  legend("bottom", inset=0.05, legend=c("Realized", "Fundamental"), col=c("#4ed9bf", "#4175b5"), bty="n", lwd=7,
         title="Community niche type", cex=1.25)
  
  par(mar=c(4.5,4.25,0.5,1))
  plot(E_real[,17] ~ fit$E_pred, ylim=range(c(E_real[,17], E_fund[,17])), type="l", lwd=7, col="#4ed9bf",
       ylab="Actual community niche (°C)", xlab="Inferred community niche (°C)", cex.lab=1.5)
  lines(E_fund[,17] ~ fit$E_pred, lwd=7, col="#4175b5")
  mtext("   e",3,-2, adj=0)
  mtext("Warmer →      ",3,-2, adj=1)
  mtext("      ← Colder",1,-2, adj=0)
  abline(a=0,b=1)

dev.off()

# fit R2
fit$R2

# E correlation coefficients
cor(E_real[,17], fit$E_pred)
cor(E_fund[,17], fit$E_pred)
