#
#
#	Help function to run community dynamics
#	Main script author: Loic Chalmandrier - lchalman@uwyo.edu, https://github.com/LoicChr
#	Co-authors (and project leaders): J. Alexander - jake.alexander@unil.ch & L. Pellissier - loic.pellissier@usys.ethz.ch

# Translated from Matlab to R by Peter Adler - peter.adler@usu.edu

# Modified by Michael Stemkovski to allow hump-shaped species niches in place of cold tolerance - m.stemkovski@gmail.com

CommunityTempDis <- function(spxp_ini, tr, seed_rain, temp, dt, temp_niche="cold_threshold"){

  N <- dim(tr)[1]
  L_land <- NROW(spxp_ini)
  spxp_int <- seed_rain%*%spxp_ini
  Biom <- rowSums(spxp_int)
  
  # growth rate modified either by distance to cold threshold or to hump shaped temp optima - MS
  if(temp_niche == "cold_threshold") growth_rate_temp <- (matrix(1,L_land,1)%*%t(tr[,2]))*(temp%*%matrix(1,1,N) - matrix(1,L_land,1)%*%t(tr[,1]))
  if(temp_niche == "temp_optima") growth_rate_temp <- (matrix(1,L_land,1) %*% t(tr[,2])) * (exp(-(1/2)*((( temp%*%matrix(1,1,N) - matrix(1,L_land,1)%*%t(tr[,5]))^2) / (matrix(1,L_land,1)%*%t(tr[,6])^2))))
  growth_rate_intra <- (matrix(1,L_land,1)%*%t(tr[,4]))*spxp_ini
  growth_rate_biom <- (matrix(1,L_land,1)%*%t(tr[,3]))*(Biom%*%matrix(1,1,N))
   
  spxp_final <- spxp_int + dt*spxp_int*(growth_rate_temp - growth_rate_intra - growth_rate_biom)
  #spxp_final <- spxp_int + 1*spxp_int*(growth_rate_temp - growth_rate_intra - growth_rate_biom)
   
  spxp_final[spxp_final<10e-12] <- 0

  return(spxp_final)
  
}

generateTemps <- function(meanT,changeT,stationary=T, burnin_yrs, baseline_yrs, warming_yrs, final_yrs){
  
  if(stationary==T){
    
    out <- meanT + c(rep(0,burnin_yrs+baseline_yrs),seq(changeT/warming_yrs,changeT,changeT/warming_yrs),
                        rep(changeT,final_yrs))
  }else{
    
    out <- meanT + seq(changeT/sim_yrs,changeT,changeT/sim_yrs)
    
  }
  
  return(out)
  
}














### working to modify CommunityTempDis
# 
# L_land = 20; Tmin = 0; Tmax = 15; Tstdev = 2; deltaT = 2
# burnin_yrs = 2000; baseline_yrs = 2000; warming_yrs = 100; final_yrs = 4000;
# stationary_periods = T
# N = 40; Gmax = 0.5; Gmin = 0.2; Lmax = 0.7; Lmin = 0.7; Cmax = 0.2; Cmin = 0.2; d = 0.01
# 
# Tmean <-seq(Tmin,Tmax,length=L_land)   # baseline mean temperatures across the landscape
# sim_yrs <- burnin_yrs + baseline_yrs + warming_yrs + final_yrs # total number of years
# 
# # generate one set of temperature time series to use in all simulations
# # rows are years, cols are sites
# temperature<-matrix(NA,sim_yrs,length(Tmean))
# for(iT in 1:length(Tmean)){
#   tmp_mean <- generateTemps(meanT=Tmean[iT],changeT=deltaT,stationary=stationary_periods, burnin_yrs, baseline_yrs, warming_yrs, final_yrs)
#   temperature[,iT] = rnorm(sim_yrs,tmp_mean,Tstdev)  # generate temperature time series
# }
# 
# image(temperature)
# plot(temperature[,20], type="l")
# lines(tmp_mean, col="red")
# 
# # set up landscape
# landscape <- cbind(1:L_land , rep(1 ,L_land),Tmean)
# plot(landscape[,3]) # mean baseline temps at sites
# 
# # 1st step: initialisation of the metacommunity 
# spxp <- array(NA, dim=c(sim_yrs, L_land, N)) # [years, sites, species]
# spxp[1,,] <- 2/N
# dt <- 0.025   
# d_ins <- d/dt # this used to be d <- d/dt. MS changed it to d_in to avoid overwriting
# 
# image(spxp[1,,]) # equal abundances of all species at all sites at year 1
# 
# # Generation of the species pool
# tr <- SpeciesPoolGen(N, Tmin, Tmax+deltaT, Gmin, Gmax, Lmin, Lmax, Cmin, Cmax)
# 
# # tr[,1]: cold tolerance threshold (Tmin_i)
# # tr[,2]: growth rate (g_i)
# # tr[,3]: interspecific competition (l_i)
# # tr[,4]: intraspecific competition (c_i)
# # if temp_niche = "temp_optima"
# # tr[,5]: temp niche mean (mu)
# # tr[,6]: temp niche deviance (sigma)
# 
# plot(tr[,1]) # cold tolerance threshold (Tmin_i)
# plot(tr[,1:2]) # cold tolerance and growth rate are correlated (g_i)
# plot(tr[,c(1,3)]) # all species equally sensitive to interspecific competition (l_i)
# plot(tr[,c(1,4)]) # all species equally sensitive to intraspecific competition (c_i)
# 
# # dispersal matrix (seeds disperse only to nearest neighbor locations)
# seed_rain <- as.matrix(dist(landscape[,1:2],diag=T,upper=T))
# seed_rain[seed_rain > 1] <- 0
# seed_rain <- seed_rain*(d_ins*dt*0.5)
# diag(seed_rain) <- 1-d_ins*dt
# seed_rain[1,1] = 1-dt*d_ins/2  
# seed_rain[L_land,L_land] = 1-dt*d_ins/2
# 
# CommunityTempDis <- function(spxp_ini, tr, seed_rain, temp, dt, temp_niche="cold_threshold"){
#   
#   N <- dim(tr)[1]
#   L_land <- NROW(spxp_ini)
#   spxp_int <- seed_rain%*%spxp_ini
#   Biom <- rowSums(spxp_int)
#   
#   if(FALSE){
#     tr <- SpeciesPoolGen(N, Tmin, Tmax+deltaT, Gmin, Gmax, Lmin, Lmax, Cmin, Cmax)
#     
#     # manually adding niche mu and sigma to test calculation
#     tr <- cbind(tr, sort(runif(N, 0, 15))) # mu
#     tr <- cbind(tr, rep(5,N)) # sigma
#   }
#   
#   # this is the modification
#   if(temp_niche == "cold_threshold") growth_rate_temp <- (matrix(1,L_land,1)%*%t(tr[,2]))*(temp%*%matrix(1,1,N) - matrix(1,L_land,1)%*%t(tr[,1]))
#   if(temp_niche == "temp_optima"){
#     # exp(-(1/2)*(((x-mu)^2)/(sigma^2))) # gaussian curve with peak value = 1
#     
#     growth_rate_temp <- (matrix(1,L_land,1) %*% t(tr[,2])) * (exp(-(1/2)*((( temp%*%matrix(1,1,N) - matrix(1,L_land,1)%*%t(tr[,5]))^2) / (matrix(1,L_land,1)%*%t(tr[,6])^2))))
#   } 
#   growth_rate_intra <- (matrix(1,L_land,1)%*%t(tr[,4]))*spxp_ini
#   growth_rate_biom <- (matrix(1,L_land,1)%*%t(tr[,3]))*(Biom%*%matrix(1,1,N))
#   
#   spxp_final <- spxp_int + dt*spxp_int*(growth_rate_temp - growth_rate_intra - growth_rate_biom)
#   
#   spxp_final[spxp_final<10e-12] <- 0
#   
#   return(spxp_final)
#   
# }
# 
# # Loop through time
# for(t in 2:sim_yrs){
#   if(FALSE){
#     t <- 2
#     spxp_ini <- spxp[t-1,,]
#     temp <- yr_temp
#   }
# 
#   
#   yr_temp <- temperature[t,]
#   
#   spxp[t,,] <- CommunityTempDis(spxp[t-1,,], tr, seed_rain, yr_temp, dt)
# } 
# 
# 
# 
# 
# 
# x <- 1:100
# mu <- 40
# sigma <- 10
# 
# gaussian <- (1/(sigma*sqrt(2*pi)))*exp(-(1/2)*(((x-mu)^2)/(sigma^2)))
# plot(gaussian)
# gaussian[40]
# 
# gaussian <- exp(-(1/2)*(((x-mu)^2)/(sigma^2)))
# plot(gaussian)
# gaussian[40]
# 
# 
# 
# 
# # math for the new niche optima formulation:
# # P_{i,j,t+1} = P'_{i,j,t} + \Delta t \times P'_{i,j,t}[g_i e^{-\frac{1}{2}(\frac{T_j-\mu_i}{\sigma_i})^2} -c_i P'_{i,j,t} - l_i \sum_{k}^{} P'_{k,j,t}]
