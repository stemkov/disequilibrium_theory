### functions to be called from other scripts


`%!in%` <- Negate(`%in%`)

# Community niche acclimating to climate
# dE/dt (as a difference equation)
# acclimation.func <- function(time, climate, lambda, E_0){
#   E <- rep(NA, length(time))
#   E[1] <- E_0
#   for(i in seq_along(time)[-1]){
#     t <- time[i]
#     E[i] <- E[i-1] - (1/lambda) * (E[i-1] - climate[i-1])
#   }
#   return(E)
# }

acclimation.func <- function(time, climate, lambda, E_0, kappa=0){
  E <- rep(NA, length(time))
  E[1] <- E_0
  for(i in seq_along(time)[-1]){
    t <- time[i]
    E[i] <- E[i-1] - (1/lambda) * (E[i-1] - climate[i-1])
  }
  return(E)
}


# Ecological response to climate and disequilibrium
# F_t
response.func <- function(climate, ecosystem, gamma = 0.02, alpha=1000, beta=1, disequil_abs = FALSE, disequil_sqr = FALSE, form=1){
  disequil <- ecosystem - climate
  if(disequil_abs) disequil <- abs(disequil)
  if(disequil_sqr) disequil <- disequil^2
  
  if(form == 0) response <- alpha + beta*climate # null linear model - no effect of disequilibrium
  # response <- alpha + (beta - gamma*disequilibrium)*climate
  if(form == 1) response <- alpha + beta*climate - gamma*disequil*climate # disequilibrium affects slope
  if(form == 2) response <- alpha + beta*climate - gamma*disequil # disequilibrium affects intercept
  if(form == 3) response <- alpha - gamma*disequil # there is no climate sensitivity, baseline response is a result of disequilibrium - to be used with multiple acclimation processes
  if(form == 4) response <- alpha - gamma*disequil*climate # there is no direct climate sensitivity, climate sensitivity is a result of disequilibrium - to be used with multiple acclimation processes
  
  if(form == 5) response <- alpha*ecosystem - gamma*disequil # baseline response is defined by ecosystem structure
  if(form == 6) response <- alpha*ecosystem + beta*climate - gamma*disequil*climate # baseline response is result of structure, and there is baseline climate sensitivity, and disequilibrium affects climate sensitivity
  
  if(form == 7) response <- alpha + (beta*climate)/(gamma*(1+disequil)) # climate sensitivity is reduced by disequilibrium. +1 so that sensitivity doesn't approach infinity if disequil is small
  if(form == 8) response <- alpha/(gamma*(1+disequil)) + beta*climate # baselines function is reduced by disequilibrium, sensitivity is unaffected
  
  if(form == 9) response <- (alpha + beta*climate) * exp(-(disequil)/gamma) # photosynthesis efficiency acclimation from Friend (2010)
  if(form == 10) response <- alpha + (beta*climate) * exp(-(disequil)/gamma) # photosynthesis efficiency acclimation from Friend (2010)
  
  if(form == 11) response <- alpha + (beta*climate) * -exp(disequil/gamma)^2 # squares the whole exp term - to be used with disequil_sqr=F - more similar to Friend 2010
  
  return(response)
}


# Data simulation function
# climate changes logistically
# ecosystem acclimates
# ecological response to climate and disequilibrium
# acclimation.response.func <- function(time, clim_assym, clim_slope, clim_0, E0,
#                                       beta, alpha,
#                                       clim_noise_sd=0, noise_vec = NA,
#                                       lambda, gamma, plot=F,
#                                       what_return = "response", climate_direct=NA, disequil_abs=FALSE, form=1){
#   
#   if(!any(is.na(climate_direct)) & length(climate_direct) > 1){
#     climate <- climate_direct # subbing in climate manually
#   } else{
#     # generating climate
#     climate <- (clim_assym-clim_0)/(1 + exp(-clim_slope*(time-median(time)))) + clim_0
#     if(!any(is.na(noise_vec))){
#       climate <- climate + noise_vec
#     } else{
#       climate <- climate + rnorm(length(time), 0, clim_noise_sd)
#     }
#   }
#   
#   ecosystem <- acclimation.func(time, climate, lambda=lambda, E_0=E0)
#   response <- response.func(climate, ecosystem, gamma=gamma, alpha, beta, disequil_abs = disequil_abs, form=form)
#   
#   if(plot){
#     par(mar=c(3,3,1,3))
#     plot(climate, type="l", lwd=2, xlab="", ylab="", xaxt="n", yaxt="n")
#     lines(ecosystem, col="red", lwd=2)
#     par(new=T)
#     plot(response, type="l", lty=2, lwd=2, axes=F, xlab="", ylab="", col="blue")
#     legend("topleft", inset=0.05, legend=c("Climate", "Community climate niche", "Ecological response"),
#            lty=c(1,1,2), col=c("black", "red", "blue"), lwd=2, cex=0.75)
#     title(ylab="Climate", line=1, cex.lab=1.4)
#     title(xlab="Time", line=1, cex.lab=1.4)
#     mtext("Response", side=4, line=1, cex=1.4)
#     
#   }
#   
#   if(what_return == "climate") return(climate)
#   if(what_return == "ecosystem") return(ecosystem)
#   if(what_return == "response") return(response)
#   if(what_return == "ce") return(data.frame(climate = climate, ecosystem = ecosystem))
#   if(what_return == "all") return(list(climate = climate, ecosystem = ecosystem, response = response))
# }

acclimation.response.func <- function(time, clim_assym, clim_slope, clim_0, E0,
                                      beta, alpha,
                                      clim_noise_sd=0, noise_vec = NA,
                                      lambda, gamma, kappa=0, plot=F, show_range=F,
                                      what_return = "response", climate_direct=NA,
                                      disequil_abs=FALSE, disequil_sqr=FALSE, form=1,
                                      legend_pos = "topleft", legend_cex = 0.75){
  
  if(!any(is.na(climate_direct)) & length(climate_direct) > 1){
    climate <- climate_direct # subbing in climate manually
  } else{
    # generating climate
    climate <- (clim_assym-clim_0)/(1 + exp(-clim_slope*(time-median(time)))) + clim_0
    if(!any(is.na(noise_vec))){
      climate <- climate + noise_vec
    } else{
      climate <- climate + rnorm(length(time), 0, clim_noise_sd)
    }
  }
  
  ecosystem <- acclimation.func(time, climate, lambda=lambda, kappa = kappa, E_0=E0)
  response <- response.func(climate, ecosystem, gamma=gamma, alpha, beta, disequil_abs = disequil_abs, disequil_sqr = disequil_sqr, form=form)
  if(plot){
    par(mar=c(3,4,1,4))
    plot(climate, type="l", lwd=2, xlab="", ylab="", xaxt="n", yaxt="n")
    lines(ecosystem, col="red", lwd=2)
    par(new=T)
    plot(response, type="l", lty=2, lwd=2, axes=F, xlab="", ylab="", col="blue")
    legend(legend_pos, inset=0.05, legend=c("Climate", "Community climate niche", "Ecological response"),
           lty=c(1,1,2), col=c("black", "red", "blue"), lwd=2, cex=legend_cex)
    title(ylab="Climate", line=1, cex.lab=1.4)
    title(xlab="Time", line=1, cex.lab=1.4)
    mtext("Response", side=4, line=1, cex=1.4)
    
    if(show_range){
      resp_range <- range(response)
      mtext(paste(round(resp_range[1]),"  "), side = 1, adj=1, line = -1, cex=0.9, col="blue")
      mtext(paste(round(resp_range[2]),"  "), side = 3, adj=1, line = -1, cex=0.9, col="blue")
      clim_range <- range(climate)
      mtext(paste("  ", round(clim_range[1])), side = 1, adj=0, line = -1, cex=0.9, col="black")
      mtext(paste("  ", round(clim_range[2])), side = 3, adj=0, line = -1, cex=0.9, col="black")
    }
  }
  if(what_return == "climate") return(climate)
  if(what_return == "ecosystem") return(ecosystem)
  if(what_return == "response") return(response)
  if(what_return == "ce") return(data.frame(climate = climate, ecosystem = ecosystem))
  if(what_return == "all") return(list(climate = climate, ecosystem = ecosystem, response = response))
}


# Fitting by simulating latent E variable
# no algebra needed
# fit.acc.model.sim <- function(F_t, C_t, t,
#                               l_guess = 50, g_guess = -0.05, a_guess=NA, b_guess=NA, E0_guess = NA,
#                               method = "Nelder-Mead", lower, upper,
#                               predict=T, plot=T,
#                               disequil_abs = FALSE, disequil_sqr=FALSE, form=1,
#                               stepwise=FALSE){
#   
#   # catching errors
#   if(method %!in% c("Nelder-Mead", "L-BFGS-B")) stop("optimization algorithm not recognized")
#   if(length(F_t) != length(C_t)) stop("Function and Climate timeseries of unequal length")
#   if(l_guess < 0) stop("no such thing as negative time-scale")
#   if(t[1] != 1 & predict) stop("if you want a plot or prediction, times should start at 1")
#   if(any(!is.integer(t)) & predict) stop("if you want a prediction, times should be integers")
#   
#   # making initial guesses about parameters from data
#   if(is.na(a_guess)) a_guess <- mean(F_t) # baseline function
#   if(is.na(b_guess))  b_guess <- mean(diff(F_t))/ mean(diff(C_t)) # baseline climate sensitivity
#   if(is.na(E0_guess)) E0_guess <- mean(C_t[1:5]) # E initial condition
#   
#   # shifting time forward 1 for autoregression
#   C_t_moved <- C_t[-1]
#   C_tmin1 <- C_t[-length(C_t)]
#   F_t_moved <- F_t[-1]
#   F_tmin1 <- F_t[-length(F_t)]
#   
#   # cost function for parameter estimation
#   cost.func <- function(pars){
#     
#     lambda <- pars[1]
#     gamma <- pars[2]
#     alpha <- pars[3]
#     beta <- pars[4]
#     E0 <- pars[5]
#     
#     E_vec <- acclimation.func(t, C_t, lambda, E0)
#     F_pred <- response.func(C_t, E_vec,
#                             gamma=gamma, alpha=alpha, beta=beta,
#                             disequil_abs = disequil_abs,
#                             disequil_sqr = disequil_sqr,
#                             form=form)
#     
#     cost <- sum((F_t - F_pred)^2)
#     
#     # for maximum likelihood... not 100% confident that I've got the formula right
#     # n <- length(F_t[-1])
#     # log_lik <- -(n/2)*log(2*pi*sigma^2) - (1/(2*sigma^2))*sum((F_t[-1] - F_pred)^2)
#     # cost <- -log_lik
#     
#     return(cost)
#   }
#   
#   # minimizing the cost function to estimate parameters
#   if(method == "Nelder-Mead") par_est <- optim(c(l_guess, g_guess, a_guess, b_guess, E0_guess), cost.func)
#   if(method == "L-BFGS-B"){
#     if(length(lower) != 5) stop("give parameter bounds if you use the L-BFGS-B algorithm")
#     par_est <- optim(c(l_guess, g_guess, a_guess, b_guess, E0_guess), cost.func, method = "L-BFGS-B",
#                      lower = lower, upper = upper) # to give reasonable bounds
#   }
#   pars <- par_est$par
#   l_est <- pars[1]
#   g_est <- pars[2]
#   a_est <- pars[3]
#   b_est <- pars[4]
#   E0_est <- pars[5]
#   
#   F_pred <- NA; E_pred <- NA; disequil_pred <- NA # NAs returned if predict isn't requested
#   if(predict){
#     # predicting based on parameter estimates
#     sim_time <- seq(head(t, 1), tail(t, 1))
#     sim <- acclimation.response.func(time = sim_time, clim_assym = NA, clim_slope = NA,
#                                      clim_0 = NA, E0 = E0_est,
#                                      beta = b_est, alpha = a_est,
#                                      lambda = l_est, gamma = g_est,
#                                      climate_direct = C_t,
#                                      plot=F, what_return = "all", disequil_abs = disequil_abs)
#     E_pred <- sim$ecosystem
#     F_pred <- sim$response
#     disequil_pred <- C_t - E_pred
#     RSS <- sum((F_t - F_pred)^2)
#     TSS <- sum((F_t - mean(F_t))^2)
#     R2 <- 1-(RSS/TSS)
#   }
#   
#   if(plot){
#     plot(C_t ~ t, type="o", lwd=2, pch=20, xlab="", ylab="", xaxt="n", yaxt="n")
#     lines(E_pred ~ t, col="red", lwd=2)
#     par(new=T)
#     plot(F_pred ~ t, type="l", lwd=2, axes=F, xlab="", ylab="", col="blue", ylim=range(c(F_pred, F_t)))
#     points(F_t ~ t, pch=20, cex=1.5, col="blue")
#     legend("topleft", inset=0.05, legend=c("Observed climate (C_t)", "Predicted community climate (E_pred)", "Predicted ecosystem function (F_pred)", "Observed ecosystem function (F_t)"),
#            lty=c(1,1,1,NA), pch=c(20,NA,NA,20), col=c("black", "red", "blue", "blue"), lwd=2, cex=0.75)
#     #legend("topleft", inset=0.05, legend=c("Climate", "Community climate"),
#     #       lty=c(1,1), col=c("black", "red"), lwd=2, cex=0.75)
#     mtext("Climate", side=2, line=1, col="red")
#     mtext("Time", side=1, line=1)
#     mtext("Function", side=4, line=1, col="blue")
#   }
#   
#   return(list(par = c(lambda = pars[1], gamma = pars[2], alpha = pars[3], beta = pars[4], E0 = pars[5]),
#               F_pred = F_pred, E_pred = E_pred, disequil_pred = disequil_pred,
#               R2 = R2))
#   
# }

fit.acc.model.sim <- function(F_t, C_t, t,
                              l_guess = 50, g_guess = -0.05, a_guess=NA, b_guess=NA, E0_guess = NA,
                              method = "Nelder-Mead", lower, upper,
                              predict=T, plot=T, legend=T, form=1,
                              disequil_abs = FALSE, disequil_sqr = FALSE, stepwise=FALSE){
  
  # catching user errors
  if(method %!in% c("Nelder-Mead", "L-BFGS-B", "BFGS")) stop("optimization algorithm not recognized")
  if(length(F_t) != length(C_t)) stop("Function and Climate timeseries of unequal length")
  if(l_guess < 0) stop("no such thing as negative acclimation time-scale")
  if(t[1] != 1 & predict) stop("if you want a plot or prediction, times should start at 1")
  if(any(!is.integer(t)) & predict) stop("if you want a prediction, times should be integers")
  
  # making initial guesses about parameters from data
  if(is.na(a_guess)) a_guess <- mean(F_t) # baseline function
  if(is.na(b_guess))  b_guess <- mean(diff(F_t))/ mean(diff(C_t)) # baseline climate sensitivity
  if(is.na(E0_guess)) E0_guess <- mean(C_t[1:5]) # E initial condition, assuming equilibrium
  
  # shifting time forward 1 for autoregression
  # C_t_moved <- C_t[-1]
  # C_tmin1 <- C_t[-length(C_t)]
  # F_t_moved <- F_t[-1]
  # F_tmin1 <- F_t[-length(F_t)]
  
  # cost function for parameter estimation
  cost.func <- function(pars){
    
    lambda <- pars[1]
    gamma <- pars[2]
    alpha <- pars[3]
    beta <- pars[4]
    E0 <- pars[5]
    
    # simulate community niche (E) curve using lambda and E0
    E_vec <- acclimation.func(t, C_t, lambda, E0) 
    # simulate ecosystem function (F) curve using E, gamma, alpha, and beta
    F_pred <- response.func(C_t, E_vec,
                            gamma=gamma, alpha=alpha, beta=beta,
                            disequil_abs = disequil_abs, disequil_sqr = disequil_sqr,
                            form=form)
    
    cost <- sum((F_t - F_pred)^2) # sum of squares
    
    # for maximum likelihood... not 100% confident that I've got the formula right
    # n <- length(F_t[-1])
    # log_lik <- -(n/2)*log(2*pi*sigma^2) - (1/(2*sigma^2))*sum((F_t[-1] - F_pred)^2)
    # cost <- -log_lik
    
    return(cost)
  }
  
  # minimizing the cost function to estimate parameters
  if(method == "Nelder-Mead") par_est <- optim(c(l_guess, g_guess, a_guess, b_guess, E0_guess), cost.func)
  if(method == "BFGS") par_est <- optim(c(l_guess, g_guess, a_guess, b_guess, E0_guess), cost.func, method="BFGS")
  if(method == "L-BFGS-B"){
    if(length(lower) != 5) stop("give parameter bounds if you use the L-BFGS-B algorithm")
    par_est <- optim(c(l_guess, g_guess, a_guess, b_guess, E0_guess), cost.func, method = "L-BFGS-B",
                     lower = lower, upper = upper) # to give reasonable bounds
  }
  pars <- par_est$par
  l_est <- pars[1]
  g_est <- pars[2]
  a_est <- pars[3]
  b_est <- pars[4]
  E0_est <- pars[5]
  
  F_pred <- NA; E_pred <- NA; disequil_pred <- NA # NAs returned if predict isn't requested
  # predicting based on parameter estimates
  if(predict){
    sim_time <- seq(head(t, 1), tail(t, 1))
    sim <- acclimation.response.func(time = sim_time, clim_assym = NA, clim_slope = NA,
                                     clim_0 = NA, E0 = E0_est,
                                     beta = b_est, alpha = a_est,
                                     lambda = l_est, gamma = g_est,
                                     climate_direct = C_t,
                                     plot=F, what_return = "all",
                                     disequil_abs = disequil_abs, disequil_sqr = disequil_sqr,
                                     form=form)
    E_pred <- sim$ecosystem
    F_pred <- sim$response
    disequil_pred <- C_t - E_pred
    RSS <- sum((F_t - F_pred)^2)
    TSS <- sum((F_t - mean(F_t))^2)
    R2 <- 1-(RSS/TSS)
  }
  
  # plot model fit
  if(plot){
    plot(C_t ~ t, type="o", lwd=2, pch=20, xlab="", ylab="")
    lines(E_pred ~ t, col="red", lwd=2)
    par(new=T)
    plot(F_pred ~ t, type="l", lwd=2, axes=F, xlab="", ylab="", col="blue", ylim=range(c(F_pred, F_t)))
    points(F_t ~ t, pch=20, cex=1.5, col="blue")
    if(legend) legend("topleft", inset=0.05, legend=c("Observed climate (C_t)", "Predicted community climate (E_pred)", "Predicted ecosystem function (F_pred)", "Observed ecosystem function (F_t)"),
                      lty=c(1,1,1,NA), pch=c(20,NA,NA,20), col=c("black", "red", "blue", "blue"), lwd=2, cex=0.75)
    #legend("topleft", inset=0.05, legend=c("Climate", "Community climate"),
    #       lty=c(1,1), col=c("black", "red"), lwd=2, cex=0.75)
    mtext("Climate & Community niche", side=2, line=2.5, col="red")
    mtext("Time", side=1, line=2.5)
    mtext("Ecosystem function", side=4, line=2.5, col="blue")
    axis(4, at=pretty(F_t), labels=pretty(F_t))
  }
  
  return(list(par = c(lambda = pars[1], gamma = pars[2], alpha = pars[3], beta = pars[4], E0 = pars[5]),
              F_pred = F_pred, E_pred = E_pred, disequil_pred = disequil_pred,
              R2 = R2))
  
}



# runs the Alexander et al. model
sim.chalmandier <- function(L_land = 20, Tmin = 0, Tmax = 15, Tstdev = 2, deltaT = 2, d = 0.01,
                            burnin_yrs = 2000, baseline_yrs = 2000, warming_yrs = 100, final_yrs = 4000,
                            stationary_periods = T,
                            N = 40, Gmax = 0.5, Gmin = 0.2, Lmax = 0.7, Lmin = 0.7, Cmax = 0.2, Cmin = 0.2,
                            temp_niche = "cold_threshold", identical_temp_noise = F, temp_direct = NA,
                            Tmin_sp = 0, Tmax_sp = 15, dt = 0.025){

  Tmean <-seq(Tmin,Tmax,length=L_land)   # baseline mean temperatures across the landscape
  sim_yrs <- burnin_yrs + baseline_yrs + warming_yrs + final_yrs # total number of years
  
  if(length(temp_direct) == 1){
    # generate one set of temperature time series to use in all simulations
    # rows are years, cols are sites
    temperature<-matrix(NA,sim_yrs,length(Tmean))
    temp_noise <- rnorm(sim_yrs, 0, Tstdev)
    for(iT in 1:length(Tmean)){
      tmp_mean <- generateTemps(meanT=Tmean[iT],changeT=deltaT,stationary=stationary_periods, burnin_yrs, baseline_yrs, warming_yrs, final_yrs)
      
      # generate temperature time series
      # identical_temp_noise gives the same noise to temperature at all sites
      # - this is more realistic than the original in which temp variation is independent at all sites
      if(identical_temp_noise){
        temperature[,iT] = tmp_mean + temp_noise
      } else{
        temperature[,iT] = rnorm(sim_yrs,tmp_mean,Tstdev)  
      }
    }
  } else{
    tmp_mean <- NA
    temperature <- temp_direct
  }
  
  
  # set up landscape
  landscape <- cbind(1:L_land , rep(1 ,L_land),Tmean)
  
  # 1st step: initialisation of the metacommunity 
  spxp <- array(NA, dim=c(sim_yrs, L_land, N)) # [years, sites, species]
  spxp[1,,] <- 2/N
  #dt <- 0.025 ########## Figure out what this is!   
  d_ins <- d/dt # this used to be d <- d/dt. MS changed it to d_ins to avoid overwriting
  
  # Generation of the species pool
  #tr <- SpeciesPoolGen(N, Tmin, Tmax+deltaT, Gmin, Gmax, Lmin, Lmax, Cmin, Cmax, temp_niche) # original
  tr <- SpeciesPoolGen(N, Tmin_sp, Tmax_sp, Gmin, Gmax, Lmin, Lmax, Cmin, Cmax, temp_niche) # MS: species niches limited to temp range pre warming
  # tr[,1]: cold tolerance threshold (Tmin_i)
  # tr[,2]: growth rate (g_i)
  # tr[,3]: interspecific competition (l_i)
  # tr[,4]: intraspecific competition (c_i)
  # if temp_niche = "temp_optima"
  # tr[,5]: temp niche mean (mu)
  # tr[,6]: temp niche deviance (sigma)
  
  # dispersal matrix (seeds disperse only to nearest neighbor locations)
  seed_rain <- as.matrix(dist(landscape[,1:2],diag=T,upper=T))
  seed_rain[seed_rain > 1] <- 0
  seed_rain <- seed_rain*(d_ins*dt*0.5)
  diag(seed_rain) <- 1-d_ins*dt
  seed_rain[1,1] = 1-dt*d_ins/2  
  seed_rain[L_land,L_land] = 1-dt*d_ins/2
  
  # Loop through time
  for(t in 2:sim_yrs){
    yr_temp <- temperature[t,]
    spxp[t,,] <- CommunityTempDis(spxp[t-1,,], tr, seed_rain, yr_temp, dt, temp_niche)
  } 
  
  return(list(tr = tr, spxp = spxp, temp = temperature, temp_mean = tmp_mean, landscape = landscape, tr = tr,
              args = list(L_land = L_land, Tmin = Tmin, Tmax = Tmax, Tstdev = Tstdev, deltaT = deltaT, burnin_yrs = burnin_yrs, baseline_yrs = baseline_yrs, warming_yrs = warming_yrs, final_yrs = final_yrs, stationary_periods = stationary_periods, N = N, Gmax = Gmax, Gmin = Gmin, Lmax = Lmax, Lmin = Lmin, Cmax = Cmax, Cmin = Cmin, d = d, sim_yrs = sim_yrs, T_sites = Tmean)))
}

