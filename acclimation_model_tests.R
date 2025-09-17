###### Parameter identifiability analysis

test.acc.model <- function(l_given, g_given, a_given, b_given, E0_given, lower=NA, upper=NA,
                           disequil_abs = FALSE, disequil_sqr = FALSE,
                           obs_error_sd = 0, obs_error_mult = 0,
                           form = 1, method="BFGS", Fmax=NA){
  
  if(FALSE){
    l_given <- lambda_given[1]
    g_given <- gamma_given[1]
    a_given <- alpha_given[1]
    b_given <- beta_given[1]
    E0_given <- E0_given[1]
    lower=NA
    upper=NA
    disequil_abs = FALSE
    disequil_sqr = TRUE
    obs_error_sd = 0
    obs_error_mult = 0.1
    form = 1
    method="BFGS"
  }
  
  time <- 1:300
  sim <- acclimation.response.func(time = time, clim_assym = 200, clim_slope = .05,
                                   clim_0 = 100, E0 = E0_given,
                                   beta = b_given, alpha = a_given,
                                   lambda = l_given, gamma = g_given,
                                   clim_noise_sd = 0, noise_vec = rnorm(300, 0, 2),
                                   plot = F, what_return = "all",
                                   disequil_abs = disequil_abs, disequil_sqr = disequil_sqr,
                                   form = form, Fmax=Fmax)
  F_sd <- sd(diff(sim$response))
  F_t <- sim$response + rnorm(300, 0, obs_error_sd) + rnorm(300, 0, obs_error_mult*F_sd)
  C_t <- sim$climate
  
  if(length(lower) != 1){
    fit <- fit.acc.model.sim(F_t = F_t, C_t = C_t, t = time, plot=F, disequil_abs = disequil_abs, disequil_sqr = disequil_sqr,
                             method = "L-BFGS-B", lower = lower, upper = upper,
                             form = form, Fmax=Fmax)
  } else{
    fit <- fit.acc.model.sim(F_t = F_t, C_t = C_t, t = time, plot=F, disequil_abs = disequil_abs, disequil_sqr = disequil_sqr,
                             form = form, method=method, Fmax=Fmax)
  }
  
  
  l_est <- as.numeric(fit$par[1])
  g_est <- as.numeric(fit$par[2])
  a_est <- as.numeric(fit$par[3])
  b_est <- as.numeric(fit$par[4])
  E0_est <- as.numeric(fit$par[5])
  
  return(c(lambda_given = l_given, lambda_est = l_est,
           gamma_given = g_given, gamma_est = g_est,
           alpha_given = a_given, alpha_est = a_est,
           beta_given = b_given, beta_est = b_est,
           E0_given = E0_given, E0_est = E0_est,
           R2 = fit$R2))
}

set.seed(0)

n_scans <- 50
lambda_given <- runif(n_scans, 5, 100)
gamma_given <- runif(n_scans, -0.1, 0.1)
alpha_given <- runif(n_scans, 500, 1500)
beta_given <- runif(n_scans, -5, 5)
E0_given <- runif(n_scans, 50, 150)

varying_all <- as.data.frame(t(mapply(test.acc.model,
                                      l_given = lambda_given,
                                      g_given = gamma_given,
                                      a_given = alpha_given,
                                      b_given = beta_given,
                                      E0_given = E0_given,
                                      MoreArgs = list(disequil_abs = FALSE,
                                                      disequil_sqr = TRUE,
                                                      method = "BFGS",
                                                      obs_error_sd = 0,
                                                      obs_error_mult = .1,
                                                      form = 1))))
# Figure S8
par(mfrow=c(3,2), mar=c(4,4,1,1))
plot(lambda_est ~ lambda_given, data=varying_all,
     ylim=c(0,100))
abline(0,1,col="red")
plot(gamma_est ~ gamma_given, data=varying_all,
     ylim=c(-0.1,0.1))
abline(0,1,col="red")
plot(alpha_est ~ alpha_given, data=varying_all,
     ylim=c(500,1500))
abline(0,1,col="red")
plot(beta_est ~ beta_given, data=varying_all,
     ylim=c(-5,5))
abline(0,1,col="red")
plot(E0_est ~ E0_given, data=varying_all,
     ylim=c(50,150))
abline(0,1,col="red")
hist(varying_all$R2, xlim=c(0,1), main="R squared")


# fit to model data figure

set.seed(0)

beta <- 1
alpha <- 1000
gamma <- 0.002
lambda <- 75
E0 <- 120
time <- c(1:300)

sim <- acclimation.response.func(time = time, clim_assym = 200, clim_slope = .05,
                                   clim_0 = 100, E0 = E0,
                                   beta = beta, alpha = alpha,
                                   lambda = lambda, gamma = gamma,
                                   clim_noise_sd = 0, noise_vec = rnorm(300, 0, 0),
                                   plot=F, what_return = "all", show_range = T,
                                   disequil_abs = FALSE, disequil_sqr = TRUE, form=1)

F_t <- sim$response + rnorm(300, 0, 100)
C_t <- sim$climate
E_t <- sim$ecosystem

# figure 7
png("figures/fit_sim_data.png", width=800, height=500, pointsize = 15)

  layout(matrix(c(1,1,1,1,1,1,2,3,4),nrow=3,ncol=3))
  
  par(mar=c(4.2,4.2,2,4.2))
  fit <- fit.acc.model.sim(F_t = F_t, C_t = C_t, t = time, disequil_abs = FALSE, disequil_sqr = TRUE,
                           method="BFGS", form=1, legend=F)
  legend("left", inset=0.05, legend=c("Climate", "Predicted E", "Predicted F", "Observed F"),
         lty=c(1,1,1,NA), pch=c(20,NA,NA,20), col=c("black", "red", "blue", "blue"), lwd=2, cex=1)
  plot(fit$disequil_pred ~ time, ylab="Inferred disequilibrium", xlab="Time", type="l", lwd=3); abline(h=0)
  plot(fit$F_pred ~ F_t, xlab="Observed F", ylab="Predicted F"); abline(0,1,col="red", lwd=2)
  plot(fit$E_pred ~ E_t, xlab="Observed E", ylab="Predicted E"); abline(0,1,col="red", lwd=2)

dev.off()


### lambda vs. dC/dt heatmap - Figure 2b

disequil.mag <- function(dC, lambda){
  climate <- 100 + dC*c(1:1000)
  sim <-  acclimation.response.func(time = 1:1000, climate_direct = climate,
                                    clim_assym = NA, clim_slope = NA, clim_0 = NA,
                                    E0 = 100, beta = 1, alpha = 1000,
                                    lambda = lambda, gamma = 1,
                                    what_return = "ce", plot=F)
  disequil <- sim$climate - sim$ecosystem
  disequil_final <- tail(disequil, 1)
  return(disequil_final)
}

dC_vals <- seq(0,1, by=.01)
lambda_vals <- seq(1,100, by=1)
par_scans <- expand.grid(dC_vals, lambda_vals)

par_scans$disequil <- mapply(disequil.mag, dC = par_scans$Var1, lambda = par_scans$Var2)
disequil_mat <- matrix(par_scans$disequil, nrow=length(dC_vals))

svg("figures/dC_lambda_heatmap.svg", width=6.5, height=5)

  par(mar=c(4.5,4.5,1,3))
  image.plot(x = dC_vals, y = lambda_vals, z = disequil_mat,
             col=colorRampPalette(c("white","orange","red"))(10),
             xlab = "Pace of climate change (degree C/time)",
             ylab= expression(paste("Timescale community adjustment ",lambda," (time)")),
             xaxt="n", yaxt="n",
             legend.shrink = 0.75)

dev.off()