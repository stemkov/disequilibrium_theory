### demonstrating climate:function relationship at different disequilibrium values

alpha <- 1000
beta <- 10
gamma <- .5
n <- 300
form <- 2
disequil_abs <- FALSE
disequil_sqr <- TRUE
C <- seq(75,125,length.out=n)
Es <- seq(80,120, by=20)
Disequils <- Es - 100
F_nodisequil <- response.func(climate = C, ecosystem = C, gamma=gamma, beta=beta, alpha=alpha, form=form, disequil_abs = disequil_abs, disequil_sqr = disequil_sqr)

cols <- colorRampPalette(c("blue","red"))(length(Es))
y_range <- range(c(response.func(climate = C, ecosystem = Es, gamma=gamma, beta=beta, alpha=alpha, form=form, disequil_abs = disequil_abs, disequil_sqr = disequil_sqr),
                   F_nodisequil <- response.func(climate = C, ecosystem = C, gamma=gamma, beta=beta, alpha=alpha, form=form, disequil_abs = disequil_abs, disequil_sqr = disequil_sqr)))

# Figure 2d
svg("figures/disequil_sensitivities_demo_sqr.svg", width=6, height=6)

  plot(NA, xlim=range(C), ylim=y_range,
       ylab="Ecosystem function", xlab="Interannual climate", yaxt="n",
       main = str2expression("F == alpha + beta*C + gamma*(C-E)*C")) # could automate ylim
  
  for(i in seq_along(Es)){
    E <- rep(Es[i], n)
    F_calc <- response.func(climate = C, ecosystem = E, gamma=gamma, beta=beta, alpha=alpha, form=form, disequil_abs = disequil_abs, disequil_sqr = disequil_sqr)
    lines(F_calc ~ C, col = cols[i], lwd=3)
    C_equil <- C[which.min(abs(C - E))]
    F_equil <- F_calc[which.min(abs(C - E))]
    if(which.min(abs(C - E)) %in% c(1,length(E))) F_equil <- NA # check if it's the first or last
    points(F_equil ~ C_equil, pch=20,cex=2, col=cols[i])
  }
  legend("topleft", inset=0.05, legend = Disequils, title="Disequilibrium", lwd=3, col=cols, cex=0.9, bty="n")
  abline(v = 100, lty=2, lwd=1)
  lines(F_nodisequil ~ C, lwd=1)

dev.off()



### icons for each functional form - for Table S1

sens.plot <- function(form, alpha, beta, gamma, E1, E2, C_mean, C_range, col1, col2, equation, disequil_abs, disequil_sqr, n=300, save_out=F, mar_vec=c(1,1,2.5,1), Fmax=NA){
  
  if(FALSE){
    form = 1
    alpha = 1000; beta = 10; gamma = 0.2;
    E1 = 90; E2 = 110; C_mean = 100; C_range = 25
    col1 = "blue"; col2 = "red";
    equation = equations[2]; disequil_abs = F; disequil_sqr=T
  }
  
  if(save_out) png(paste("figures/icon_",form,"_2.png",sep=""), width=200, height=140)
  par(mar=mar_vec)
  # print(C_mean)
  # print(C_range)
  # print(n)
  C <- seq(C_mean - C_range, C_mean + C_range, length.out=n)
  Es <- c(E1, E2)
  Disequils <- Es - 100
  F_nodisequil <- response.func(climate = C, ecosystem = C, gamma=gamma, beta=beta, alpha=alpha, form=form, disequil_abs = disequil_abs, disequil_sqr = disequil_sqr, Fmax=Fmax)
  
  cols <- c(col1, col2)
  y_range <- range(c(response.func(climate = C, ecosystem = Es, gamma=gamma, beta=beta, alpha=alpha, form=form, disequil_abs = disequil_abs, disequil_sqr = disequil_sqr, Fmax=Fmax),
                     F_nodisequil))
  
  plot(NA, xlim=range(C), ylim=y_range,
       ylab="", xlab="", yaxt="n", xaxt="n")
  
  for(i in seq_along(Es)){
    E <- rep(Es[i], n)
    F_calc <- response.func(climate = C, ecosystem = E, gamma=gamma, beta=beta, alpha=alpha, form=form, disequil_abs = disequil_abs, disequil_sqr = disequil_sqr, Fmax=Fmax)
    lines(F_calc ~ C, col = cols[i], lwd=4)
    C_equil <- C[which.min(abs(C - E))]
    F_equil <- F_calc[which.min(abs(C - E))]
    if(which.min(abs(C - E)) %in% c(1,length(E))) F_equil <- NA # check if it's the first or last
    points(F_equil ~ C_equil, pch=20,cex=2, col=cols[i])
  }
  abline(v = 100, lty=2, lwd=2)
  lines(F_nodisequil ~ C, lwd=2)
  mtext("C ", side=1, line=.5, adj=1, cex=1.5) # C and F outside box
  mtext("F", side=3, line=-1.5, at = -20, cex=1.5)
  mtext("0", side=1, line=.45, at=-10, cex=1.5)
  if(equation != "") mtext(str2expression(equation),side=3, line=0.45)
  
  if(save_out) dev.off()
}

equations <- c("F == alpha + beta*C",
               "F == alpha + beta*C + gamma*(C-E)^2*C",
               "F == alpha + beta*C + gamma*(C-E)^2",
               "F == alpha + gamma*(C-E)^2",
               "F == alpha + gamma*(C-E)^2*C",
               "F == alpha*E + gamma*(C-E)^2",
               "F == alpha*E + beta*C + gamma*(C-E)^2*C",
               "F == alpha + (beta*C)/(gamma*(C-E)^2)",
               "(alpha + beta*C) * e^(-(C-E)^2/gamma)^2")

alphas <- rep(1000,9)
gammas <- rep(0.001,9); gammas[3] <- 0.2; gammas[6] <- 10; gammas[7] <- .05; gammas[8] <- 1; gammas[9] <- 2000

# adding nonlinear equilibrium forms
nonlin_equations <- c("F == alpha + kappa * (1 - e^(-beta*C)) - gamma*(C-E)^2",
                      "F == alpha + kappa * (1 - e^(-beta * C + gamma * (C-E)^2 * C))",
                      "F == alpha + kappa * (1 - e^(-beta * C + gamma * (C-E)^2))")
nonlin_alphas <- c(1000,1000,1000)
nonlin_gammas <- c(0.1, 0.00000125, 0.00005)
nonlin_betas <- c(0.02, 0.02, 0.02)
nonlin_Fmaxs <- rep(1500, 3)

equations <- c(equations, nonlin_equations)
gammas <- c(gammas, nonlin_gammas)
alphas <- c(alphas, nonlin_alphas)
betas <- c(rep(10, 9), nonlin_betas)
Fmaxs <- c(rep(NA, 9), nonlin_Fmaxs)

forms <- c(0:7, 9, 1003, 1001, 1002)

forms_table <- data.table(form = forms,
                          alpha = alphas, beta = betas, gamma = gammas,
                          E1 = 75, E2 = 125, C_mean = 100, C_range = 100,
                          col1 = "blue", col2 = "red",
                          equation = equations, disequil_abs = F, disequil_sqr=T,
                          Fmax = Fmaxs)

par(mar=c(1,1,2.5,1),mfrow=c(3,4)) # set save_out to FALSE to plot all icons here

forms_table[, sens.plot(form = form, alpha = alpha, beta = beta, gamma = gamma,
                        E1 = E1, E2 = E2, C_mean = C_mean, C_range = C_range,
                        col1 = col1, col2 = col2, equation = equation,
                        disequil_abs = disequil_abs, disequil_sqr = disequil_sqr,
                        Fmax = Fmax,
                        mar_vec=c(1.5,1.5,0.25,0.25),
                        save_out = TRUE),
            by=.(form)]
