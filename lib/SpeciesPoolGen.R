#
#
# Help function to generate traits within the species pool
#	Main script author: Loïc Chalmandrier - lchalman@uwyo.edu, https://github.com/LoicChr
#	Co-authors (and project leaders): J. Alexander - jake.alexander@unil.ch & L. Pellissier - loic.pellissier@usys.ethz.ch

# Translated from Matlab to R by Peter Adler - peter.adler@usu.edu

# Modified by Michael Stemkovski to allow hump-shaped species niches in place of cold tolerance - m.stemkovski@gmail.com

SpeciesPoolGen <- function(N,Tmin, Tmax, Gmin, Gmax, Lmin, Lmax, Cmin, Cmax, temp_niche="cold_threshold"){
  
  if(temp_niche == "cold_threshold"){
    
    tr = matrix(0, N, 4) 
    
    # Competition intra
    tr[,1] <- sort(runif(N)*(Tmax-Tmin) + Tmin) # Tmin_i
    tr[,4] <- runif(N)*(Cmax-Cmin) + Cmin # c_i
    
    # Growth rate
    slope <- (Gmax-Gmin)/(Tmax-Tmin)
    tr[,2] <- tr[,1]*slope+Gmin - slope*Tmin # g_i
    
    # Biomass tolerance
    slope <- -(Lmax-Lmin)/(Tmax-Tmin)
    tr[,3] <- tr[,1]*slope+Lmax - slope*Tmin # l_i
  }
  
  if(temp_niche == "temp_optima"){
    
    tr = matrix(0, N, 6) 
    
    tr[,1] <- sort(runif(N)*(Tmax-Tmin) + Tmin) # Tmin_i, now unused
    tr[,5] <- sort(runif(N)*(Tmax-Tmin) + Tmin) # mu_i, currently same as Tmin_i
    tr[,6] <- rep((Tmax-Tmin)/10, N) # sigma_i, one tenth of whole temperature range, same for all species
    
    # Competition intra
    tr[,4] <- runif(N)*(Cmax-Cmin) + Cmin # c_i
    
    # Growth rate
    slope <- (Gmax-Gmin)/(Tmax-Tmin)
    tr[,2] <- tr[,5]*slope+Gmin - slope*Tmin # g_i # now correlated with temp optima rather than cold threshold
    
    # Biomass tolerance
    slope <- -(Lmax-Lmin)/(Tmax-Tmin)
    tr[,3] <- tr[,1]*slope+Lmax - slope*Tmin # l_i
    
  }
  
  return(tr)

}
  


