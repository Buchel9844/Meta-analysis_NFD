# The script is originally written in python, we run python code from within R
# make sure that python is installed on the system used
# This is not part of the install.packages("reticulate") command!
library(reticulate)
np <- import("numpy",convert=F)
# loads the relevant python code
source_python("Code/numerical_NFD.py")



# create the differential equation system
n_spec <- 2 # number of species in the system, must be an integer



# Lotka Volterra model
new_F <- function(N,mu,A){
  return(mu - A%*%N)
}


LV_NDF_Jurg <- matrix(NA,nrow(Data),8)
LV_NDF_Jurg <- as.data.frame(LV_NDF_Jurg)
names(LV_NDF_Jurg ) <- c("NDi","NDj","NOi",
                         "NOj","FDi","FDj",
                         "ci","cj")


for (i in 1:nrow(Data)){ # start of the loop
  
  
  mu <- c(data=c( 1, 1)) # intrinsic growth rate
  A <- matrix(data=c(Data$LV_a_ii[i],Data$LV_a_ji[i],
                     Data$LV_a_ij[i],Data$LV_a_jj[i]),
              nrow=n_spec, ncol=n_spec) # interaction matrix
  
  if (any(is.na(c(mu,A)))) next  # When there is a value equal to NA, the computation of NDF for that line is stopped
  if (any(A[1,1]<=0)) next
  if (any(A[2,2]<=0)) next
  
  
  # compute relevant parameters with python
  # the parameter `from_R = TRUE` changes data types from R to python
  pars <- NFD_model(new_F,n_spec,args=list(mu,A), from_R = TRUE)
  
  
  LV_NDF_Jurg$NDi[i] <- pars$ND[1]
  LV_NDF_Jurg$NDj[i] <- pars$ND[2]
  
  
  LV_NDF_Jurg$NOi[i] <- pars$NO[1]
  LV_NDF_Jurg$NOj[i] <- pars$NO[2]
  
  LV_NDF_Jurg$FDi[i] <- pars$FD[1]
  LV_NDF_Jurg$FDj[i] <- pars$FD[2]
  
  LV_NDF_Jurg$ci[i] <- pars$c[3]
  LV_NDF_Jurg$cj[i] <- pars$c[2]
  
} # end of the loop 
