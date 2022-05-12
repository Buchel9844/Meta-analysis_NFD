# @author: J.W.Spaak
# Example how to compute the ND and FD for a given differential equation setting

# The script is originally written in python, we run python code from within R
# make sure that python is installed on the system used
# This is not part of the install.packages("reticulate") command!
library(reticulate)
np <- import("numpy",convert=F)
# loads the relevant python code
source_python("Code/numerical_NFD.py")



# create the differential equation system
n_spec <- 2 # number of species in the system, must be an integer
set.seed(0) # set random seed for reproduce ability


# Annual Plant 2 model
APM2 <- function(N,s,g,lam,A){
  #print(c(1/N , (log((1-g)*s + (g*lam*exp(A%*%(g*N)))))))
  return(log((1-g)*s + (g*lam*exp(A%*%(g*1000*N)))))
  
} # log to transform the discrete time into continuous


AP2_NDF_Jurg <- matrix(NA,nrow(Data),8)
AP2_NDF_Jurg <- as.data.frame(AP2_NDF_Jurg)
names(AP2_NDF_Jurg )<- c("NDi","NDj","NOi",
                         "NOj","FDi","FDj",
                         "ci","cj")


for (i in 1:nrow(Data)){ # start of the loop
  lam<-c(Data$AP_lambda_i_exp[i],Data$AP_lambda_j_exp[i]) # viable seedt
  s <- c(Data$AP_s_i_exp[i],Data$AP_s_j_exp[i]) #seed survivak
  g <- c(Data$AP_g_i_exp[i],Data$AP_g_j_exp[i]) #germination
  A <- matrix(data=c(Data$AP_a_ii_exp[i],Data$AP_a_ji_exp[i],Data$AP_a_ij_exp[i],Data$AP_a_jj_exp[i]),
              nrow=n_spec, ncol=n_spec) # interaction matrix
  
  if (any(is.na(c(lam,s,g,A)))) next  # When there is a value equal to NA, the computation of NDF for that line is stopped
  if (any(g*lam<=1)) next
  
  N_star <- ((log(1-(1-g)*s)/(lam*g))/A) # not needed for lotka volterra
  N_star[2,1] <- N_star[1,1]
  N_star[1,2] <- N_star[2,2] # 
  N_star <- np$array(N_star) 
  N_star$setflags(write =TRUE)
  pars0<-list(N_star= N_star) # 
  
  #print(N_star)
  # compute relevant parameters with python
  # the parameter `from_R = TRUE` changes data types from R to python
  pars <- NFD_model(APM2, n_spec,args=list(s,g,lam,A), from_R = TRUE,
                    pars=pars0)

  
  AP2_NDF_Jurg$NDi[i] <- pars$ND[1]
  AP2_NDF_Jurg$NDj[i] <- pars$ND[2]
  
  
  AP2_NDF_Jurg$NOi[i] <- pars$NO[1]
  AP2_NDF_Jurg$NOj[i] <- pars$NO[2]
  
  AP2_NDF_Jurg$FDi[i] <- pars$FD[1]
  AP2_NDF_Jurg$FDj[i] <- pars$FD[2]
  
  AP2_NDF_Jurg$ci[i] <- pars$c[3]
  AP2_NDF_Jurg$cj[i] <- pars$c[2]
  
} # end of the loop 

