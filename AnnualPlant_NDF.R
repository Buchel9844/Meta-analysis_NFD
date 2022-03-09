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


# Annual Plant model
APM <- function(N,s,g,lam,A){
 # print (c(N,s,g,lam,A))
  # print(log((1-g)*s + (g*lam)/(1+ A%*%(g*N))))
  return(log((1-g)*s + (g*lam)/(1+ A%*%(g*N))))

} # log to transform the discrete time into continuous


AP_NDF_Jurg <- matrix(NA,nrow(Data),8)
AP_NDF_Jurg <- as.data.frame(AP_NDF_Jurg)
names(AP_NDF_Jurg) <- c("NDi","NDj","NOi",
                        "NOj","FDi","FDj",
                        "ci","cj")


for (i in 1:nrow(Data)){ # start of the loop
lam<-c(Data$AP_lambda_i[i],Data$AP_lambda_j[i]) # viable seedt
s <- c(Data$AP_s_i[i],Data$AP_s_j[i]) #seed survivak
g <- c(Data$AP_g_i[i],Data$AP_g_j[i]) #germination
A <- matrix(data=c(Data$AP_a_ii[i],Data$AP_a_ji[i],Data$AP_a_ij[i],Data$AP_a_jj[i]),
            nrow=n_spec, ncol=n_spec) # interaction matrix

if (any(is.na(c(lam,s,g,A)))) next  # When there is a value equal to NA, the computation of NDF for that line is stopped
if (any(g*lam <= 1)) next



N_star <-(((lam*g)/(1-(1-g)*s)-1)/(A*g))
N_star[2,1] <- N_star[1,1] 
N_star[1,2] <- N_star[2,2] 
N_star <- np$array(N_star) 
N_star$setflags(write =TRUE)
pars0<-list(N_star= N_star) 
        

# compute relevant parameters with python
# the parameter `from_R = TRUE` changes data types from R to python
try({pars <- NFD_model(APM, n_spec,args=list(s,g,lam,A), from_R = TRUE,
                  pars=pars0)

AP_NDF_Jurg$NDi[i] <- pars$ND[1]
AP_NDF_Jurg$NDj[i] <- pars$ND[2]


AP_NDF_Jurg$NOi[i] <- pars$NO[1]
AP_NDF_Jurg$NOj[i] <- pars$NO[2]

AP_NDF_Jurg$FDi[i] <- pars$FD[1]
AP_NDF_Jurg$FDj[i] <- pars$FD[2]

AP_NDF_Jurg$ci[i] <- pars$c[3]
AP_NDF_Jurg$cj[i] <- pars$c[2]
  })

                          } # end of the loop 












