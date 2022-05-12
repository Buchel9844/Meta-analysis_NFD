
#==================================================================================================
# GODOY'S (2014) DEFINITION - APM
#==================================================================================================

#First compute the demographic ratio
Data_Full$Dem.Ratio <- (((Data_Full$AP_lambda_j*Data_Full$AP_g_j)/(1-(1-Data_Full$AP_g_j)*Data_Full$AP_s_j)) - 1) / 
  (((Data_Full$AP_lambda_i*Data_Full$AP_g_i)/(1-(1-Data_Full$AP_g_i)*Data_Full$AP_s_i)) - 1)



Data_Full$Godoy.ND <- 1 - sqrt((Data_Full$AP_a_ij * Data_Full$AP_a_ji) / (Data_Full$AP_a_ii * Data_Full$AP_a_jj)) 
Data_Full$Godoy.FDi <- Data_Full$Dem.Ratio * (sqrt((Data_Full$AP_a_jj * Data_Full$AP_a_ji) / (Data_Full$AP_a_ii * Data_Full$AP_a_ij)))



#==================================================================================================
# CARROLL'S (2011) DEFINITION - NO MODEL
#==================================================================================================

#----- Compute the monoculture equilibrium density of the res_ident (Nj) 
#----- and the invas_ion growth rate (ri) ----
#-1-for a LV model
Data_Full$Nj_Stars <- (1 / Data_Full$LV_a_jj)
Data_Full$ri_Stars <- (1 - (Data_Full$LV_a_ij / Data_Full$LV_a_jj))
Data_Full$Ni_Stars <- 1 / Data_Full$LV_a_ii
Data_Full$rj_Stars <- 1 - (Data_Full$LV_a_ji / Data_Full$LV_a_ii)
Data_Full$f0i <- 1
Data_Full$f0j <- 1

#-2-for a annual plant  model (Godoy)
index<-Data_Full$Model.use == 'APM'
Data_Full$Nj_Stars[index] <- ((Data_Full$AP_lambda_j - 1) / Data_Full$AP_a_jj)[index]
Data_Full$ri_Stars[index] <- (log(Data_Full$AP_lambda_i / ( 1 + (Data_Full$AP_a_ij / Data_Full$AP_a_jj) * (Data_Full$AP_lambda_j - 1))))[index]
Data_Full$Ni_Stars[index] <- ((Data_Full$AP_lambda_i - 1) / Data_Full$AP_a_ii)[index]
Data_Full$rj_Stars[index] <- (log(Data_Full$AP_lambda_j / (1 + (Data_Full$AP_a_ji / Data_Full$AP_a_ii) * (Data_Full$AP_lambda_i - 1))))[index]
Data_Full$f0i[index] <- (log(Data_Full$AP_lambda_i) * 1)[index]
Data_Full$f0j[index] <- (log(Data_Full$AP_lambda_j) * 1)[index]
sum(index)

#-3- for a annual plant  model 2 (Bimler)
index1 <- Data_Full$Model.use == 'AP2'
Data_Full$Nj_Stars[index1] <- ((Data_Full$AP_lambda_j - 1) / Data_Full$AP_a_jj_exp)[index1]
Data_Full$ri_Stars[index1] <- (log(Data_Full$AP_lambda_i / (exp(Data_Full$AP_a_ij_exp / Data_Full$AP_a_jj_exp) * log(Data_Full$AP_lambda_j))))[index1]
Data_Full$Ni_Stars[index1] <- ((Data_Full$AP_lambda_i - 1) / Data_Full$AP_a_ii_exp)[index1]
Data_Full$rj_Stars[index1] <- (log(Data_Full$AP_lambda_j / (exp(Data_Full$AP_a_ji_exp / Data_Full$AP_a_ii_exp) * log(Data_Full$AP_lambda_i))))[index1]
Data_Full$f0i[index1] <- (log(Data_Full$AP_lambda_i) * 1)[index1]
Data_Full$f0j[index1] <- (log(Data_Full$AP_lambda_j) * 1)[index1]



#----- Compute Sens_itivity -----
Data_Full$Si_Stars <- (Data_Full$f0i - Data_Full$ri_Stars) / Data_Full$f0i
Data_Full$Sj_Stars <- (Data_Full$f0j - Data_Full$rj_Stars) / Data_Full$f0j

index2 <- Data_Full$Initial.Definition == 'Carroll'
Data_Full$Si_Stars[index2] <- ((Data_Full$g_ii - Data_Full$g_ij) / Data_Full$g_ii)[index2]
Data_Full$Sj_Stars[index2] <- ((Data_Full$g_jj - Data_Full$g_ji) / Data_Full$g_jj)[index2]

#### Warning message:In Ops.factor(Data_Full$g_ii, Data_Full$g_ij) : '-' not meaningful for factors
# the s_i and s_j produce with the above method are not the same than the one
# in the papers
#Data_Full$Si_Stars[index2] <- Data_Full$s_i[index2]
#Data_Full$Sj_Stars[index2] <- Data_Full$s_j[index2]




# Computation of the sens_itivities from the NDF, for Verseglou , Li  and Tan 
index8 <- (Data_Full$Initial.Definition =="Carroll" & (Data_Full$First.Author == 'Li'|Data_Full$First.Author == 'Tan'))

Data_Full$Si_Stars[index8] <- ((1 - Data_Full$ND)* Data_Full$FD)[index8]
Data_Full$Sj_Stars[index8] <- ((1 - Data_Full$ND)/ Data_Full$FD)[index8]


# to be able to calculate the geometric mean and standard deviation, 
#there is a need to create a matrix to inverse the column and the row

S<-data.frame(Data_Full$Sj_Stars,Data_Full$Si_Stars)
tS<-as.data.frame(t(S))
x<-data.matrix(tS,rownames.force = NA)

#----- Computing the geometric mean and the standard geometric deviation ----

S['M'] <- apply(x, 2, geometric.mean)

#S['V'] <- exp(apply(log(x), 2, var)/2)
# computation of NFD
Data_Full['Carroll.ND'] <- (1-(S['M']))
#Data_Full['Carroll.ND'] <- 1 - sqrt(Data_Full$Sj_Stars*Data_Full$Si_Stars)

#Data_Full['Carroll.FD'] <- 1-(1/(S['V']))
Data_Full$Carroll.FDi <- NA
for(i in 1:nrow(Data_Full)){
  if (is.na(Data_Full$Si_Stars[i])|is.na(Data_Full$Sj_Stars[i])) next
if (Data_Full$Sj_Stars[i] < Data_Full$Si_Stars[i]){
  Data_Full$Carroll.FDi[i] <- 1 - 1/sqrt(Data_Full$Si_Stars[i]/Data_Full$Sj_Stars[i])
 }
  else{
    Data_Full$Carroll.FDi[i] <- 1- 1/sqrt(Data_Full$Sj_Stars[i]/Data_Full$Si_Stars[i])
  }
}


#----- Copy the NFD of paper us_ing Carrol's def ----
# all the paper us_ing carroll's definition g_ive their NFD but don't always g_ive the sens_itivity matrix
#Data_Full$Carroll.ND [Data_Full$Model.use == "none" ] <- Data_Full$ND [Data_Full$Model.use == "none"]
#Data_Full$Carroll.FD [Data_Full$Model.use == "none" ] <- Data_Full$FD [Data_Full$Model.use == "none"]


#==================================================================================================
# Zhao'S (2016) DEFINITION - NO MODEL
#==================================================================================================



Data_Full$Zhao.ND  <- ((Data_Full$rj_Stars + Data_Full$ri_Stars)/2) 
#Data_Full$Zhao.FDi <- log10( Data_Full$Ni_Stars / Data_Full$Nj_Stars)
#Data_Full$Zhao.FDj <- log10( Data_Full$Nj_Stars / Data_Full$Ni_Stars)

n <- data.frame(x=log10( Data_Full$Ni_Stars / Data_Full$Nj_Stars) ,y=(log10( Data_Full$Ni_Stars / Data_Full$Nj_Stars))) 
Data_Full$Zhao.FDi <- apply(n , 1, max)

index5 <- Data_Full$First.Author== 'Zhao'
Data_Full$Zhao.ND [index5]  <- (((Data_Full$ND_Zhao)-1)/2) [index5] 
Data_Full$Zhao.FDi [index5]  <- Data_Full$FD_Zhao [index5] 

index6<-Data_Full$Initial.Definition== 'Carroll'
Data_Full$Zhao.ND [index6] <- ((Data_Full$g_ij + Data_Full$g_ji)/2) [index6]

#==================================================================================================
# Chesson'S (2000) DEFINITION - LOTKA VOLTERRA MODEL
#==================================================================================================

Data_Full$Chesson00.ND <- 1 - sqrt((Data_Full$LV_a_ij * Data_Full$LV_a_ji) / (Data_Full$LV_a_ii * Data_Full$LV_a_jj)) 
Data_Full$Chesson00.FDi <- sqrt( (Data_Full$LV_a_jj * Data_Full$LV_a_ji) / (Data_Full$LV_a_ii * Data_Full$LV_a_ij) ) 


#==================================================================================================
# BIMLER'S (2018) DEFINITION - APM 2
#==================================================================================================

Data_Full$Bimler.ND<- 1 - ( (exp(Data_Full$AP_a_ij_exp) * exp(Data_Full$AP_a_ji_exp) ) / ( exp(Data_Full$AP_a_ii_exp) * exp(Data_Full$AP_a_jj_exp) ) )
Data_Full$Bimler.FDi <- ( exp(Data_Full$AP_a_jj_exp) * exp(Data_Full$AP_a_ji_exp) ) / (exp(Data_Full$AP_a_ii_exp) * exp(Data_Full$AP_a_ij_exp) ) 


#==================================================================================================
# ADLER'S (2007) DEFINITION - APM
#==================================================================================================

Data_Full$Adler.ND <- log( Data_Full$AP_lambda_j / ( 1 + (Data_Full$AP_a_ij / Data_Full$AP_a_jj) * (Data_Full$AP_lambda_j - 1)))
Data_Full$Adler.FDi <- log( Data_Full$AP_lambda_j / Data_Full$AP_lambda_i )
index5<-Data_Full$Initial.Definition == 'Bimler'
Data_Full$Adler.FDi [index5] <-(1 * 0)[index5]
index4<-Data_Full$Initial.Definition == 'Chesson 00'
Data_Full$Adler.FDi [index4] <-(1 * 0)[index4]


