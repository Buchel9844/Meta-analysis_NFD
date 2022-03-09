#######################################################################################################-
## Pairwise coexistence analysis ----
##
## to accompany 

## Contains:
## 1.0 Compute the three definitions
## 2.0 Computation of Spaak definition for the paper using Carroll's definition

#######################################################################################################----
## 1.0 Compute the three definitions
#######################################################################################################----
# need to install xcode: in your terminal put "xcode-select --install"
install.packages("reticulate")
library(reticulate)
use_condaenv("NFD_meta") # pandas numpy scipy
# loads the relevant python code
source_python("Code/numerical_NFD.py") # function NFD_model
source_python("Code/compute_NFD_metaanalysis.py")

NFD_metaanalysis.initial  <- read.csv("Data/NFD_meta.csv")
str(NFD_metaanalysis.initial)
NFD_metaanalysis_FD <- gather(NFD_metaanalysis.initial,"FC","FZ","FS" , key = "F.definition", 
                              value = "FD", na.rm = FALSE)

NFD_metaanalysis_ND <- gather(NFD_metaanalysis.initial,"NC_i","NZ_i","NS_i" , key = "N.definition", 
                              value = "ND", na.rm = FALSE)

number.row.NFD_metaanalysis <- nrow(NFD_metaanalysis.initial)

NFD_metaanalysis <- subset(NFD_metaanalysis_ND,select=c("First.Author","Year","Initial.Definition",
                                                     "Model.use","lab.OR.field","treatment",
                                                     "Kingdom","species.i","species.j","S_i","S_j",
                                                     "r_j","r_i","comp_outcome","ND"))

NFD_metaanalysis$Definition <- rep(c("Carroll","Zhao","Spaak"),each=number.row.NFD_metaanalysis)

NFD_metaanalysis$FD <- NFD_metaanalysis_FD$FD
nrow(NFD_metaanalysis)
#######################################################################################################----
## 2.0 Computation of Spaak definition for the paper using Carroll's definition
#######################################################################################################----
#==================================================================================================
# 2.1 Creation of data_frame with the sensitivity Si and Sj
#==================================================================================================

# number of point 
npoints <- 10

# Creation of a subset of the Data_Full with the "Carroll's data"

SpaakSensitivity <- subset(NFD_metaanalysis.initial, NFD_metaanalysis.initial$Initial.Definition == "Carroll",
                           select=c("First.Author","Year","Initial.Definition",
                           "Model.use","lab.OR.field","treatment",
                           "Kingdom","species.i","species.j","S_i","S_j",
                           "r_j","r_i","comp_outcome"))
names(SpaakSensitivity)

# We need to transform the factor into character, to be able to transpose them in an other data_frame
SpaakSensitivity %>% mutate_if(is.factor, as.character) -> SpaakSensitivity

# Creation of a new dataframe with NFD according to Spaak
names.SpaakSensitivity.NFD <- c(names(SpaakSensitivity),
                                 "ND","FD",
                                 "Coex.Carroll","Coex.Spaak")
ncolumn <- length(names.SpaakSensitivity.NFD)

SpaakSensitivity.NFD <- data.frame(matrix(nrow=0, ncol = ncolumn))
names(SpaakSensitivity.NFD ) <- names.SpaakSensitivity.NFD 

#========================================================================================================
# 2.2 ForLoop which defines the range of NFD and than choose three points on the line within this range
#========================================================================================================

for (i in 1:nrow(SpaakSensitivity)) { # for each row 
  if (is.na(SpaakSensitivity$S_i[i]) |is.na(SpaakSensitivity$S_j[i]) |
      SpaakSensitivity$S_i[i]== 0 | SpaakSensitivity$S_j[i] == 0) next
  # randomly sample NO between S_i, |S_i|<|NO_i|<|S_j|
  S <- c(abs(SpaakSensitivity$S_i[i]),abs(SpaakSensitivity$S_j[i]))
  nd <-  runif(npoints, min=min(S), max=max(S))

  fi <- 1-nd/(S[1])
  fj <- 1-nd/(S[2]) 
  ndi <- 1- sign(SpaakSensitivity$S_i[i])*nd 
  ndj <- 1- sign(SpaakSensitivity$S_j[i])*nd 
  
  # construction of a matrix
  NFD.to.add <- as.data.frame(matrix(nrow=(((npoints))), ncol = ncolumn))
  names(NFD.to.add) <-  names.SpaakSensitivity.NFD
  # Copy of the information proper to the line for eah of the points
  for (name.col in names(SpaakSensitivity)){
    NFD.to.add[1:10,name.col] <- SpaakSensitivity[i,name.col]
  }
  # Copy of the NFD approximation for the inferior competitor
  
  NFD.to.add$ND <- ndi
  for ( n in 1:10){
  NFD.to.add$FD[n] <- max(fi[n],fj[n])
  }
 
  # Addition of a column, that will tell us if the Competitive Outcome according to spaak 
  # match the competitive outcome according to Carroll
  NFD.to.add$Coex.Carroll <- "Competition"
  NFD.to.add$Coex.Spaak <- "Competition"
  
  for (n in 1:npoints){
    if(NFD.to.add$S_j[n] < 1 &
       NFD.to.add$S_i[n] < 1){ 
      NFD.to.add$Coex.Carroll[n] <- "Coexistence"}
    
    if(NFD.to.add$FD[n] < NFD.to.add$ND[n]){ 
      NFD.to.add$Coex.Spaak[n] <- "Coexistence"}
    if(1 < NFD.to.add$ND[n]){ 
      NFD.to.add$Coex.Spaak[n] <- "Coexistence"}
  }
  
  # Bind the dataframe NFD.to.add with the general one SpaakSensitivity.NFD
  SpaakSensitivity.NFD <- rbind(SpaakSensitivity.NFD, NFD.to.add)
  
} # end big loop

#Remove the value where spaak and carroll don't agree on the comp.out, to certify the right approximation of the NFD
Spaak.fr.Carroll  <- subset(SpaakSensitivity.NFD,
                          #  SpaakSensitivity.NFD$Coex.Carroll == SpaakSensitivity.NFD$Coex.Spaak,
                            select=c(names(SpaakSensitivity),
                                     "ND","FD"))
Spaak.fr.Carroll$Definition <- "Spaak"
#========================================================================================================
# 2.3 Gather the data under same data.frame
#========================================================================================================


NFD_metaanalysis$ND[NFD_metaanalysis$Initial.Definition== "Zhao" &
                    NFD_metaanalysis$Definition == "Zhao" &
                    NFD_metaanalysis$First.Author == "Zhao"] <- NFD_meta$ND_Zhao[NFD_meta$Initial.Definition =="Zhao"]
NFD_metaanalysis$FD[NFD_metaanalysis$Initial.Definition== "Zhao" &
                      NFD_metaanalysis$Definition == "Zhao" &
                      NFD_metaanalysis$First.Author == "Zhao"] <- NFD_meta$FD_Zhao[NFD_meta$Initial.Definition =="Zhao"]



NFD_metaanalysis <- bind_rows(NFD_metaanalysis,Spaak.fr.Carroll)

# Reorder the column Kingdom 
kd <- c("Phytoplankton","Bacteria","Annual plant","Perennial plant")
NFD_metaanalysis <- arrange(transform(NFD_metaanalysis,
                             Kingdom=factor(Kingdom,levels=kd)),Kingdom)
# Reorder the column Definition
def <- c("Spaak","Carroll","Zhao")
NFD_metaanalysis <- arrange(transform(NFD_metaanalysis,
                             Definition=factor(Definition,levels=def)),Definition)

#######################################################################################################----
## 3.0 Preparation data.frame for graphics
#######################################################################################################----
#========================================================================================================
# 3.1 Setting number of points per subset, the percentage for each interaction and transparensy
#========================================================================================================
#------ Number of point------
# create a data.frame with the number of point per graph
npoint <- data.frame(matrix(ncol=6, nrow=0)) # new data_frame

names(npoint) <- c("Number","Kingdom","Definition",
                   "facilitation","positiv.freq.dep","Neg.freq.dep") 
# Forloop to count the number of NFD for each subset 
for (n in def) {
  for (i in  kd){
    npoint.specific <- data.frame(matrix(ncol=6, nrow=1))
    npoint.fd <- c(NFD_metaanalysis$FD[NFD_metaanalysis$Kingdom == i & NFD_metaanalysis$Definition == n])
    npoint.nd  <- c(NFD_metaanalysis$ND[NFD_metaanalysis$Kingdom == i & NFD_metaanalysis$Definition == n])
    
    facilitation <- c(NFD_metaanalysis$FD[NFD_metaanalysis$Kingdom == i & NFD_metaanalysis$Definition == n &
                                   NFD_metaanalysis$ND > 1])
    Priorityeff <- c(NFD_metaanalysis$FD[NFD_metaanalysis$Kingdom == i & NFD_metaanalysis$Definition == n &
                                  NFD_metaanalysis$ND < 0])
    npoint.fd  <- npoint.fd[!is.na(npoint.fd)]
    npoint.nd  <- npoint.fd[!is.na(npoint.nd)]
    npoint.nfd <- min(length(npoint.fd),length(npoint.nd))
    facilitation <- facilitation[!is.na(facilitation)]
    Priorityeff <- Priorityeff[!is.na(Priorityeff)]
    npoint.specific[1,] <- c(npoint.nfd, i, n,
                             length(facilitation)/npoint.nfd,
                             length(Priorityeff)/npoint.nfd,
                             (npoint.nfd-length(facilitation)-length(Priorityeff))/npoint.nfd)
    npoint <- rbind(npoint,npoint.specific)
  }
}
names(npoint) <- c("Number","Kingdom","Definition",
                   "Facilitation","Positiv.freq.dep","Negativ.freq.dep")

# Adjust the nombre of point based 
npoint$Number <- as.numeric(npoint$Number)
# Spaak and zhao have FD specific to species, not community, so they have twice the number of point 
#npoint$Number[npoint$Definition == "Spaak"] <-
#  ceiling(npoint$Number[ npoint$Definition == "Spaak"]) 
#npoint$Number[npoint$Definition == "Zhao"] <-
# ceiling(npoint$Number[ npoint$Definition == "Zhao"]) 

# For spaak, compute for carroll's data, 
# the NFD are stock under Xs and Ys as they have 10 times the weight of the other NFD
# So as 10 points correspond to 1 community, we need to divided the number of point by 10
npoint$Number[npoint$Kingdom == "Bacteria" & npoint$Definition == "Spaak"] <-
  ceiling(npoint$Number[npoint$Kingdom == "Bacteria" & npoint$Definition == "Spaak"]/10) 
npoint$Number[npoint$Kingdom == "Phytoplankton" & npoint$Definition == "Spaak"] <-
  ceiling(npoint$Number[npoint$Kingdom == "Phytoplankton" & npoint$Definition == "Spaak"]/10) 


# creation od a data_frame with the data_frame npoint to input on the graph 
dat_text <- data.frame(
  label = as.character(npoint$Number),
  Kingdom   = npoint$Kingdom,
  Definition   = npoint$Definition
)
#------ Data frame with the percentage of interactions------
df_interaction <- gather(npoint, c("Facilitation","Positiv.freq.dep","Negativ.freq.dep"),
                         key = "interaction",
                         value="percentage")

#------ Set the transparency------
# to have the same amount of color on eah graph we set the transparency of the dot 
# according to the respectiv number of point on each of them
ntransparency <- data.frame(matrix(ncol=3, nrow=0))
for (n in def) {
  for (i in  kd){
    npoint.specific <- data.frame(matrix(ncol=3, nrow=1))
    npoint.fd <- c(NFD_metaanalysis$FD[NFD_metaanalysis$Kingdom == i & NFD_metaanalysis$Definition == n])
    npoint.nd  <- c(NFD_metaanalysis$ND[NFD_metaanalysis$Kingdom == i & NFD_metaanalysis$Definition == n])
    npoint.fd  <- npoint.fd[!is.na(npoint.fd)]
    npoint.nd  <- npoint.fd[!is.na(npoint.nd)]
    npoint.nfd <- min(length(npoint.fd),length(npoint.nd))
    # if the number is bigger than 1000, the transparency will be too high, so we limit it at 1000
    if (npoint.nfd >1000) {npoint.nfd <- 1000} 
    npoint.specific[1,] <- c(npoint.nfd, i, n)
    ntransparency <- rbind(ntransparency,npoint.specific)
  }
}

names(ntransparency) <- c("Number","Kingdom","Definition")
ntransparency$Number <- as.numeric(ntransparency$Number)
# we define a column in the data_frame NFD_metaanalysis which defines the transparency
for (n in def) {
  for (i in  kd){
    if (ntransparency$Number[ntransparency$Kingdom == i & ntransparency$Definition == n] == 0) next
    
    NFD_metaanalysis$transp[NFD_metaanalysis$Kingdom == i & NFD_metaanalysis$Definition == n] <- 
      (1*min(npoint$Number[1:6]))/ntransparency$Number[ntransparency$Kingdom == i & ntransparency$Definition == n]
  }
}


#========================================================================================================
# 3.2 Setting the outcome of competition for each data points 
#========================================================================================================
# addition of two columns in the data_frame NFD_metaanalysis
NFD_metaanalysis$Condition <- NA
NFD_metaanalysis$Interaction <- NA
NFD_metaanalysis$Color<- NA

# For Spaak 
for ( i in 1:nrow(NFD_metaanalysis)){
  if(NFD_metaanalysis$Definition[i] == "Spaak"){
    NFD_metaanalysis$Condition[i] <- "Competitive Exclusion"
    NFD_metaanalysis$Interaction[i] <- "Negative frequency dependance"
    NFD_metaanalysis$Color[i] <- "red"
    if (is.na(NFD_metaanalysis$ND[i]) | is.na(NFD_metaanalysis$FD[i])) next # no need to change the comp.out as the data won't be represented
    if (NFD_metaanalysis$ND[i] > 0 & 
        NFD_metaanalysis$ND[i] > NFD_metaanalysis$FD[i]){
      NFD_metaanalysis$Condition[i] <- "Coexistence"
      NFD_metaanalysis$Color[i] <- "chartreuse"
    }
    if (NFD_metaanalysis$ND[i] < 0 & 
        abs(NFD_metaanalysis$ND[i]) > NFD_metaanalysis$FD[i]){
      NFD_metaanalysis$Condition[i] <- "Priority Effect"
      NFD_metaanalysis$Color[i] <- "dodgerblue"
    }
    if (NFD_metaanalysis$ND[i] > 1){
      NFD_metaanalysis$Interaction[i] <- "Facilitation"
    }
    if (NFD_metaanalysis$ND[i] < 0){
      NFD_metaanalysis$Interaction[i] <- "Positive frequency dependance"
    }
  }
}

# For Carroll
for ( i in 1:nrow(NFD_metaanalysis)){
  if(NFD_metaanalysis$Definition[i] == "Carroll"){
    NFD_metaanalysis$Condition[i] <- "Competitive Exclusion"
    NFD_metaanalysis$Color[i] <- "red"
    if (is.na(NFD_metaanalysis$S_i[i]) | is.na(NFD_metaanalysis$S_j[i])) next # no need to change the comp.out as the data won't be represented
    if (NFD_metaanalysis$S_i[i] < 1 & 
        NFD_metaanalysis$S_j[i] < 1){
      NFD_metaanalysis$Condition[i] <- "Coexistence"
      NFD_metaanalysis$Color[i] <- "chartreuse"
    }
    if (NFD_metaanalysis$S_i[i] > 1 & 
        NFD_metaanalysis$S_j[i] > 1){
      NFD_metaanalysis$Condition[i] <- "Priority Effect"
      NFD_metaanalysis$Color[i] <- "dodgerblue"
    }
  }
}

# For zhao
for ( i in 1:nrow(NFD_metaanalysis)){
  if(NFD_metaanalysis$Definition[i] == "Zhao"){
    NFD_metaanalysis$Condition[i] <- "Competitive Exclusion"
    NFD_metaanalysis$Color[i] <- "red"
    if (is.na(NFD_metaanalysis$r_i[i]) | is.na(NFD_metaanalysis$r_j[i])) next # no need to change the comp.out as the data won't be represented
    if (NFD_metaanalysis$r_i[i] > 0 & 
        NFD_metaanalysis$r_j[i] > 0){
      NFD_metaanalysis$Condition[i] <- "Coexistence"
      NFD_metaanalysis$Color[i] <- "chartreuse"
    }
    if (NFD_metaanalysis$r_i[i] < 0 & 
        NFD_metaanalysis$r_j[i] < 0){
      NFD_metaanalysis$Condition[i] <- "Priority Effect"
      NFD_metaanalysis$Color[i] <- "dodgerblue"
    }
    if(NFD_metaanalysis$First.Author[i] == "Zhao"){
      NFD_metaanalysis$Condition[NFD_metaanalysis$First.Author == "Zhao"] <- c(NFD_meta$Zhao_Coex[NFD_meta$Initial.Definition =="Zhao"])
      if(NFD_metaanalysis$Condition[i] < 0) {
        NFD_metaanalysis$Condition[i] <- "Competitive Exclusion"
        NFD_metaanalysis$Color[i] <- "red"
      }
      else {
        NFD_metaanalysis$Condition[i] <- "Coexistence"
        NFD_metaanalysis$Color[i] <- "chartreuse"
      }
      }
  }
}




NFD_metaanalysis  <- arrange(transform(NFD_metaanalysis,
                              Condition=factor(Condition,
                                               levels=c("Coexistence",
                                                        "Competitive Exclusion",
                                                        "Priority Effect"
                                               ))),Condition)




