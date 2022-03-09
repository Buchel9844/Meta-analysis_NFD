#######################################################################################################----
# Chessons 2000
#######################################################################################################----
#========================================================================================================
# Adler & al, 2018
#========================================================================================================

# First part of the code is made by Chesson & al, 2018-Wrapper.r
# the only modification was to add the column "Focal.Form2" into the table 'effectD
#After running the CompRegress_wrapper.r file, we used the data.frame effectD produce in the global environment
#ATTNETION, you need to change de stewd to your local file datapluscode in the wrapper

#----Wrapper.r file, to run----  

# this script calls other scripts to compete all stages of analysis
rm(list=ls())

# required packages
req_libs <- c("lme4","tidyr", "dplyr", "ggplot2","ggthemes","viridis","gridExtra",
              "grid","gtable","cowplot","texreg")

# use this function to check if each package is on the local machine
# if a package is installed, it will be loaded
# if any are not, the missing package(s) will be installed and loaded
package_check <- lapply(req_libs, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

# import data
rawD <- read.csv("Data/Adler2018/CompetitionRegressionData071517.csv")

# clean up data
source("Data/Adler2018/CompRegress_cleanData.r")  # this returns a "cleanD" data frame
rm(rawD)  # use cleanD from here on

# scale the competition coefficients
source("Data/Adler2018/CompRegress_scaleAlphas.r")  # this returns "effect" and "response" data frames


#----New script----


Data_Adler <- effectD # copy the effect data.frame from the wrapper file
str(Data_Adler)

Data_Adler$Intra.species <- as.factor(Data_Adler$Intra.species) # transform character type into factor, 
Data_Adler$Target.species <- as.factor(Data_Adler$Target.species) #otherwise the loop don't work
Data_Adler$treatment <- as.character(Data_Adler$treatment) 
Data_Adler$Response2 <- as.character(Data_Adler$Response2) 

Data_Adler$alpha.inter <- (-1*Data_Adler$alpha.inter)
Data_Adler$alpha.intra <- (-1*Data_Adler$alpha.intra)

species_adler <- c(levels(as.factor(Data_Adler$Intra.species))) # create a vector with all the species 


# when treatment is NA, replace it by the response variable measured
# that way e know the difference btw the competition coefficient values
for(i in 1:nrow(Data_Adler)){
  if(is.na(Data_Adler$treatment[i])){
    Data_Adler$treatment[i] <- Data_Adler$Response2[i]
  }
}
Data_Adler <- subset(Data_Adler, select=c(First.Author, Year, treatment,
                                          lab.OR.field,Focal.Form2,
                                          Intra.species, Target.species,
                                          alpha.intra,alpha.inter))
names(Data_Adler) <-  c("First.Author","Year","treatment","lab.OR.field","Focal.Form2",
                        "Intra.species","Target.species",
                        "alpha.intra.i","alpha.inter.i")

# rearrange the data.frame ti have the inter-coeff when target species is now focal (alpha_ij)
Data_Adler_inv <- Data_Adler %>%
  group_by(First.Author, Year, treatment) %>% 
  arrange(treatment,Target.species) %>% 
  select(First.Author, Year, treatment,
         lab.OR.field,Focal.Form2,
         Intra.species, Target.species,
         alpha.intra.i,alpha.inter.i) %>% 
  ungroup(First.Author, Year, treatment)

names(Data_Adler_inv) <-  c("First.Author","Year","treatment","lab.OR.field","Focal.Form2",
                            "Target.species","Intra.species",
                            "alpha.intra.j","alpha.inter.j")

# group the coeff in one data.frame
Data_Adler_Full <- merge(Data_Adler, Data_Adler_inv)

Data_Adler_Full$Intra.species  <- as.character(Data_Adler_Full$Intra.species)
Data_Adler_Full$Target.species <- as.character(Data_Adler_Full$Target.species)
Data_Adler_Full <- subset(Data_Adler_Full, Data_Adler_Full$Intra.species != Data_Adler_Full$Target.species)

# To avoid the repetition of interaction, when specie1=i and specie2=j, 
#you want to avoid the repetition when specie1=j and specie2=i

CopyFull <- Data_Adler_Full 
SaveFull <- data.frame(matrix(ncol=ncol(Data_Adler_Full), nrow=0))

for (i in species_adler){
  if (nrow(CopyFull[CopyFull$Intra.species==i,]) <= 0) next
  
  SaveFull <- rbind(SaveFull, CopyFull[CopyFull$Intra.species==i,])
  CopyFull <- CopyFull[-which(CopyFull$Intra.species==i| CopyFull$Target.species==i),]
}

# when all coefficient are negative, we can change the sign of all of them

Data_Adler <- subset(SaveFull, select=c("First.Author","Year","treatment","lab.OR.field",
                                     "Intra.species","Target.species","Focal.Form2",
                                     "alpha.intra.j","alpha.inter.j",
                                     "alpha.intra.i","alpha.inter.i"))
names(Data_Adler) <- c("First.Author","Year",'treatment',"lab.OR.field",
                      "species.i","species.j","Kingdom",
                       "LV_a_jj","LV_a_ij","LV_a_ii","LV_a_ji")  

Data_Adler$Kingdom <- as.character(Data_Adler$Kingdom)
for (n in 1:nrow(Data_Adler)){
     if (Data_Adler$Kingdom[n] == "annual"){
       Data_Adler$Kingdom[n] <- "Annual plant"}
   }
 
 for (i in c("woody" ,"perennial.herb")){
   for (n in 1:nrow(Data_Adler)){
     if (Data_Adler$Kingdom[n] == i){
       Data_Adler$Kingdom[n] <- "Perennial plant"
     }
       
   }
 }  

Data_Adler$Year <- as.character(Data_Adler$Year)
Data_Adler$Model.use <- "Lotka Volterra"
Data_Adler$Initial.Definition <- "Chesson"


#========================================================================================================
# Verseglou & al, 2018
#========================================================================================================
Vers_df <- read.csv("Data/Coeff_int_Verseglou2018.csv")
Vers_df$focal_species <- as.character(Vers_df$focal_species)
spVers <- Vers_df$focal_species
names(Vers_df) <- c("species.i",spVers)

Data_Vers <- gather(Vers_df,spVers,
                    key="species.j", value="LV_a_ij",na.rm = FALSE,factor_key=TRUE)

Data_Vers <- Data_Vers[order(Data_Vers$species.i),]
Data_Vers_2 <-  Data_Vers[order(Data_Vers$species.j),]
Data_Vers$LV_a_ji <- Data_Vers_2$LV_a_ij
# intra coeff and s,g, lambda

Data_Vers_intra <- subset(Data_Vers,Data_Vers$species.i == Data_Vers$species.j,
                          select=c(species.i, LV_a_ij))

names(Data_Vers_intra) <- c("species.i","LV_a_ii")
Data_Vers <- merge(Data_Vers,Data_Vers_intra)
names(Data_Vers_intra) <- c("species.j","LV_a_jj")
Data_Vers <- merge(Data_Vers,Data_Vers_intra)

Data_Vers <- merge(Data_Vers,Data_Vers_intra)

# Delete when specie1 = i and j 
Data_Vers <- subset(Data_Vers, Data_Vers$species.i != Data_Vers$species.j) 

# To avoid the repetition of interaction, when specie1=i and specie2=j, you want to avoid the repetition when specie1=j and specie2=i
CopyVers<-Data_Vers
SaveVers<-data.frame(matrix(ncol=ncol(Data_Vers), nrow=0))
for (i in spVers){
  SaveVers<-rbind(SaveVers, CopyVers[CopyVers$species.i==i,])
  CopyVers<-CopyVers[-which(CopyVers$species.i==i | CopyVers$species.j==i),]
}

Data_Vers <- SaveVers
#replace the acronyms by real names

spVers_long <- c("Plantago lanceolata", "Rumex acetosa", "Prunella vulgaris", "Fragaria vesca",
                 "Stellaria graminea", "Agrostis capillaris", "Phleum pratense",
                 "Cynosyrus cristatus")
spVers_long <- spVers_long[order(spVers_long)]
Data_Vers$species.i <- as.character(Data_Vers$species.i)
Data_Vers$species.j <- as.character(Data_Vers$species.j)
for (i in 1:length(spVers)){
  Data_Vers$species.i[Data_Vers$species.i == spVers[i]] <- spVers_long[i]
  Data_Vers$species.j[Data_Vers$species.j == spVers[i]] <- spVers_long[i]
}

Data_Vers$Kingdom <- "Perennial plant"
Data_Vers$Model.use <- "Lotka Volterra"
Data_Vers$Initial.Definition <- "Chesson"
Data_Vers$First.Author <- "Verseglou"
Data_Vers$Year <- "2018"
Data_Vers$treatment <- "Control"
Data_Vers$lab.OR.field<- "lab"



#========================================================================================================
# Armitage & al, 2019
#========================================================================================================
Data_Armi <- read.csv("Data/Coeff_int_Armitage2019.csv",sep=",")
Data_Armi$Kingdom <- "Annual plant"
Data_Armi$Model.use <- "Lotka Volterra"
Data_Armi$Initial.Definition <- "Chesson"
Data_Armi$First.Author <- "Armitage"
Data_Armi$Year <- "2019"
Data_Armi$treatment <- "NA"
Data_Armi$lab.OR.field<- "lab"

#========================================================================================================
# Ke & al, 2019
#========================================================================================================

Data_Ke <- read.csv("Data/Coeff_int_Ke2019.csv",sep=",")
Data_Ke$Kingdom <- "Annual plant"
Data_Ke$Model.use <- "Lotka Volterra"
Data_Ke$Initial.Definition <- "Chesson"
Data_Ke$First.Author <- "Ke"
Data_Ke$Year <- "2019"
Data_Ke$lab.OR.field<- "lab"


#######################################################################################################----
# Adler & al, 2007
#######################################################################################################----
#========================================================================================================
# Siefert & al, 2018
#========================================================================================================

Sief_df  <- read.csv("Data/Coeff_int_Siefert2018.csv")

spSief <- levels(as.factor(Sief_df$Focal.species))
Sief_df$Focal.species <- as.character(Sief_df$Focal.species)

Data_Sief<- gather(Sief_df, spSief , key = "species.j", value = "AP_a_ij")
names(Data_Sief) <- c("species.i","treatment","AP_lambda_i",
                    "species.j","AP_a_ij")

Data_Sief <-Data_Sief[order(Data_Sief$species.i),]
Data_Sief_inter <-  Data_Sief[order(Data_Sief$species.j),]
Data_Sief$AP_a_ji <- Data_Sief_inter$AP_a_ij
Data_Sief$AP_lambda_j <- Data_Sief_inter$AP_lambda_i

# intra coeff and s,g, lambda

Data_Sief_intra <- subset(Data_Sief,Data_Sief$species.i == Data_Sief$species.j,
                          select=c("species.i","treatment","AP_a_ij"))
names(Data_Sief_intra) <- c("species.i","treatment","AP_a_ii")
Data_Sief <- merge(Data_Sief,Data_Sief_intra)
names(Data_Sief_intra) <- c("species.j","treatment","AP_a_jj")
Data_Sief <- merge(Data_Sief,Data_Sief_intra)

# Delete when specie1 = i and j 
Data_Sief <- subset(Data_Sief, Data_Sief$species.i != Data_Sief$species.j) 

# To avoid the repetition of interaction, when specie1=i and specie2=j, you want to avoid the repetition when specie1=j and specie2=i
CopySief<-Data_Sief
SaveSief<-data.frame(matrix(ncol=ncol(Data_Sief), nrow=0))
for (i in spSief){
  SaveSief<-rbind(SaveSief, CopySief[CopySief$species.i==i,])
  CopySief<-CopySief[-which(CopySief$species.i==i | CopySief$species.j==i),]
}

Data_Sief <- SaveSief
#replace the acronyms by real names

spSief_long <- c("Trifolium barbigerum","Trifolium bifidum",
                 "Trifolium macraei","Trifolium microdon")

spSief_long <- spSief_long[order(spSief_long)]
Data_Sief$species.i <- as.character(Data_Sief$species.i)
Data_Sief$species.j <- as.character(Data_Sief$species.j)
for (i in 1:length(spSief)){
  Data_Sief$species.i[Data_Sief$species.i == spSief[i]] <- spSief_long[i]
  Data_Sief$species.j[Data_Sief$species.j == spSief[i]] <- spSief_long[i]
}

Data_Sief$Kingdom <- "Annual plant"
Data_Sief$Model.use <- "APM"
Data_Sief$AP_s_i <- 0
Data_Sief$AP_s_j <- 0
Data_Sief$AP_g_i <- 1
Data_Sief$AP_g_j <- 1
Data_Sief$Initial.Definition <- "Adler"
Data_Sief$First.Author <- "Siefert"
Data_Sief$Year <- "2018"
Data_Sief$treatment <- "Control"
Data_Sief$lab.OR.field<- "lab"

#######################################################################################################----
# Carroll 2011
#######################################################################################################----
#========================================================================================================
# Venail & al, 2014
#========================================================================================================

Ven_df <- read.csv("Data/Growth_rates_Venail2014.csv")
spVen <- levels(as.factor(Ven_df$Resident))

#Mean of the replicats

Ven_df<- unite(Ven_df, Resident, Invader,col = "InteractionSp", sep = "-")
# GI is the growth rate of the invader
GI <- tapply(Ven_df$Growth.Invader,
             Ven_df$InteractionSp, mean,na.rm=TRUE)

# GIM is the growth rate of the invader in monoculture
GIM <- tapply(Ven_df$Growth.Invader.in.monoculture,
              Ven_df$InteractionSp, mean,na.rm=TRUE)

# Data frame with the  GR for j ( g_jj and g_ji)
Data_Ven <- data.frame(GI,GIM)
names(Data_Ven)<- c("g_ji","g_jj")
Data_Ven$species <- rownames(Data_Ven)
Data_Ven <- remove_rownames(Data_Ven)


Data_Ven  <- separate(Data_Ven, species,
                      into = c("species.i", "species.j"))

# Addition the GR for i ( g_ii and g_ij)

Data_Ven_2 <- Data_Ven[order(Data_Ven$species.j),]
Data_Ven$g_ii <- Data_Ven_2$g_jj
Data_Ven$g_ij <- Data_Ven_2$g_ji


# Remove communities with the same species twice
Data_Ven <- subset(Data_Ven, Data_Ven[,1] != Data_Ven[,2])

# Remove the similar communities
CopyVen <- Data_Ven
SaveVen <- data.frame(matrix(ncol=ncol(Data_Ven), nrow=0))
for (i in spVen ){
  if (i == "") next
  SaveVen <- rbind(SaveVen, CopyVen[CopyVen$species.i==i,])
  CopyVen <- CopyVen[-which(CopyVen$species.i==i | CopyVen$species.j==i),]
}

Data_Ven <- SaveVen
Data_Ven$Kingdom <- "Phytoplankton"
Data_Ven$Model.use <- "None"
Data_Ven$Initial.Definition <- "Carroll"
Data_Ven$First.Author <- "Venail"
Data_Ven$Year <- "2014"
Data_Ven$treatment <- "Control"
Data_Ven$lab.OR.field<- "lab"

#========================================================================================================
# Tan & al,2016
#========================================================================================================
Tan_df <- read.csv("Data/Growth_rates_Tan2016.csv")
Tan_df$cfd <-  as.numeric(Tan_df$cfd)
Tan_df$nd <-  as.numeric(Tan_df$nd)

spTan <- levels(as.factor(Tan_df$Competitor))
Data_Tan <- data.frame(species.i=c("Pseudomonas fluorescens SBW25"),
                       species.j=levels(as.factor(Tan_df$Competitor)))
Tan_df$Competitor <- as.character(Tan_df$Competitor)
Data_Tan$species.j <- as.character(Data_Tan$species.j)

for (sp in spTan){
  Data_Tan$ND[Data_Tan$species.j ==sp] <- Tan_df$nd[Tan_df$Competitor == sp][1]
  Data_Tan$FD[Data_Tan$species.j ==sp] <- Tan_df$cfd[Tan_df$Competitor == sp][1]
}


spTan_long <- c("Aerogenes hydrophilia","Bacillus cereus","Bacillus pumilus","control",
                "Pseudomonas fluorescens","Pseudomonas putida","Serratia marcescens")

for (i in 1:length(spTan)){
  Data_Tan$species.j[Data_Tan$species.j == spTan[i]] <-spTan_long[i]
}

Data_Tan$Kingdom <- "Bacteria"
Data_Tan$Model.use <- "None"
Data_Tan$Initial.Definition <- "Carroll"
Data_Tan$First.Author <- "Tan"
Data_Tan$Year <- "2016"
Data_Tan$treatment <- "Control"
Data_Tan$lab.OR.field<- "lab"

#========================================================================================================
# Narwani & al, 2013
#========================================================================================================

Narw_df <- read.csv("Data/Sensitivity_Narwani2013.csv")

#implementation of the invasion growth rate
Data_Narwa <- data.frame(Narw_df$Invaded, Narw_df$Invader, Narw_df$Invader.growth.rate,
                           Narw_df$Monoculture.growth.rate..Invader.,
                           Narw_df$Sensitivity.w.negatives)
Data_Narwa <- remove_missing(Data_Narwa)
names(Data_Narwa)<-c("species.i","species.j","g_ji","g_jj","S_i")

Data_Narwa_inv <- Data_Narwa[order(Data_Narwa$species.j),]
Data_Narwa$g_ij <- Data_Narwa_inv$g_ji
Data_Narwa$S_j <- Data_Narwa_inv$S_i
Data_Narwa$g_ii <- Data_Narwa_inv$g_jj

# Remove communities with the same species twice
Data_Narwa <- subset(Data_Narwa, Data_Narwa[,1] != Data_Narwa[,2])

# Remove the similar communities
CopyNar <- data.frame(Data_Narwa)
SaveNar <- data.frame(matrix(ncol=ncol(Data_Narwa), nrow=0))
spNar <- levels(as.factor(Data_Narwa$species.j))
for (i in spNar){
  if (i == "") next
  SaveNar <-rbind(SaveNar, CopyNar[CopyNar$species.i==i,])
  CopyNar <-CopyNar[-which(CopyNar$species.i==i | CopyNar$species.j==i),]
}

Data_Narwa <- SaveNar
Data_Narwa$Kingdom <- "Phytoplankton"
Data_Narwa$Model.use <- "None"
Data_Narwa$Initial.Definition <- "Carroll"
Data_Narwa$First.Author <- "Narwani"
Data_Narwa$Year <- "2013"
Data_Narwa$treatment <- "Control"
Data_Narwa$lab.OR.field<- "lab"

Data_Narwa$species.i <- as.character(Data_Narwa$species.i)
Data_Narwa$species.j <- as.character(Data_Narwa$species.j)

spNar_long <- c("Chlorella sorokiniana","Coelastrum microporum", "Cosmarium turpinii",
                "Elakatothrix viridis", "Scenedesmus acuminatus", "Selenastrum capricornutum",
                "Staurastrum punctulatum","Tetraedron minimum.")
for (i in 1:length(spNar)){
  Data_Narwa$species.i[Data_Narwa$species.i == spNar[i]] <- spNar_long[i]
  Data_Narwa$species.j[Data_Narwa$species.j == spNar[i]] <- spNar_long[i]
}


#========================================================================================================
# Ocampo & al, 2018
#========================================================================================================
Data_Ocam <- data.frame(species.i=c("Azolla pinnata"),species.j=c("Azolla rubra"))
Data_Ocam$g_ii <- 0.037605634	
Data_Ocam$g_jj <- 0.023239437	
Data_Ocam$g_ij <- 0.022112676	
Data_Ocam$g_ji <- -0.004225352
Data_Ocam$Kingdom <- "Annual plant"
Data_Ocam$Kingdom <- "Annual plant"
Data_Ocam$Model.use <- "None"
Data_Ocam$Initial.Definition <- "Carroll"
Data_Ocam$First.Author <- "Ocampo"
Data_Ocam$Year <- "2018"
Data_Ocam$treatment <- "Control"
Data_Ocam$lab.OR.field<- "lab"

#========================================================================================================
# Li & al, 2019
#========================================================================================================
Li_df <- read.csv("Data/NFD_Li2018.csv")

spLi_nat <-levels(as.factor(Li_df$Native))
spLi_inv <-levels(as.factor(Li_df$Invader))


Data_Li <- data.frame(species.i=Li_df$Native,
                      species.j=Li_df$Invader,
                      ND =Li_df$Niche.differences..ND.,
                      FD =Li_df$Relative.fitness.differences..RFD.)

Data_Li$Kingdom <- "Bacteria"
Data_Li$Model.use <- "None"
Data_Li$Initial.Definition <- "Carroll"
Data_Li$First.Author <- "Li"
Data_Li$Year <- "2019"
Data_Li$treatment <- "Control"
Data_Li$lab.OR.field<- "lab"

spLi_long_nat <- c("Novosphingobium sediminicola","Novosphingobium sediminis",
                   "Flectobacillus roseus" ,"Vogesella indigofera " ,
                   "Aeromonas jandaei","Rhodococcus maanshanensis",
                   "Leptothrix discophora" ,"Paenibacillus cineris")
spLi_long_inv <- c("Staphylococcus pasteuri","Bacillus cereus","Serratia marcescens")

Data_Li$species.i <- as.character(Data_Li$species.i)
Data_Li$species.j <- as.character(Data_Li$species.j)

for (i in 1:length(spLi_nat)){
  Data_Li$species.i[Data_Li$species.i == spLi_nat[i]] <- spLi_long_nat[i]
}
for (i in 1:length(spLi_inv)){
  Data_Li$species.j[Data_Li$species.j == spLi_inv[i]] <- spLi_long_inv[i]
}

#========================================================================================================
# Jackrel & al, 2020
#========================================================================================================

Data_Jack <- read.csv("Data/Sensitivity_Jackrel2020.csv",sep=",")
Data_Jack$Kingdom <- "Phytoplankton"
Data_Jack$Model.use <- "None"
Data_Jack$treatment <- "NA"
Data_Jack$Initial.Definition <- "Carroll"
Data_Jack$First.Author <- "Jackrel"
Data_Jack$Year <- "2020"
Data_Jack$lab.OR.field<- "lab"

#========================================================================================================
# Grainger & al, 2020
#========================================================================================================
Data_Grainger<- read.csv("Data/NFD_Grainger2019.csv",sep=",")

Data_Grainger$Kingdom <- "Bacteria"
Data_Grainger$Model.use <- "None"
Data_Grainger$Initial.Definition <- "Carroll"
Data_Grainger$First.Author <- "Grainger"
Data_Grainger$Year <- "2019"
Data_Grainger$lab.OR.field<- "lab"

#========================================================================================================
# Gallego & al, 2019
#========================================================================================================
Data_Gall<- read.csv("Data/NFD_Gallego2019.csv",sep=",")

Data_Gall$Kingdom <- "Phytoplankton"
Data_Gall$species.i <- "Phytoplankton"
Data_Gall$Model.use <- "None"
Data_Gall$treatment <- "NA"
Data_Gall$Initial.Definition <- "Carroll"
Data_Gall$First.Author <- "Gallego"
Data_Gall$Year <- "2019"
Data_Gall$lab.OR.field<- "lab"

#========================================================================================================
# Jia & al, 2020
#========================================================================================================


Data_Jia<- read.csv("Data/Sensitivity_Jia2020.csv",sep=",")

Data_Jia$Kingdom <- "Phytoplankton"
Data_Jia$Model.use <- "None"
Data_Jia$treatment <- "NA"
Data_Jia$Initial.Definition <- "Carroll"
Data_Jia$First.Author <- "Jia"
Data_Jia$Year <- "2020"
Data_Jia$lab.OR.field<- "lab"

#######################################################################################################----
# Godoy & Levine 2014
#######################################################################################################----
#========================================================================================================
# Godoy & Levine, 2014
#========================================================================================================

GodLev_df <- read.csv("Data/Coeff_int_GodLev2014.csv")
GodLev_df$focal_species <- as.character(GodLev_df$focal_species )
spGodLev <- GodLev_df$focal_species 
names(GodLev_df) <- c("focal_species",spGodLev)

GodLev_df<- gather(GodLev_df,spGodLev,
                      key="species.j", value="AP_a_ij",na.rm = FALSE,factor_key=TRUE)

GodLev_df_lambda <- read.csv("Data/Vital_rates_GodLev2014.csv")
Data_GodLev <- merge(GodLev_df,GodLev_df_lambda)
names(Data_GodLev) <- c("species.i","species.j","AP_a_ij",
                        "AP_lambda_i","AP_g_i","AP_s_i")

Data_GodLev <-Data_GodLev[order(Data_GodLev$species.i),]
Data_GodLev_inter <-  Data_GodLev[order(Data_GodLev$species.j),]
Data_GodLev$AP_a_ji <- Data_GodLev_inter$AP_a_ij
Data_GodLev$AP_g_j <- Data_GodLev_inter$AP_g_i
Data_GodLev$AP_s_j <- Data_GodLev_inter$AP_s_i
Data_GodLev$AP_lambda_j <- Data_GodLev_inter$AP_lambda_i

# intra coeff
Data_GodLev_intra <- subset(Data_GodLev,Data_GodLev$species.i == Data_GodLev$species.j,
                            select = c(species.i,AP_a_ij))

Data_GodLev_intra$species.i <- as.character(Data_GodLev_intra$species.i)
names(Data_GodLev_intra) <- c("species.i","AP_a_ii")
Data_GodLev <- merge(Data_GodLev,Data_GodLev_intra)
names(Data_GodLev_intra) <- c("species.j","AP_a_jj")
Data_GodLev <- merge(Data_GodLev,Data_GodLev_intra)

# Delete when specie1 = i and j 
Data_GodLev <- subset(Data_GodLev, Data_GodLev$species.i != Data_GodLev$species.j) 

# To avoid the repetition of interaction, when specie1=i and specie2=j, you want to avoid the repetition when specie1=j and specie2=i
CopyGodLev<-Data_GodLev
SaveGodLev<-data.frame(matrix(ncol=ncol(Data_GodLev), nrow=0))
for (i in spGodLev){
  SaveGodLev<-rbind(SaveGodLev, CopyGodLev[CopyGodLev$species.i==i,])
  CopyGodLev<-CopyGodLev[-which(CopyGodLev$species.i==i | CopyGodLev$species.j==i),]
}

Data_GodLev <- SaveGodLev

#replace the acronyms by real names
Data_GodLev$Kingdom <- "Annual plant"
Data_GodLev$Model.use <- "APM"
Data_GodLev$Initial.Definition <- "Godoy"
Data_GodLev$First.Author <- "Godoy & Levine"
Data_GodLev$Year <- "2014"
Data_GodLev$treatment <- "Control"
Data_GodLev$lab.OR.field<- "field"

#========================================================================================================
# Godoy & al,2014
#========================================================================================================
God_df <- read.csv("Data/Coeff_int_Godoy2014.csv")
God_df$focal_species <- as.character(God_df$focal_species)
spGod <- God_df$focal_species
names(God_df) <- c("species.i",spGod)

Data_God <- gather(God_df,spGod,
                    key="species.j", value="AP_a_ij",na.rm = FALSE,factor_key=TRUE)

Data_God <-Data_God[order(Data_God$species.i),]
Data_God_inter <-  Data_God[order(Data_God$species.j),]
Data_God$AP_a_ji <- Data_God_inter$AP_a_ij

# intra coeff and s,g, lambda
God_df_lambad <- read.csv("Data/Vital_rates_Godoy2014.csv")
names(God_df_lambad) <- c("species.i","AP_lambda_i","AP_s_i","AP_g_i")
Data_God_intra <- subset(Data_God,Data_God$species.i == Data_God$species.j)
names(Data_God_intra) <- c("species.i","species.j","AP_a_ii","AP_a_jj")

Data_God_intra  <- merge(God_df_lambad,Data_God_intra)

Data_God_intra$species.i <- as.character(Data_God_intra$species.i)
Data_God_intra$species.j <- as.character(Data_God_intra$species.j)
Data_God_intra <- data.frame(species.i = rep(c(Data_God_intra$species.i), each=nrow(Data_God_intra)),
                             species.j = rep(c(Data_God_intra$species.j), time=nrow(Data_God_intra)),
                             AP_a_ii = rep(c(Data_God_intra$AP_a_ii), each=nrow(Data_God_intra)),
                             AP_a_jj = rep(c(Data_God_intra$AP_a_jj), time=nrow(Data_God_intra)),
                             AP_lambda_i = rep(c(Data_God_intra$AP_lambda_i), each=nrow(Data_God_intra)),
                             AP_s_i = rep(c(Data_God_intra$AP_s_i), each=nrow(Data_God_intra)),
                             AP_g_i = rep(c(Data_God_intra$AP_g_i ), each=nrow(Data_God_intra)),
                             AP_lambda_j = rep(c(Data_God_intra$AP_lambda_i), time=nrow(Data_God_intra)),
                             AP_s_j = rep(c(Data_God_intra$AP_s_i), time=nrow(Data_God_intra)),
                             AP_g_j = rep(c(Data_God_intra$AP_g_i ), time=nrow(Data_God_intra)))


Data_God <- merge(Data_God,Data_God_intra)

# Delete when specie1 = i and j 
Data_God <- subset(Data_God, Data_God$species.i != Data_God$species.j) 

# To avoid the repetition of interaction, when specie1=i and specie2=j, you want to avoid the repetition when specie1=j and specie2=i
CopyGo<-Data_God
SaveGo<-data.frame(matrix(ncol=ncol(Data_God), nrow=0))
for (i in spGod){
  SaveGo<-rbind(SaveGo, CopyGo[CopyGo$species.i==i,])
  CopyGo<-CopyGo[-which(CopyGo$species.i==i | CopyGo$species.j==i),]
}

Data_God <- SaveGo
#replace the acronyms by real names

spGod_long <- c("Agoseris heterophylla","Agoseris retrorsa","Centaurea melitensis",
                  "Hemizonia congesta","Lasthenia californica","Micropus californicus",
                  "Amsinckia menziesii","Silene gallica","Euphorbia peplus",
                  "Lotus purshianus","Lotus wrangeliensis","Medicago polymorpha",
                  "Erodium botrys","Erodium cicutarim","Geranium carolinianum ",
                  "Salvia columbariae","Anagallis arvensis","Plantago erecta",
                  "Navarretia atractyloides","Navarretia jaredii","Clarkia bottae",
                  "Clarkia purpurea")
spGod_long <- spGod_long[order(spGod_long)]
Data_God$species.i <- as.character(Data_God$species.i)
Data_God$species.j <- as.character(Data_God$species.j)
for (i in 1:length(spGod)){
  Data_God$species.i[Data_God$species.i == spGod[i]] <- spGod_long[i]
  Data_God$species.j[Data_God$species.j == spGod[i]] <- spGod_long[i]
}

Data_God$Kingdom <- "Annual plant"
Data_God$Model.use <- "APM"
Data_God$Initial.Definition <- "Godoy"
Data_God$First.Author <- "Godoy"
Data_God$Year <- "2014"
Data_God$treatment <- "Control"
Data_God$lab.OR.field<- "field"

#========================================================================================================
# Germain & al, 2016
#========================================================================================================
Ger_df <- read.csv("Data/Coeff_int_Germain2016.csv")
Ger_df_g_s <- read.csv("Data/Vital_rates_Germain2016.csv")

Ger_df <- merge(Ger_df,Ger_df_g_s)

Ger_df$species.i <- as.character(Ger_df$species.i)
Ger_df$species.j <- as.character(Ger_df$species.j)
Ger_df$species.k <- as.character(Ger_df$species.k)
Ger_df$traitment <- as.character(Ger_df$traitment)

Data_Germ <- data.frame(species.i = c(Ger_df$species.i, Ger_df$species.i),
                     species.j = c(Ger_df$species.j, Ger_df$species.k),
                     treatment = c(Ger_df$traitment, Ger_df$traitment),
                     AP_lambda_i = c(Ger_df$Lamb_i, Ger_df$Lamb_i),
                     AP_lambda_j = c(Ger_df$Lamb_j, Ger_df$Lamb_k),
                     AP_a_ii = c(Ger_df$alpha_ii, Ger_df$alpha_ii),
                     AP_a_jj = c(Ger_df$alpha_jj, Ger_df$alpha_kk),
                     AP_a_ij = c(Ger_df$alpha_ij, Ger_df$alpha_ik),
                     AP_a_ji = c(Ger_df$alpha_ji, Ger_df$alpha_ki),
                     AP_g_i = c(Ger_df$gi, Ger_df$gi),
                     AP_g_j = c(Ger_df$gj, Ger_df$gk),
                     AP_s_i = c(Ger_df$si, Ger_df$si),
                     AP_s_j = c(Ger_df$sj, Ger_df$sk)
                     )

Data_Germ$Kingdom <- "Annual plant"
Data_Germ$Model.use <- "APM"
Data_Germ$Initial.Definition <- "Godoy"
Data_Germ$First.Author <- "Germain"
Data_Germ$Year <- "2016"
Data_Germ$lab.OR.field<- "field"

#========================================================================================================
# Rey & al, 2017
#========================================================================================================
Rey_df <- read.csv("Data/Coeff_int_Rey2017.csv")
str(Rey_df)
Data_Rey <- bind_rows(subset(Rey_df,Rey_df$Species=="Brachypodium distachyon",select=c(Humid,Dry)),
                      subset(Rey_df,Rey_df$Species=="Brachypodium hybridum",select=c(Humid,Dry)))
Data_Rey <- data.frame(t(Data_Rey))
Data_Rey <- remove_rownames(Data_Rey)

names(Data_Rey) <- c("AP_g_i","AP_lambda_i","AP_a_ii","AP_a_ij",
                     "AP_g_j","AP_lambda_j","AP_a_jj","AP_a_ji")

Data_Rey$species.i <- "Brachypodium distachyon"
Data_Rey$species.j <- "Brachypodium hybridum"
Data_Rey$AP_s_i <- 0
Data_Rey$AP_s_j <- 0
Data_Rey$Kingdom <- "Annual plant"
Data_Rey$Model.use <- "APM"
Data_Rey$Initial.Definition <- "Godoy"
Data_Rey$First.Author <- "Rey"
Data_Rey$Year <- "2017"
Data_Rey$treatment <- c("Humid","Dry")
Data_Rey$lab.OR.field<- "field"

#========================================================================================================
# Cardinaux & al, 2018
#========================================================================================================
Car_df <- read.csv("Data/Coeff_int_Cardinaux2018.csv")
spCar <- levels(as.factor(Car_df$focal_species))
Car_df$species.j <-  "Plantago alpina"
Car_df$species.k <-  "Poa alpina"
for (i in 1:nrow(Car_df)){
if(Car_df$focal_species_origin[i] =="High"){
  Car_df$species.j[i] <-  "Plantago lanceolata"
  Car_df$species.k[i] <-  "Poa trivialis"
}
}
Car_df$focal_species <- as.character(Car_df$focal_species)
Car_df$soil_viota_origin <- as.character(Car_df$soil_viota_origin)
Data_Car <- data_frame(species.i = c(Car_df$focal_species, Car_df$focal_species),
                       species.j = c(Car_df$species.j, Car_df$species.k),
                       treatment = c(Car_df$soil_viota_origin,Car_df$soil_viota_origin),
                       AP_lambda_i = c(Car_df$lamb_i, Car_df$lamb_i),
                       AP_a_ii = c(Car_df$a_ii, Car_df$a_ii),
                       AP_a_ij = c(Car_df$a_ij, Car_df$a_ik))
Data_Car <- Data_Car[order(Data_Car$species.i),]
Data_Car_inter <- Data_Car[order(Data_Car$species.j),]
Data_Car$AP_a_jj <- Data_Car_inter$AP_a_ii
Data_Car$AP_a_ji <- Data_Car_inter$AP_a_ij
Data_Car$AP_lambda_j <- Data_Car_inter$AP_lambda_i

# Delete when specie1 = i and j 
Data_Car <- subset(Data_Car, Data_Car$species.i != Data_Car$species.j) 

# To avoid the repetition of interaction, when specie1=i and specie2=j, you want to avoid the repetition when specie1=j and specie2=i
CopyCar<-Data_Car
SaveCar<-data.frame(matrix(ncol=ncol(Data_Car), nrow=0))
for (i in spCar){
  SaveCar<-rbind(SaveCar, CopyCar[CopyCar$species.i==i,])
  CopyCar<-CopyCar[-which(CopyCar$species.i==i | CopyCar$species.j==i),]
}

Data_Car <- SaveCar
Data_Car$AP_g_i <- 1
Data_Car$AP_s_i <- 0
Data_Car$AP_g_j <- 1
Data_Car$AP_s_j <- 0
Data_Car$Kingdom <- "Perennial plant"
Data_Car$Model.use <- "APM"
Data_Car$Initial.Definition <- "Godoy"
Data_Car$First.Author <- "Cardinaux"
Data_Car$Year <- "2018"
Data_Car$treatment <- "Control"
Data_Car$lab.OR.field<- "field"
#========================================================================================================
# Hallett & al, 2018
#========================================================================================================

Hal_df <- read.csv("Data/Coeff_int_Hallet2018.csv")

Hal_df$focal_species <- as.character(Hal_df$focal_species )
spHal <- Hal_df$focal_species 
names(Hal_df) <- c("species.i","AP_lambda_i","AP_beta_i",spHal)
Hal_df <- subset(Hal_df, select=c("species.i","AP_lambda_i",spHal))
Data_Hal <- gather(Hal_df,spHal,
                   key="species.j", value="AP_a_ij",na.rm = FALSE,factor_key=TRUE)

Data_Hal <-Data_Hal[order(Data_Hal$species.i),]
Data_Hal_inter <-  Data_Hal[order(Data_Hal$species.j),]
Data_Hal$AP_a_ji <- Data_Hal_inter$AP_a_ij

# intra coeff and s,g, lambda

Data_Hal_intra <- subset(Data_Hal,Data_Hal$species.i == Data_Hal$species.j)
names(Data_Hal_intra) <- c("species.i","AP_lambda_i",
                           "species.j","AP_a_ii","AP_a_jj")

Data_Hal_intra$species.i <- as.character(Data_Hal_intra$species.i)
Data_Hal_intra$species.j <- as.character(Data_Hal_intra$species.j)
Data_Hal_intra <- data.frame(species.i = rep(c(Data_Hal_intra$species.i), each=nrow(Data_Hal_intra)),
                             species.j = rep(c(Data_Hal_intra$species.j), time=nrow(Data_Hal_intra)),
                             AP_a_ii = rep(c(Data_Hal_intra$AP_a_ii), each=nrow(Data_Hal_intra)),
                             AP_a_jj = rep(c(Data_Hal_intra$AP_a_jj), time=nrow(Data_Hal_intra)),
                             AP_lambda_i = rep(c(Data_Hal_intra$AP_lambda_i), each=nrow(Data_Hal_intra)),
                             AP_lambda_j = rep(c(Data_Hal_intra$AP_lambda_i), times=nrow(Data_Hal_intra)))



Data_Hal <- merge(Data_Hal,Data_Hal_intra)

# Delete when specie1 = i and j 
Data_Hal <- subset(Data_Hal, Data_Hal$species.i != Data_Hal$species.j) 

# To avoid the repetition of interaction, when specie1=i and specie2=j, you want to avoid the repetition when specie1=j and specie2=i
CopyHal<-Data_Hal
SaveHal<-data.frame(matrix(ncol=ncol(Data_Hal), nrow=0))
for (i in spHal){
  SaveHal<-rbind(SaveHal, CopyHal[CopyHal$species.i==i,])
  CopyHal<-CopyHal[-which(CopyHal$species.i==i | CopyHal$species.j==i),]
}

Data_Hal <- SaveHal
#replace the acronyms by real names

spHal_long <- c("Bromus hordeaceus","Calycadenia multiglandulosa",
                "Lasthenia californica","Microseris douglasii",
                "Plantago erecta","Vulpia microstachys")
spHal_long <- spHal_long[order(spHal_long)]
Data_Hal$species.i <- as.character(Data_Hal$species.i)
Data_Hal$species.j <- as.character(Data_Hal$species.j)
for (i in 1:length(spHal)){
  Data_Hal$species.i[Data_Hal$species.i == spHal[i]] <- spHal_long[i]
  Data_Hal$species.j[Data_Hal$species.j == spHal[i]] <- spHal_long[i]
}

Data_Hal$AP_g_i <- 1
Data_Hal$AP_s_i <- 0
Data_Hal$AP_g_j <- 1
Data_Hal$AP_s_j <- 0
Data_Hal$Kingdom <- "Annual plant"
Data_Hal$Model.use <- "APM"
Data_Hal$Initial.Definition <- "Godoy"
Data_Hal$First.Author <- "Hallet"
Data_Hal$Year <- "2018"
Data_Hal$treatment <- "Control"
Data_Hal$lab.OR.field<- "field"

#========================================================================================================
# Lanuza & al, 2018
#========================================================================================================

Lanuz_df <- read.csv("Data/Coeff_int_Lanuza2018.csv")

Lanuz_df$focal_species <- as.character(Lanuz_df$focal_species )
spLanuz <- Lanuz_df$focal_species 
names(Lanuz_df) <- c("species.i","AP_lambda_i",spLanuz)

Data_Lanuz <- gather(Lanuz_df,spLanuz,
                            key="species.j", value="AP_a_ij",na.rm = FALSE,factor_key=TRUE)

Data_Lanuz <-Data_Lanuz[order(Data_Lanuz$species.i),]
Data_Lanuz_inter <-  Data_Lanuz[order(Data_Lanuz$species.j),]
Data_Lanuz$AP_a_ji <- Data_Lanuz_inter$AP_a_ij
Data_Lanuz$AP_lambda_j <- Data_Lanuz_inter$AP_lambda_i

# intra coeff and s,g, lambda

Data_Lanuz_intra <- subset(Data_Lanuz,Data_Lanuz$species.i == Data_Lanuz$species.j)
names(Data_Lanuz_intra) <- c("species.i","AP_lambda_i","species.j","AP_a_ii","AP_a_jj")
Data_Lanuz_intra$species.i <- as.character(Data_Lanuz_intra$species.i)
Data_Lanuz_intra$species.j <- as.character(Data_Lanuz_intra$species.j)
Data_Lanuz_intra <- data.frame(species.i = rep(c(Data_Lanuz_intra$species.i), each=nrow(Data_Lanuz_intra)),
                                      species.j = rep(c(Data_Lanuz_intra$species.j), time=nrow(Data_Lanuz_intra)),
                                      AP_a_ii = rep(c(Data_Lanuz_intra$AP_a_ii), each=nrow(Data_Lanuz_intra)),
                                      AP_a_jj = rep(c(Data_Lanuz_intra$AP_a_jj), time=nrow(Data_Lanuz_intra)),
                                      AP_lambda_i = rep(c(Data_Lanuz_intra$AP_lambda_i), each=nrow(Data_Lanuz_intra)))
                                      #AP_s_i = rep(c(Data_Lanuz_intra$AP_s_i), each=nrow(Data_Lanuz_intra)),
                                      #AP_g_i = rep(c(Data_Lanuz_intra$AP_g_i ), each=nrow(Data_Lanuz_intra))
                               


Data_Lanuz <- merge(Data_Lanuz,Data_Lanuz_intra)

# Delete when specie1 = i and j 
Data_Lanuz <- subset(Data_Lanuz, Data_Lanuz$species.i != Data_Lanuz$species.j) 

# To avoid the repetition of interaction, when specie1=i and specie2=j, you want to avoid the repetition when specie1=j and specie2=i
CopyLanuz<-Data_Lanuz
SaveLanuz<-data.frame(matrix(ncol=ncol(Data_Lanuz), nrow=0))
for (i in spLanuz){
  SaveLanuz<-rbind(SaveLanuz, CopyLanuz[CopyLanuz$species.i==i,])
  CopyLanuz<-CopyLanuz[-which(CopyLanuz$species.i==i | CopyLanuz$species.j==i),]
}

Data_Lanuz <- SaveLanuz
#replace the acronyms by real names

Data_Lanuz$AP_g_i <- 1
Data_Lanuz$AP_s_i <- 0
Data_Lanuz$AP_g_j <- 1
Data_Lanuz$AP_s_j <- 0
Data_Lanuz$Kingdom <- "Annual plant"
Data_Lanuz$Model.use <- "APM"
Data_Lanuz$Initial.Definition <- "Godoy"
Data_Lanuz$First.Author <- "Lanuza"
Data_Lanuz$Year <- "2018"
Data_Lanuz$treatment <- "Control"
Data_Lanuz$lab.OR.field<- "field"

#========================================================================================================
# Matias & al, 2018
#========================================================================================================
Matias_df <- read.csv("Data/Coeff_int_Matias2018.csv")
Matias_df$species.i <- as.character(Matias_df$species.i)
spMatias <- Matias_df$species.i 
names(Matias_df) <- c("species.i",spMatias)

Matias_df<- gather(Matias_df,spMatias,
                   key="species.j", value="AP_a_ij",na.rm = FALSE,factor_key=TRUE)

Matias_df_lambda <- read.csv("Data/Vital_rates_Matias2018.csv")
Data_Matias <- left_join(Matias_df,Matias_df_lambda,by=c("species.i"))

names(Data_Matias) <- c("species.i","species.j","AP_a_ij",
                        "AP_lambda_i","AP_g_i","AP_s_i")

Data_Matias <-Data_Matias[order(Data_Matias$species.i),]
Data_Matias_inter <-  Data_Matias[order(Data_Matias$species.j),]
Data_Matias$AP_a_ji <- Data_Matias_inter$AP_a_ij
Data_Matias$AP_g_j <- Data_Matias_inter$AP_g_i
Data_Matias$AP_s_j <- Data_Matias_inter$AP_s_i
Data_Matias$AP_lambda_j <- Data_Matias_inter$AP_lambda_i

# intra coeff
Data_Matias_intra <- subset(Data_Matias,Data_Matias$species.i == Data_Matias$species.j,
                            select = c(species.i,AP_a_ij))

Data_Matias_intra$species.i <- as.character(Data_Matias_intra$species.i)
names(Data_Matias_intra) <- c("species.i","AP_a_ii")
Data_Matias <- merge(Data_Matias,Data_Matias_intra)
names(Data_Matias_intra) <- c("species.j","AP_a_jj")
Data_Matias <- merge(Data_Matias,Data_Matias_intra)

# Delete when specie1 = i and j 
Data_Matias <- subset(Data_Matias, Data_Matias$species.i != Data_Matias$species.j) 

# To avoid the repetition of interaction, when specie1=i and specie2=j, you want to avoid the repetition when specie1=j and specie2=i
CopyMatias<-Data_Matias
SaveMatias<-data.frame(matrix(ncol=ncol(Data_Matias), nrow=0))
for (i in spMatias){
  SaveMatias<-rbind(SaveMatias, CopyMatias[CopyMatias$species.i==i,])
  CopyMatias<-CopyMatias[-which(CopyMatias$species.i==i | CopyMatias$species.j==i),]
}

Data_Matias <- SaveMatias
#replace the acronyms by real names

Data_Matias$Kingdom <- "Annual plant"
Data_Matias$Model.use <- "APM"
Data_Matias$Initial.Definition <- "Godoy"
Data_Matias$First.Author <- "Matias"
Data_Matias$Year <- "2018"
Data_Matias$treatment <- "Control"
Data_Matias$lab.OR.field<- "field"

#========================================================================================================
# Petry & al, 2018
#========================================================================================================

Petry_df <- read.csv("Data/Coeff_int_Petry2018.csv")

spPetry <- c(as.character(Petry_df$species.i))

Petry_df_demo <-read.csv("Data/Vital_rates_Petry2018.csv")
Petry_df_demo$AP_lambda_i <- Petry_df_demo$lambda * (1 - Petry_df_demo$prop)
Petry_df_demo <- subset(Petry_df_demo,select=c("species" ,"s","g","eta.ant" ))
names(Petry_df_demo) <- c("species.i","AP_s_i","AP_g_i","AP_lambda_i")
Petry_df <- merge(Petry_df,Petry_df_demo)

Data_Petry <- gather(Petry_df,spPetry,
                     key="species.j", value="AP_a_ij",na.rm = FALSE,factor_key=TRUE)

# intra coeff and s,g, lambda

Data_Petry_intra <- subset(Data_Petry,Data_Petry$species.i == Data_Petry$species.j)

Data_Petry_intra$species.i <- as.character(Data_Petry_intra$species.i)
Data_Petry_intra$species.j <- as.character(Data_Petry_intra$species.j)
Data_Petry_intra <- data.frame(species.i = rep(c(Data_Petry_intra$species.i), each=nrow(Data_Petry_intra)),
                                      species.j = rep(c(Data_Petry_intra$species.j), time=nrow(Data_Petry_intra)),
                                      AP_a_ii = rep(c(Data_Petry_intra$AP_a_ij), each=nrow(Data_Petry_intra)))
                                      

Data_Petry <- merge(Data_Petry,Data_Petry_intra)
                                      
Data_Petry <-Data_Petry[order(Data_Petry$species.i),]
Data_Petry_inter <-  Data_Petry[order(Data_Petry$species.j),]
Data_Petry$AP_a_ji <- Data_Petry_inter$AP_a_ij
Data_Petry$AP_g_j <- Data_Petry_inter$AP_g_i
Data_Petry$AP_s_j <- Data_Petry_inter$AP_s_i
Data_Petry$AP_lambda_j <- Data_Petry_inter$AP_lambda_i
Data_Petry$AP_a_jj <- Data_Petry_inter$AP_a_ii

# Delete when specie1 = i and j 
Data_Petry <- subset(Data_Petry, Data_Petry$species.i != Data_Petry$species.j) 

# To avoid the repetition of interaction, when specie1=i and specie2=j, you want to avoid the repetition when specie1=j and specie2=i
CopyPetry<-Data_Petry
SavePetry<-data.frame(matrix(ncol=ncol(Data_Petry), nrow=0))
for (i in spPetry){
  SavePetry<-rbind(SavePetry, CopyPetry[CopyPetry$species.i==i,])
  CopyPetry<-CopyPetry[-which(CopyPetry$species.i==i | CopyPetry$species.j==i),]
}

Data_Petry <- SavePetry
#replace the acronyms by real names

spPetry_long <- c("Acmispon wrangelianus","Amsinckia menziesii","Centaurea melitensis",
                  "Clarkia purpurea","Erodium botrys",
                  "Erodium cicutarium","Euphorbia spathulata",
                  "Lasthenia californica","Plantago erecta",
                  "Salvia columbariae","Uropappus lindleyi")
spPetry_long <- spPetry_long[order(spPetry_long)]
Data_Petry$species.i <- as.character(Data_Petry$species.i)
Data_Petry$species.j <- as.character(Data_Petry$species.j)
for (i in 1:length(spPetry)){
  Data_Petry$species.i[Data_Petry$species.i == spPetry[i]] <- spPetry_long[i]
  Data_Petry$species.j[Data_Petry$species.j == spPetry[i]] <- spPetry_long[i]
}

Data_Petry$Kingdom <- "Annual plant"
Data_Petry$Model.use <- "APM"
Data_Petry$Initial.Definition <- "Godoy"
Data_Petry$First.Author <- "Petry"
Data_Petry$Year <- "2018"
Data_Petry$treatment <- "Control"
Data_Petry$lab.OR.field <- "field"
Data_Petry$AP_g.lambda_i <-  Data_Petry$AP_g_i*Data_Petry$AP_lambda_i
Data_Petry$AP_g.lambda_j <-  Data_Petry$AP_g_j*Data_Petry$AP_lambda_j


#========================================================================================================
# Wainwright & al, 2018
#========================================================================================================
Wain_df  <- read.csv("Data/Coeff_int_Wainwirght2018.csv")
Wain_df <- subset(Wain_df, select=c("Reserve.treatment.combination","Focal","Parameter","MLE.estimate"))
spWain <- levels(as.factor(Wain_df$Focal))

Wain_df <- spread(Wain_df, Parameter, MLE.estimate)# spread the data.frame

names(Wain_df) <- c("treatment","species.i", "A", "H","AP_a_ii","Other", "T", "W") # dataframe for a_ij
Wain_df <- subset(Wain_df,select=c("treatment","AP_a_ii","species.i", "A", "H", "T", "W"))

Wain_df <- gather(Wain_df, 'A','H','T','W', key = "species.j", value = "AP_a_ij")

#add the lambda
Wain_lambda <- read.csv("Data/Lambdas_Wainwright2018.csv")
names(Wain_lambda) <- c("treatment",spWain)
Wain_lambda <- gather(Wain_lambda,spWain,key = "species.i", value = "AP_lambda_i")
Wain_df <- merge(Wain_df,Wain_lambda)
# add the si and gi
Wain_df_g_s  <- read.csv("Data/Vital_rates_Wainwright2018.csv")
Wain_df_g_s  <- gather(Wain_df_g_s , c("B","P"),key= "treatment",value = "value") # spread the data.frame
Wain_df_g_s  <- spread(Wain_df_g_s , parameter,value) # spread the data.frame
Wain_df_g_s <- rbind(Wain_df_g_s,Wain_df_g_s)
Wain_df_g_s  <- Wain_df_g_s[order(Wain_df_g_s$species.i),]
Wain_df_g_s$treatment <- c(rep(c("BC","PC","BW","PW"),times=(nrow(Wain_df_g_s)/4)))
names(Wain_df_g_s) <- c("species.i","treatment","AP_g_i","AP_s_i")
Wain_df <- merge(Wain_df,Wain_df_g_s)
Wain_df_g_s_2 <-Wain_df_g_s
names(Wain_df_g_s_2) <- c("species.j","treatment","AP_g_j","AP_s_j")
Wain_df <- merge(Wain_df,Wain_df_g_s_2)


Data_Wain  <-Wain_df[order(Wain_df$species.i),]
Data_Wain_inter <-  Data_Wain[order(Data_Wain$species.j),]
Data_Wain$AP_a_ji <- Data_Wain_inter$AP_a_ij
Data_Wain$AP_a_jj <- Data_Wain_inter$AP_a_ii
Data_Wain$AP_a_jj <- Data_Wain_inter$AP_a_ii
Data_Wain$AP_lambda_j <- Data_Wain_inter$AP_lambda_i
Data_Wain$AP_g_j <- Data_Wain_inter$AP_g_i
Data_Wain$AP_s_j <- Data_Wain_inter$AP_s_i

# Delete when specie1 = i and j 
Data_Wain <- subset(Data_Wain, Data_Wain$species.i != Data_Wain$species.j) 

# To avoid the repetition of interaction, when specie1=i and specie2=j, you want to avoid the repetition when specie1=j and specie2=i
CopyWain<-Data_Wain
SaveWain<-data.frame(matrix(ncol=ncol(Data_Wain), nrow=0))
for (i in spWain){
  SaveWain<-rbind(SaveWain, CopyWain[CopyWain$species.i==i,])
  CopyWain<-CopyWain[-which(CopyWain$species.i==i | CopyWain$species.j==i),]
}

Data_Wain <- SaveWain
#replace the acronyms by real names

spWain_long <- c("Arctotheca calendula","Hypochaeris glabra",
                 "Trachymene cyanopetala","Waitzia acuminata")

spWain_long <- spWain_long[order(spWain_long)]
Data_Wain$species.i <- as.character(Data_Wain$species.i)
Data_Wain$species.j <- as.character(Data_Wain$species.j)
for (i in 1:length(spWain)){
  Data_Wain$species.i[Data_Wain$species.i == spWain[i]] <- spWain_long[i]
  Data_Wain$species.j[Data_Wain$species.j == spWain[i]] <- spWain_long[i]
}

Data_Wain$Kingdom <- "Annual plant"
Data_Wain$Model.use <- "APM"
Data_Wain$Initial.Definition <- "Godoy"
Data_Wain$First.Author <- "Wainwright"
Data_Wain$Year <- "2018"
Data_Wain$lab.OR.field<- "field"

#========================================================================================================
# Blackford & al, 2020
#========================================================================================================


Data_Black <- read.csv("Data/Coeff_int_Blackford2020.csv",sep=",")

Data_Black <- as.data.frame(t(as.data.frame(apply(Data_Black,2,median))))
Data_Black <-remove_rownames(Data_Black)
Data_Black[,3:8] <- as.numeric(Data_Black[,3:8])

Data_Black$Kingdom <- "Annual plant"
Data_Black$Model.use <- "APM"
Data_Black$treatment <- "NA"
Data_Black$Initial.Definition <- "Godoy"
Data_Black$First.Author <- "Blackford"
Data_Black$Year <- "2020"

Data_Black$lab.OR.field<- "lab"

#######################################################################################################----
# Zhao & al, 2016
#######################################################################################################----
Zhao_df <- read.csv("Data/NFD_Zhao2016.csv")

Data_Zhao <- data.frame(species.i=Zhao_df$strain.pair.id,
                        species.j=Zhao_df$microcosm.pair.id,
                        treatment= Zhao_df$sympatry.allopatry,
                        ND_Zhao = Zhao_df$niche.overlap,
                        FD_Zhao = Zhao_df$fitness.diff,
                        Zhao_Coex =  Zhao_df$strength.coexistence)

Data_Zhao$Kingdom <- "Bacteria"
Data_Zhao$Model.use <- "None"
Data_Zhao$Initial.Definition <- "Zhao"
Data_Zhao$First.Author <- "Zhao"
Data_Zhao$Year <- "2016"
Data_Zhao$lab.OR.field<- "lab"

#######################################################################################################----
# Bimler & al, 2018
#######################################################################################################----
Bim_df_lamb <- read.csv("Data/Lambdas_Bimler2018.csv")
Bim_df_g_s <- read.csv("Data/Vital_rates_Bimler2018.csv")
Bim_df <- read.csv("Data/Coeff_int_Bimler2018.csv")
names(Bim_df) <- c("species.i","species.j","AP_a_ii_exp",
                   "AP_a_ij_exp","AP_a_jj_exp","AP_a_ji_exp","treatment")
spBim <- levels(as.factor(Bim_df$species.i))

#add the lamb
Bim_df_lambi <- Bim_df_lamb
names(Bim_df_lambi) <- c("species.i","AP_lambda_i_exp")
Bim_df <- merge(Bim_df,Bim_df_lambi)

Bim_df_lambj <- Bim_df_lamb
names(Bim_df_lambj) <- c("species.j","AP_lambda_j_exp")
Bim_df <- merge(Bim_df,Bim_df_lambj)

#add the survival rate
Bim_df_g_s$resv <- as.character(Bim_df_g_s$resv)
Bim_df_g_s$resv[Bim_df_g_s$resv== "K"] <- "Bendering"
Bim_df_g_s$resv[Bim_df_g_s$resv== "P"] <- "Perenjori"

Bim_df_g_si <- Bim_df_g_s
names(Bim_df_g_si) <- c("species.i","treatment","AP_g_i_exp","AP_s_i_exp")
Bim_df <- merge(Bim_df,Bim_df_g_si)

Bim_df_g_sj <- Bim_df_g_s
names(Bim_df_g_sj) <- c("species.j","treatment","AP_g_j_exp","AP_s_j_exp")
Bim_df <- merge(Bim_df,Bim_df_g_sj)

Data_Bim  <- Bim_df
# Delete when specie1 = i and j 
Data_Bim <- subset(Data_Bim, Data_Bim$species.i != Data_Bim$species.j) 

# To avoid the repetition of interaction, when specie1=i and specie2=j, you want to avoid the repetition when specie1=j and specie2=i
CopyBim<-Data_Bim
SaveBim<-data.frame(matrix(ncol=ncol(Data_Bim), nrow=0))
for (i in spBim){
  SaveBim<-rbind(SaveBim, CopyBim[CopyBim$species.i==i,])
  CopyBim<-CopyBim[-which(CopyBim$species.i==i | CopyBim$species.j==i),]
}

Data_Bim <- SaveBim
#replace the acronyms by real names

spBim_long <- c("Arctotheca calendula","Hypochaeris glabra",
                "Trachymene cyanopetala","Waitzia acuminata")

spBim_long <- spBim_long[order(spBim_long)]
Data_Bim$species.i <- as.character(Data_Bim$species.i)
Data_Bim$species.j <- as.character(Data_Bim$species.j)
for (i in 1:length(spBim)){
  Data_Bim$species.i[Data_Bim$species.i == spBim[i]] <- spBim_long[i]
  Data_Bim$species.j[Data_Bim$species.j == spBim[i]] <- spBim_long[i]
}

Data_Bim$Kingdom <- "Annual plant"
Data_Bim$Model.use <- "AP2"
Data_Bim$Initial.Definition <- "Bimler"
Data_Bim$First.Author <- "Bimler"
Data_Bim$Year <- "2018"
Data_Bim$treatment <- "Control"
Data_Bim$lab.OR.field <- "field"

#######################################################################################################----
# Zhang & al, 2019
#######################################################################################################----

Zhang_df <- read.csv("Data/Coeff_int_Zhang2019.csv")

Zhang_df.lambda<- subset(Zhang_df,par=="intercept", select=c("spp_target","spp_comp","Estimate"))
Zhang_df.ii<-subset(Zhang_df,par=="alphaii", select=c("spp_target","spp_comp","Estimate"))
Zhang_df.ij<-subset(Zhang_df,par=="alphaij", select=c("spp_target","spp_comp","Estimate"))
Zhang_df <- right_join(Zhang_df.lambda,Zhang_df.ii, by=c("spp_target","spp_comp"))
Zhang_df <-right_join(Zhang_df,Zhang_df.ij, by=c("spp_target","spp_comp"))

head(Zhang_df)
names(Zhang_df) <- c("species.i","species.j","AP_lambda_i","AP_a_ii","AP_a_ij")

Data_Zhang  <-Zhang_df[order(Zhang_df$species.i),]
Data_Zhang_inter <-  Data_Zhang[order(Data_Zhang$species.j),]
Data_Zhang$AP_a_ji <- Data_Zhang_inter$AP_a_ij
Data_Zhang$AP_a_jj <- Data_Zhang_inter$AP_a_ii
Data_Zhang$AP_lambda_j <- Data_Zhang_inter$AP_lambda_i


# Delete when specie1 = i and j 
Data_Zhang <- subset(Data_Zhang, Data_Zhang$species.i != Data_Zhang$species.j) 

# To avoid the repetition of interaction, when specie1=i and specie2=j, you want to avoid the repetition when specie1=j and specie2=i
spZhang <- levels(as.factor(Data_Zhang$species.i))
CopyZhang<-Data_Zhang
SaveZhang<-data.frame(matrix(ncol=ncol(Data_Zhang), nrow=0))
for (i in spZhang){
  SaveZhang<-rbind(SaveZhang, CopyZhang[CopyZhang$species.i==i,])
  CopyZhang<-CopyZhang[-which(CopyZhang$species.i==i | CopyZhang$species.j==i),]
}

Data_Zhang <- SaveZhang


Data_Zhang$Kingdom <- "Annual plant"
Data_Zhang$Model.use <- "APM"
Data_Zhang$Initial.Definition <- "Saavedra"
Data_Zhang$First.Author <- "Zhang"
Data_Zhang$Year <- "2019"
Data_Zhang$lab.OR.field<- "mesocosm"
#######################################################################################################----
# Spaak & De laender, 2020
#######################################################################################################----
Data_Spaak <- data.frame(
ND	= 0.2057,
FD	=0.2057,
r_j =	0.0511668,
fj_0	=0.22729,
S_i	=0.816253208,
r_i	=0.02148,
S_j	=0.774883189,
fi_0 =0.1169,
Kingdom = "Phytoplankton",
species.i = "Phytoplankton",
Model.use = "None",
Initial.Definition = "Spaak",
First.Author = "Spaak",
Year = "2020",
lab.OR.field= "lab"
)

#######################################################################################################----
# All data
#######################################################################################################----

df <-c(Data_Adler,Data_Bim,Data_Lanuz,Data_Car,Data_God,Data_GodLev,
    Data_Hal,Data_Li,Data_Lanuz,
    Data_Narwa,Data_Ocam,Data_Rey,
    Data_Sief,Data_Tan,Data_Ven,
    Data_Vers,Data_Wain,Data_Zhao)

# NFD_meta  is a 'data.frame':	781 obs. of  44 variables:
NFD_meta <- full_join(Data_Bim,Data_Lanuz) #1-2
NFD_meta <- full_join(NFD_meta,Data_Adler) #3
NFD_meta <- full_join(NFD_meta,Data_Petry) #4
NFD_meta <- full_join(NFD_meta,Data_Car) #5
NFD_meta <- full_join(NFD_meta,Data_God) #
NFD_meta <- full_join(NFD_meta,Data_GodLev) #
NFD_meta <- full_join(NFD_meta,Data_Hal) #
NFD_meta <- full_join(NFD_meta,Data_Li) #
NFD_meta <- full_join(NFD_meta,Data_Narwa) 
NFD_meta <- full_join(NFD_meta,Data_Ocam)
NFD_meta <- full_join(NFD_meta,Data_Rey) #
NFD_meta <- full_join(NFD_meta,Data_Sief) #
NFD_meta <- full_join(NFD_meta,Data_Tan) #
NFD_meta <- full_join(NFD_meta,Data_Ven) #
NFD_meta <- full_join(NFD_meta,Data_Wain) #
NFD_meta <- full_join(NFD_meta,Data_Zhao) #
NFD_meta <- full_join(NFD_meta,Data_Vers) #
NFD_meta <- full_join(NFD_meta,Data_Germ) #
NFD_meta <- full_join(NFD_meta,Data_Armi) #
NFD_meta <- full_join(NFD_meta,Data_Ke) #
NFD_meta <- full_join(NFD_meta,Data_Jack ) #
NFD_meta <- full_join(NFD_meta, Data_Jia) #
NFD_meta <- full_join(NFD_meta,Data_Grainger) #
NFD_meta <- full_join(NFD_meta,Data_Matias) #
NFD_meta <- full_join(NFD_meta,Data_Black) #
NFD_meta <- full_join(NFD_meta,Data_Zhang) #
NFD_meta <- full_join(NFD_meta,Data_Gall) 
NFD_meta <- full_join(NFD_meta,Data_Spaak) #


# check the data is full, you should have 38 authors
authors <- as.factor(NFD_meta$First.Author)
levels(as.factor(authors))
str(NFD_meta)
write_csv(NFD_meta,file.path("Data","NFD_meta.csv"), 
          na = "NA", append = F,
          col_names =TRUE)

str(NFD_meta)

NFD_meta.recap <- aggregate(cbind(count = species.i) ~ First.Author + Year + Kingdom + Initial.Definition 
                            + Model.use + lab.OR.field ,
                            NFD_meta, FUN = function(x){NROW(x)})
NFD_meta.recap$unique.com <- NA
for ( n in levels(as.factor(NFD_meta.recap$First.Author))){
  g <- NFD_meta[NFD_meta$First.Author==n,c("species.j","species.i")]
  g <- unique(g)
   NFD_meta.recap$unique.com[NFD_meta.recap$First.Author ==n] <- nrow(g)
}
meta_analysis.adler <- levels(as.factor(Data_Adler$First.Author))
sum.meta_analysis.adler <- sum(NFD_meta.recap[which(NFD_meta.recap$First.Author %in% meta_analysis.adler ),"count"])
sum.com.meta_analysis.adler <- sum(NFD_meta.recap[which(NFD_meta.recap$First.Author %in% meta_analysis.adler ),"unique.com"])

NFD_meta.recap <- NFD_meta.recap[which(!NFD_meta.recap$First.Author %in% levels(as.factor(Data_Adler$First.Author))),]
                                              
NFD_meta.recap <- add_row(NFD_meta.recap,First.Author = "Adler", Kingdom = "Perennial plant",
        Initial.Definition = "Chesson 00", 
        Year=2018,
        count= sum.meta_analysis.adler,
        unique.com = sum.com.meta_analysis.adler)

NFD_meta.recap <- add_row(NFD_meta.recap,First.Author = "Spaak", Kingdom = "Phytoplankton",
                          Initial.Definition = "Spaak", 
                          Year=2021,
                          Model.use= "none",
                          lab.OR.field= "lab",
                          count= 1,
                          unique.com = 1)

NFD_meta.recap <- arrange(NFD_meta.recap, Kingdom)

NFD_meta.recap <- add_row(NFD_meta.recap,First.Author = "Total", 
                          Kingdom = "",
                          Initial.Definition = "", 
                          count= nrow(NFD_meta),
                          unique.com = sum(NFD_meta.recap$unique.com))

write_csv(NFD_meta.recap,file.path("Results","NFD_meta.recap.csv"), 
          na = "NA", append = F,
          col_names =TRUE)

npoint.kd <- data.frame(kingdom= kd,communities=rep(NA,times=4))
for (i in  kd){
  npoint.kd[which(npoint.kd$kingdom==i),"communities"] <- length(NFD_meta$Kingdom[which(NFD_meta$Kingdom == i)])
}
npoint.kd[length(kd)+1,"communities"] <- sum(npoint.kd$communities)
npoint.kd[length(kd)+1,"kingdom"] <- "total"
