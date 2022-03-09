library("tidyverse")
options(dplyr.summarise.inform = FALSE)
library("DescTools")
library("entropy")
library(readr)

setwd("~/project_statistical_test")


detach_package <- function(pkg, character.only = FALSE)
{
  if(!character.only)
  {
    pkg <- deparse(substitute(pkg))
  }
  search_item <- paste("package", pkg, sep = ":")
  while(search_item %in% search())
  {
    detach(search_item, unload = TRUE, character.only = TRUE)
  }
}
detach_package("plyr", TRUE)


# Loading the data set
NFD_clust <- read_csv("Data_for_Javier_22_02_25.csv")

# Filtering the dataset
NFD_clust <- NFD_clust %>% 
  # at least the probability of one process should be positive
  filter(prob_prio + prob_fac + prob_neg_freq>=0) %>% 
  # at least the probability of one cluster should be positive
  filter(low_ND_clust  + high_ND_clust   + var_clust>=0) %>%
  # the growth methods should be either observation, space or time
  filter(Growth_method == "observation" | Growth_method == "space" | Growth_method =="time" )

# I define study as an string joining the author and the year o pubblication
NFD_clust$study <- paste0(NFD_clust$First.Author, '-', NFD_clust$Year) 

# Since a unique study uses the "exponential annual plant model", we join them with
# the studies with "annual plant model"
NFD_clust[NFD_clust$Model.use == "Exponential Annual plant", ]$Model.use = "Annual plant model"
# In the dataset, some rows have a typo regarging the model, so we correct it
NFD_clust[NFD_clust$Model.use == "Annual plant  model", ]$Model.use = "Annual plant model"


NFD_clust <- NFD_clust %>% 
  mutate(ID = row_number(),
         Kingdom = as.factor(Kingdom),
         comp_outcome = as.factor(comp_outcome),
         Model.use = as.factor(Model.use),
         lab.OR.field = as.factor(lab.OR.field),
         Growth_method = as.factor(Growth_method),
         Sympatric = as.factor(Sympatric)
  )

#levels of the different facotrs
kingdoms       <- levels(NFD_clust$Kingdom);           n_kingdoms <- length(kingdoms)
comp_outcomes  <- levels(NFD_clust$comp_outcome);      n_comp_outcomes <- length(comp_outcomes)
clusters       <- c("low_ND", "high_ND", "variation"); n_clusters <- length(clusters)
processes      <- c("fac", "prior", "neg_FD");         n_processes <- length(processes)
Models         <- levels(NFD_clust$Model.use);         n_models <- length(Models)
ExpSettings    <- levels(NFD_clust$lab.OR.field);      n_ExpSettings <- length(ExpSettings)
Growth_methodS <- levels(NFD_clust$Growth_method);     n_Growth_methodS <- length(Growth_methodS)
SympatricS     <- levels(NFD_clust$Sympatric);         n_SympatricS  <- length(SympatricS)

# number of random iteration
n_shuffle = 5E2

freqs0 <- array(rep(0, n_kingdoms*n_clusters*n_comp_outcomes*n_processes*n_models*n_ExpSettings*n_Growth_methodS*n_Growth_methodS), #n_kingdoms*n_clusters*n_outcomes
                dim = c(n_kingdoms, n_clusters, n_comp_outcomes, n_processes, n_models, n_ExpSettings, n_Growth_methodS, n_SympatricS),
                dimnames = list(kingdoms, clusters, comp_outcomes, processes, Models, ExpSettings, Growth_methodS, SympatricS))

Freqs_real <- list(freqs0)[rep(1,n_shuffle)]




for (i in 1:n_shuffle){
  if (i %% 1 == 0){print(paste0(i, '-', Sys.time())) }
  data_i <- NFD_clust
  data_i <- data_i %>% mutate(prob_pos_freq = 1.  - prob_fac - prob_neg_freq)
  
  # We randomly assign each species pair to a cluster and a process, based on the probabilities
  data_i$process_factor = "neg_FD"
  data_i$cluster_factor = "variation"
  
  z_0 <- runif(nrow(data_i), 0, 1)
  z_1 <- runif(nrow(data_i), 0, 1)
  data_i[z_0 <= data_i$prob_pos_freq,]$process_factor = "pos_FD"
  data_i[z_1 <= data_i$prob_fac/(1 - data_i$prob_pos_freq +0.0001) & z_0 > data_i$prob_pos_freq, ]$process_factor = "fac"
  
  z_0 <- runif(nrow(data_i), 0, 1)
  z_1 <- runif(nrow(data_i), 0, 1)
  
  data_i[z_0 <= data_i$low_ND_clust,]$cluster_factor = "low_ND"
  data_i[z_1 <= data_i$high_ND_clust/(1 - data_i$low_ND_clust +0.0001) & z_0 > data_i$low_ND_clust, ]$cluster_factor = "high_ND"
  
  data_i <- data_i %>%
    mutate(cluster_factor = as.factor(cluster_factor),
           process_factor = as.factor(process_factor)
    )
  
  
  # Some species pairs are studied multiple times in some studies. We randomly select one of the rows.
  data_0 <- data_i %>%
    group_by(study,species.i,species.j) %>%
    sample_n(1)
  
  data_s <- data_0 %>% ungroup()
  
  Freqs_real[[i]] <- xtabs(count ~ study + Kingdom + cluster_factor + comp_outcome + process_factor + Model.use + lab.OR.field + Growth_method + Sympatric, 
                           data=data_s %>%
                             group_by(.dots=c("study", "Kingdom", "cluster_factor", "comp_outcome", "process_factor", "Model.use", "lab.OR.field", "Growth_method", "Sympatric"), .drop=FALSE) %>%
                             dplyr::summarise(count= n()))
  
  
}

save.image("AnalysisShuffle_meta_clusters_vfinal.RData")

