library(Matrix)
library(metafor)
library(plyr)

library(tidyverse)
library(readr)
library(ggplot2)

library(extrafont)
loadfonts(device = "win")



setwd("~/project_statistical_test")

# Loading the data
load("~/project_statistical_test/AnalysisShuffle_meta_clusters_vfinal") 


# levels of the different factors
OUTCOMES   = levels(data_0$comp_outcome) 
CLUSTERS   = levels(data_0$cluster_factor) 
KINGDOMS   = levels(data_0$Kingdom) 
STUDIES    = unique(data_0$study)
PROCESSES  = levels(data_0$process_factor) 
MODELS     = levels(data_0$Model.use) 
SETTINGS   = levels(data_0$lab.OR.field) 
GROWTHMETS = levels(data_0$Growth_method) 
SYMPATRICS = levels(data_0$Sympatric) 

# Compute the average frequency among each of the factors for the different shuffles.
Real_mu_freq <- apply(simplify2array(Freqs_real),1:9,mean)


# For each study we store the proportions of the different factors
Freq_study_2 = data.frame()
for (study in STUDIES){
  test = c(study);  names_test = c("study")
  
  # frequency of each kingdom
  for (kingdom in KINGDOMS){
    test = append(test, sum(as.numeric(Real_mu_freq[study, kingdom, , , , , , , ])))
    names_test = append(names_test, kingdom) }
  
  # frequency of each cluster
  for (clust in CLUSTERS){
    test = append(test, sum(as.numeric(Real_mu_freq[study, , clust , , , , , , ])))
    names_test = append(names_test, clust) }
  
  # frequency of each outcome
  for (outcome in OUTCOMES){
    test = append(test, sum(as.numeric(Real_mu_freq[study, , ,outcome,  , , , , ])))
    names_test = append(names_test, outcome) }
  
  # frequency of each process
  for (process in PROCESSES){
    test = append(test, sum(as.numeric(Real_mu_freq[study, , , ,process, ,  , , ])))
    names_test = append(names_test, process)
  }
  
  # frequency of each model
  for (model_i in MODELS){
    test = append(test, sum(as.numeric(Real_mu_freq[study, , , , , model_i,  , , ])))
    names_test = append(names_test, model_i)
  }
  
  # frequency of each setting
  for (setting_i in SETTINGS){
    test = append(test, sum(as.numeric(Real_mu_freq[study, , , , , ,setting_i, , ])))
    names_test = append(names_test, setting_i)
  }
  
  # frequency of each growth methods
  for (gm_i in GROWTHMETS){
    test = append(test, sum(as.numeric(Real_mu_freq[study, , , , , , ,gm_i, ])))
    names_test = append(names_test, gm_i)
  }
  
  # frequency of each sympatric
  for (sympa_i in SYMPATRICS){
    test = append(test, sum(as.numeric(Real_mu_freq[study, , , , , , , ,sympa_i])))
    names_test = append(names_test, sympa_i)
  }
  
  # frequency of each combination clust-outcome
  for (clust in CLUSTERS){
    for (outcome in OUTCOMES){
      test = append(test, sum(as.numeric(Real_mu_freq[study, ,clust ,outcome,  , , , , ])))
      names_test = append(names_test,  paste(clust, outcome, sep="_")) } }
  
  # frequency of each combination clust-process
  for (clust in CLUSTERS){
    for (process in PROCESSES){
      test = append(test, sum(as.numeric(Real_mu_freq[study, ,clust , , process, ,  , , ])))
      names_test = append(names_test,  paste(clust, process, sep="_")) } }
  
  # frequency of each combination kingdom-cluster
  for (kingdom in KINGDOMS){
    for (clust in CLUSTERS){
      test = append(test, sum(as.numeric(Real_mu_freq[study,kingdom ,clust , , , , , , ])))
      names_test = append(names_test,  paste(kingdom, clust, sep="_")) } }
  
  # frequency of each combination kingdom-outcome
  for (kingdom in KINGDOMS){
    for (outcome in OUTCOMES){
      test = append(test, sum(as.numeric(Real_mu_freq[study,kingdom , ,outcome,  , , , , ])))
      names_test = append(names_test,  paste(kingdom, outcome, sep="_")) } }
  
  # frequency of each combination clust-model
  for (clust in CLUSTERS){
    for (model_i in MODELS){
      test = append(test, sum(as.numeric(Real_mu_freq[study, ,clust , , , model_i,  , , ])))
      names_test = append(names_test,  paste(clust, model_i, sep="_")) } }
  
  # frequency of each combination clust-setting
  for (clust in CLUSTERS){
    for (setting_i in SETTINGS){
      test = append(test, sum(as.numeric(Real_mu_freq[study, ,clust , , ,  ,setting_i, , ])))
      names_test = append(names_test,  paste(clust, setting_i, sep="_")) } }
  
  # frequency of each combination clust-gowthmethod
  for (clust in CLUSTERS){
    for (gm_i in GROWTHMETS){
      test = append(test, sum(as.numeric(Real_mu_freq[study, ,clust , , ,  , ,gm_i, ])))
      names_test = append(names_test,  paste(clust, gm_i, sep="_")) } }
  
  # frequency of each combination clust-sympatric
  for (clust in CLUSTERS){
    for (sympa_i in SYMPATRICS){
      test = append(test, sum(as.numeric(Real_mu_freq[study, ,clust , , ,  , , ,sympa_i])))
      names_test = append(names_test,  paste(clust, sympa_i, sep="_")) } }
  
  new_row <- data.frame(t(test)); names(new_row) <- names_test; Freq_study_2 = rbind(Freq_study_2, new_row)
}
rm(new_row, test, kingdom, clust, outcome, study, process, model_i, setting_i, gm_i, sympa_i)

# convert columns to numbers
for (i_col in seq(2,length(names_test))){
  Freq_study_2[[names_test[i_col]]] = as.numeric(Freq_study_2[[names_test[i_col]]]) }
rm(i_col, names_test)

# convert numbers of occurrences of each kingdoms to probabilities (so \sum_i p(x_i) = 1)
Freq_study_2 <- Freq_study_2 %>% 
  mutate(
    frac_AnnualP = `Annual plant` / (`Annual plant` + Bacteria + `Perennial plant` + Phytoplankton),
    frac_Bact    = Bacteria    / (`Annual plant` + Bacteria + `Perennial plant` + Phytoplankton),
    frac_PerennP = `Perennial plant` / (`Annual plant` + Bacteria + `Perennial plant` + Phytoplankton),
    frac_Phyto   = Phytoplankton   / (`Annual plant` + Bacteria + `Perennial plant` + Phytoplankton),
    is_AnnualP = `Annual plant` == pmax(`Annual plant`, Bacteria, `Perennial plant`, Phytoplankton),
    is_Bact    = Bacteria    == pmax(`Annual plant`, Bacteria, `Perennial plant`, Phytoplankton),
    is_PerennP = `Perennial plant` == pmax(`Annual plant`, Bacteria, `Perennial plant`, Phytoplankton),
    is_Phyto   = Phytoplankton   == pmax(Bacteria, `Perennial plant`, Phytoplankton)
    )

# For each study we determine the ecological group (or kingodom) based on the proportions of the groups 
test_kingdoms <- c()
for (study_i in STUDIES){
  Freq_study_2_i <- Freq_study_2 %>% filter(study == study_i)
  test_king = ''
  if (Freq_study_2_i$frac_AnnualP==1){ test_king <- paste(test_king, "Annual plant", sep = '')    }
  if (Freq_study_2_i$frac_Bact==1   ){ test_king <- paste(test_king, "Bacteria", sep = '')        }
  if (Freq_study_2_i$frac_PerennP==1){ test_king <- paste(test_king, "Perennial plant", sep = '') }
  if (Freq_study_2_i$frac_Phyto==1  ){ test_king <- paste(test_king, "Phytoplankton", sep = '')   }
  test_kingdoms <- c(test_kingdoms, test_king)
}
Freq_study_2$max_king = test_kingdoms

# We order the dataset based on the studied ecological group
Freq_study_2 <- arrange(Freq_study_2, max_king)


# count the number of times for each study that the outcome does not occur
# and the fraction of species pairs with different outcomes, processes, and clusters
Freq_study_2 <- Freq_study_2 %>% mutate(number_data = low_ND + high_ND + variation)
Freq_study_2 <- Freq_study_2 %>% mutate(
  No_Coex = number_data - Coexistence, 
  No_Prio = number_data - `Priority effect`, 
  No_Excl  = number_data - `Competitive exclusion`,
  
  frac_Coex =   Coexistence/(Coexistence + `Competitive exclusion` + `Priority effect`),
  frac_Excl =   `Competitive exclusion`/(Coexistence + `Competitive exclusion` + `Priority effect`),
  frac_Prior = `Priority effect`/(Coexistence + `Competitive exclusion` + `Priority effect`),
  
  frac_lowND =  low_ND/(low_ND + high_ND + variation),
  frac_highND = high_ND/(low_ND + high_ND + variation),
  frac_variation = variation/(low_ND + high_ND + variation),
  
  frac_negFD = neg_FD/(neg_FD + fac + pos_FD),
  frac_fac = fac/(neg_FD + fac + pos_FD),
  frac_posFD = pos_FD/(neg_FD + fac + pos_FD)
  )


# based on the fraction of occurrences, we select the population model...
Freq_study_2$model = ""
Freq_study_2[Freq_study_2$`Annual plant model` > 0,]$model = "Annual plant model"
Freq_study_2[Freq_study_2$`Lotka Volterra` > 0,]$model = "Lotka Volterra"
Freq_study_2[Freq_study_2$None > 0,]$model = "None"

# ...experimental setting
Freq_study_2$setting = ""
Freq_study_2[Freq_study_2$lab > 0,]$setting = "lab"
Freq_study_2[Freq_study_2$Field > 0,]$setting = "Field"
Freq_study_2[Freq_study_2$Greenhouse  > 0,]$setting = "Greenhouse "

# ... growth method
Freq_study_2$growthmethod = ""
Freq_study_2[Freq_study_2$observation > 0,]$growthmethod = "observation"
Freq_study_2[Freq_study_2$space > 0,]$growthmethod = "space"
Freq_study_2[Freq_study_2$time > 0,]$growthmethod = "time"

# ... sympatric
Freq_study_2$sympatric = "False"
Freq_study_2[Freq_study_2$`TRUE` > 0,]$sympatric = "True"

# ...we convert all of them to factors
Freq_study_2 <- Freq_study_2 %>% mutate(
  model = as.factor(model), setting = as.factor(setting),
  growthmethod = as.factor(growthmethod), sympatric = as.factor(sympatric),
  max_king = as.factor(max_king)
)


# Performing the meta-analyses


VT = "LS" # type of sampling variances to calculate. "LS" (default: large sample approximation). 
MEAS = "PR" #which effect size or outcome measure should be calculated. "PR" for the raw proportion,


##### CLUSTER low ND
dat_lowND <- escalc(
  measure = MEAS, xi = get("low_ND"), ni = get("number_data")*1.0001,
  data = Freq_study_2, add = 0.5, to = "only0", drop00 = TRUE, vtype=VT,
  var.names=c("y", "v"),
  add.measure=TRUE, append=TRUE, replace=TRUE, digits=4)

# low_ND vs kingdom 
res_lowND <- rma(yi = y, vi = v, mods = ~ factor(max_king), method = "ML", data = dat_lowND)
predictions <- predict(res_lowND, newmods=cbind(c(0,1,0,0),c(0,0,1,0),c(0,0,0,1)), addx=TRUE)
predictions$king <- levels(Freq_study_2$max_king)
predictions$king1 <- levels(Freq_study_2$max_king)
predictions$clust <- c("low n. d.", "low n. d.", "low n. d.", "low n. d.")
predictions$outc <- c("any", "any", "any", "any")
predictions_lowND_kingdom <- data.frame(predictions) %>% select(pred, se, ci.lb, ci.ub, cr.lb, cr.ub, king, king1, clust, outc)

# low_ND vs model
res_lowND <- rma(yi = y, vi = v, mods = ~ factor(model), method = "ML", data = dat_lowND)
predictions <- predict(res_lowND, newmods=cbind(c(0,1,0),c(0,0,1)), addx=TRUE)
predictions$model <- levels(Freq_study_2$model)
predictions$clust <- c("low n. d.", "low n. d.", "low n. d.")
predictions$outc <- c("any", "any", "any")
predictions_lowND_model <- data.frame(predictions) %>% select(pred, se, ci.lb, ci.ub, cr.lb, cr.ub, model, clust, outc)

# low_ND vs setting
res_lowND <- rma(yi = y, vi = v, mods = ~ factor(setting), method = "ML", data = dat_lowND)
predictions <- predict(res_lowND, newmods=cbind(c(0,1,0),c(0,0,1)), addx=TRUE)
predictions$setting <- levels(Freq_study_2$setting)
predictions$clust <- c("low n. d.", "low n. d.", "low n. d.")
predictions$outc <- c("any", "any", "any")
predictions_lowND_setting <- data.frame(predictions) %>% select(pred, se, ci.lb, ci.ub, cr.lb, cr.ub, setting, clust, outc)

# low_ND vs growth method
res_lowND <- rma(yi = y, vi = v, mods = ~ factor(growthmethod), method = "ML", data = dat_lowND)
predictions <- predict(res_lowND, newmods=cbind(c(0,1,0),c(0,0,1)), addx=TRUE)
predictions$growthmethod <- levels(Freq_study_2$growthmethod)
predictions$clust <- c("low n. d.", "low n. d.", "low n. d.")
predictions$outc <- c("any", "any", "any")
predictions_lowND_growthmeth <- data.frame(predictions) %>% select(pred, se, ci.lb, ci.ub, cr.lb, cr.ub, growthmethod, clust, outc)

# low_ND vs sympatric
res_lowND <- rma(yi = y, vi = v, mods = ~ factor(sympatric), method = "ML", data = dat_lowND)
predictions <- predict(res_lowND, newmods=c(0,1), addx=TRUE)
predictions$sympatric <- levels(Freq_study_2$sympatric)
predictions$clust <- c("low n. d.", "low n. d.")
predictions$outc <- c("any", "any")
predictions_lowND_sympa <- data.frame(predictions) %>% select(pred, se, ci.lb, ci.ub, cr.lb, cr.ub, sympatric, clust, outc)

# low_ND alone
res_lowND <- rma(yi = y, vi = v, method = "ML", data = dat_lowND)
predictions <- predict(res_lowND, addx=TRUE)
predictions$clust <- c("low n. d.")
predictions_lowND_nothing <- data.frame(predictions) %>% select(pred, se, ci.lb, ci.ub, cr.lb, cr.ub, clust)


##### CLUSTER high ND
dat_highND <- escalc(
  measure = MEAS, xi = get("high_ND"), ni = get("number_data")*1.0001,
  data = Freq_study_2, add = 0.5, to = "only0", drop00 = TRUE, vtype=VT,
  var.names=c("y", "v"),
  add.measure=TRUE, append=TRUE, replace=TRUE, digits=4)

# high_ND vs kingdom 
res_highND <- rma(yi = y, vi = v, mods = ~ factor(max_king), method = "ML", data = dat_highND)
predictions <- predict(res_highND, newmods=cbind(c(0,1,0,0),c(0,0,1,0),c(0,0,0,1)), addx=TRUE)
predictions$king <- levels(Freq_study_2$max_king)
predictions$king1 <- levels(Freq_study_2$max_king)
predictions$clust <- c("high n. d.", "high n. d.", "high n. d.", "high n. d.")
predictions$outc <- c("any", "any", "any", "any")
predictions_highND_kingdom <- data.frame(predictions) %>% select(pred, se, ci.lb, ci.ub, cr.lb, cr.ub, king, king1, clust, outc)

# high_ND vs model
res_highND <- rma(yi = y, vi = v, mods = ~ factor(model), method = "ML", data = dat_highND)
predictions <- predict(res_highND, newmods=cbind(c(0,1,0),c(0,0,1)), addx=TRUE)
predictions$model <- levels(Freq_study_2$model)
predictions$clust <- c("high n. d.", "high n. d.", "high n. d.")
predictions$outc <- c("any", "any", "any")
predictions_highND_model <- data.frame(predictions) %>% select(pred, se, ci.lb, ci.ub, cr.lb, cr.ub, model, clust, outc)

# high_ND vs setting
res_highND <- rma(yi = y, vi = v, mods = ~ factor(setting), method = "ML", data = dat_highND)
predictions <- predict(res_highND, newmods=cbind(c(0,1,0),c(0,0,1)), addx=TRUE)
predictions$setting <- levels(Freq_study_2$setting)
predictions$clust <- c("high n. d.", "high n. d.", "high n. d.")
predictions$outc <- c("any", "any", "any")
predictions_highND_setting <- data.frame(predictions) %>% select(pred, se, ci.lb, ci.ub, cr.lb, cr.ub, setting, clust, outc)

# high_ND vs growth method
res_highND <- rma(yi = y, vi = v, mods = ~ factor(growthmethod), method = "ML", data = dat_highND)
predictions <- predict(res_highND, newmods=cbind(c(0,1,0),c(0,0,1)), addx=TRUE)
predictions$growthmethod <- levels(Freq_study_2$growthmethod)
predictions$clust <- c("high n. d.", "high n. d.", "high n. d.")
predictions$outc <- c("any", "any", "any")
predictions_highND_growthmeth <- data.frame(predictions) %>% select(pred, se, ci.lb, ci.ub, cr.lb, cr.ub, growthmethod, clust, outc)

# high_ND vs sympatric
res_highND <- rma(yi = y, vi = v, mods = ~ factor(sympatric), method = "ML", data = dat_highND)
predictions <- predict(res_highND, newmods=c(0,1), addx=TRUE)
predictions$sympatric <- levels(Freq_study_2$sympatric)
predictions$clust <- c("high n. d.", "high n. d.")
predictions$outc <- c("any", "any")
predictions_highND_sympa <- data.frame(predictions) %>% select(pred, se, ci.lb, ci.ub, cr.lb, cr.ub, sympatric, clust, outc)

# high ND alone
res_highND <- rma(yi = y, vi = v, method = "ML", data = dat_highND)
predictions <- predict(res_highND, addx=TRUE)
predictions$clust <- c("high n. d.")
predictions_highND_nothing <- data.frame(predictions) %>% select(pred, se, ci.lb, ci.ub, cr.lb, cr.ub, clust)


##### CLUSTER variation
dat_var <- escalc(
  measure = MEAS, xi = get("variation"), ni = get("number_data")*1.0001,
  data = Freq_study_2, add = 0.5, to = "only0", drop00 = TRUE, vtype=VT,
  var.names=c("y", "v"),
  add.measure=TRUE, append=TRUE, replace=TRUE, digits=4)

# variation vs kingdom 
res_varND <- rma(yi = y, vi = v, mods = ~ factor(max_king), method = "ML", data = dat_var)
predictions <- predict(res_varND, newmods=cbind(c(0,1,0,0),c(0,0,1,0),c(0,0,0,1)), addx=TRUE)
predictions$king <- levels(Freq_study_2$max_king)
predictions$king1 <- levels(Freq_study_2$max_king)
predictions$clust <- c("variation", "variation", "variation", "variation")
predictions$outc <- c("any", "any", "any", "any")
predictions_var_kingdom <- data.frame(predictions) %>% select(pred, se, ci.lb, ci.ub, cr.lb, cr.ub, king, king1, clust, outc)

# variation vs model
res_varND <- rma(yi = y, vi = v, mods = ~ factor(model), method = "ML", data = dat_var)
predictions <- predict(res_varND, newmods=cbind(c(0,1,0),c(0,0,1)), addx=TRUE)
predictions$model <- levels(Freq_study_2$model)
predictions$clust <- c("variation", "variation", "variation")
predictions$outc <- c("any", "any", "any")
predictions_var_model <- data.frame(predictions) %>% select(pred, se, ci.lb, ci.ub, cr.lb, cr.ub, model, clust, outc)

# variation vs setting
res_varND <- rma(yi = y, vi = v, mods = ~ factor(setting), method = "ML", data = dat_var)
predictions <- predict(res_varND, newmods=cbind(c(0,1,0),c(0,0,1)), addx=TRUE)
predictions$setting <- levels(Freq_study_2$setting)
predictions$clust <- c("variation", "variation", "variation")
predictions$outc <- c("any", "any", "any")
predictions_var_setting <- data.frame(predictions) %>% select(pred, se, ci.lb, ci.ub, cr.lb, cr.ub, setting, clust, outc)

# variation vs growth method
res_varND <- rma(yi = y, vi = v, mods = ~ factor(growthmethod), method = "ML", data = dat_var)
predictions <- predict(res_varND, newmods=cbind(c(0,1,0),c(0,0,1)), addx=TRUE)
predictions$growthmethod <- levels(Freq_study_2$growthmethod)
predictions$clust <- c("variation", "variation", "variation")
predictions$outc <- c("any", "any", "any")
predictions_var_growthmeth <- data.frame(predictions) %>% select(pred, se, ci.lb, ci.ub, cr.lb, cr.ub, growthmethod, clust, outc)

# variation vs sympatric
res_varND <- rma(yi = y, vi = v, mods = ~ factor(sympatric), method = "ML", data = dat_var)
predictions <- predict(res_varND, newmods=c(0,1), addx=TRUE)
predictions$sympatric <- levels(Freq_study_2$sympatric)
predictions$clust <- c("variation", "variation")
predictions$outc <- c("any", "any")
predictions_var_sympa <- data.frame(predictions) %>% select(pred, se, ci.lb, ci.ub, cr.lb, cr.ub, sympatric, clust, outc)

# variation alone
res_varND <- rma(yi = y, vi = v, method = "ML", data = dat_var)
predictions <- predict(res_varND, addx=TRUE)
predictions$clust <- c("variation")
predictions_var_nothing <- data.frame(predictions) %>% select(pred, se, ci.lb, ci.ub, cr.lb, cr.ub,  clust)


### PREDICTIONS: CLUSTER (KINGDOM, MODEL, SETTING)
predictions_kingdom <- rbind(predictions_lowND_kingdom, predictions_highND_kingdom, predictions_var_kingdom)
predictions_setting <- rbind(predictions_lowND_setting, predictions_highND_setting, predictions_var_setting)
predictions_model   <- rbind(predictions_lowND_model,   predictions_highND_model,   predictions_var_model)
predictions_growthmeth   <- rbind(predictions_lowND_growthmeth,   predictions_highND_growthmeth,   predictions_var_growthmeth)
predictions_sympa   <- rbind(predictions_lowND_sympa,   predictions_highND_sympa,   predictions_var_sympa)
predictions_nothing   <- rbind(predictions_lowND_nothing,   predictions_highND_nothing,   predictions_var_nothing)


predictions_kingdom <- predictions_kingdom %>% mutate(king = as.factor(king), king1 = as.factor(king1), clust = as.factor(clust))
predictions_setting <- predictions_setting %>% mutate(setting = as.factor(setting), clust = as.factor(clust))
predictions_model <- predictions_model %>% mutate(model = as.factor(model), clust = as.factor(clust))
predictions_growthmeth <- predictions_growthmeth %>% mutate(growthmethod = as.factor(growthmethod), clust = as.factor(clust))
predictions_sympa <- predictions_sympa %>% mutate(sympatric = as.factor(sympatric), clust = as.factor(clust))



# All individual studies together

dat_lowND_2 <- dat_lowND
dat_lowND_2$king1 = dat_lowND_2$max_king
dat_lowND_2$pred = dat_lowND_2$y
dat_lowND_2$clust <- "low n. d."

dat_highND_2 <- dat_highND
dat_highND_2$king1 = dat_highND_2$max_king
dat_highND_2$pred = dat_highND_2$y
dat_highND_2$clust <- "high n. d."

dat_var_2 <- dat_var
dat_var_2$king1 = dat_var_2$max_king
dat_var_2$pred = dat_var_2$y
dat_var_2$clust <- "variation"


dat_all <- rbind(dat_lowND_2, dat_highND_2, dat_var_2)
dat_all <- dat_all %>% mutate(clust = factor(clust))

# renaming the clusters
dat_all$clust <- mapvalues(dat_all$clust, from = c("high n. d.", "low n. d.", "variation"),
                                to = c("High ND", "Low ND", "Variable"))
predictions_kingdom$clust <- mapvalues(predictions_kingdom$clust, from = c("high n. d.", "low n. d.", "variation"),
                                       to = c("High ND", "Low ND", "Variable"))
predictions_setting$clust <- mapvalues(predictions_setting$clust, from = c("high n. d.", "low n. d.", "variation"),
                                       to = c("High ND", "Low ND", "Variable"))
predictions_model$clust <- mapvalues(predictions_model$clust, from = c("high n. d.", "low n. d.", "variation"),
                                       to = c("High ND", "Low ND", "Variable"))
predictions_growthmeth$clust <- mapvalues(predictions_growthmeth$clust, from = c("high n. d.", "low n. d.", "variation"),
                                     to = c("High ND", "Low ND", "Variable"))
predictions_sympa$clust <- mapvalues(predictions_sympa$clust, from = c("high n. d.", "low n. d.", "variation"),
                                     to = c("High ND", "Low ND", "Variable"))
predictions_nothing$clust <- mapvalues(predictions_nothing$clust, from = c("high n. d.", "low n. d.", "variation"),
                                       to = c("High ND", "Low ND", "Variable"))
# reorder the clusters
dat_all <- dat_all %>% mutate(clust = factor(clust, levels = c("Variable", "Low ND", "High ND")))
predictions_kingdom <- predictions_kingdom %>% mutate(clust = factor(clust, levels = c("Variable", "Low ND", "High ND")))
predictions_setting <- predictions_setting %>% mutate(clust = factor(clust, levels = c("Variable", "Low ND", "High ND")))
predictions_model <- predictions_model %>% mutate(clust = factor(clust, levels = c("Variable", "Low ND", "High ND")))
predictions_growthmeth <- predictions_growthmeth %>% mutate(clust = factor(clust, levels = c("Variable", "Low ND", "High ND")))
predictions_sympa <- predictions_sympa %>% mutate(clust = factor(clust, levels = c("Variable", "Low ND", "High ND")))
predictions_nothing <- predictions_nothing %>% mutate(clust = factor(clust, levels = c("Variable", "Low ND", "High ND")))

# # renaming the clusters
# dat_all$clust <- mapvalues(dat_all$clust, from = c("high n. d.", "low n. d.", "variation"), 
#                                 to = c("A", "B", "C"))
# predictions_kingdom$clust <- mapvalues(predictions_kingdom$clust, from = c("high n. d.", "low n. d.", "variation"), 
#                                        to = c("A", "B", "C"))
# 
# predictions_setting$clust <- mapvalues(predictions_setting$clust, from = c("high n. d.", "low n. d.", "variation"), 
#                                        to = c("A", "B", "C"))
# 
# predictions_model$clust <- mapvalues(predictions_model$clust, from = c("high n. d.", "low n. d.", "variation"), 
#                                        to = c("A", "B", "C"))
# predictions_growthmeth$clust <- mapvalues(predictions_growthmeth$clust, from = c("high n. d.", "low n. d.", "variation"), 
#                                      to = c("A", "B", "C"))
# predictions_sympa$clust <- mapvalues(predictions_sympa$clust, from = c("high n. d.", "low n. d.", "variation"), 
#                                      to = c("A", "B", "C"))
# predictions_nothing$clust <- mapvalues(predictions_nothing$clust, from = c("high n. d.", "low n. d.", "variation"), 
#                                        to = c("A", "B", "C"))
# 
# # reorder the clusters
# dat_all <- dat_all %>% mutate(clust = factor(clust, levels = c("C", "B", "A")))
# predictions_kingdom <- predictions_kingdom %>% mutate(clust = factor(clust, levels = c("C", "B", "A")))
# predictions_setting <- predictions_setting %>% mutate(clust = factor(clust, levels = c("C", "B", "A")))
# predictions_model <- predictions_model %>% mutate(clust = factor(clust, levels = c("C", "B", "A")))
# predictions_growthmeth <- predictions_growthmeth %>% mutate(clust = factor(clust, levels = c("C", "B", "A")))
# predictions_sympa <- predictions_sympa %>% mutate(clust = factor(clust, levels = c("C", "B", "A")))
# predictions_nothing <- predictions_nothing %>% mutate(clust = factor(clust, levels = c("C", "B", "A")))

# Renaming factors
# ... kingdom
dat_all$king1 <- mapvalues(
  dat_all$king1, from = c("Annual plant", "Bacteria", "Perennial plant", "Phytoplankton"), 
  to = c("Annual plant", "Bacteria/yeast", "Perennial plant", "Phytoplankton"))
predictions_kingdom$king1 <- mapvalues(
  predictions_kingdom$king1, from = c("Annual plant", "Bacteria", "Perennial plant", "Phytoplankton"), 
  to = c("Annual plant", "Bacteria/yeast", "Perennial plant", "Phytoplankton"))
# ... experimental setting
dat_all$setting <- mapvalues(
  dat_all$setting, from = c("Field", "Greenhouse ", "lab"), 
  to = c("Field", "Greenhouse", "Lab"))
predictions_setting$setting <- mapvalues(
  predictions_setting$setting, from = c("Field", "Greenhouse ", "lab"), 
  to = c("Field", "Greenhouse", "Lab"))
# ... growth method
dat_all$growthmethod <- mapvalues(
  dat_all$growthmethod, from = c("observation", "space", "time"), 
  to = c("Observation", "Space", "Time"))
predictions_growthmeth$growthmethod <- mapvalues(
  predictions_growthmeth$growthmethod, from = c("observation", "space", "time"), 
  to = c("Observation", "Space", "Time"))
# ... population model
dat_all$model <- mapvalues(
  dat_all$model, from = c("Annual plant model", "Lotka Volterra", "None"), 
  to = c("Annual pl. mod.", "Lotka Volterra", "None"))
predictions_model$model <- mapvalues(
  predictions_model$model, from = c("Annual plant model", "Lotka Volterra", "None"), 
  to = c("Annual pl. mod.", "Lotka Volterra", "None"))



# plot panel A: proportion of the clusters
plotA <- ggplot() +
  geom_point(data = dat_all,
             aes(y = as.numeric(factor(clust)) - 0.15 * ( 0.5 + 1), x=pred, size = number_data),
             alpha = I(0.2)) +
  geom_point(data = predictions_nothing,
             aes(y = as.numeric(factor(clust)) , x=pred),
             size = I(3)) +
  geom_errorbarh(data = predictions_nothing,
                 aes(y = as.numeric(factor(clust)) , x=pred, xmin=ci.lb, xmax=ci.ub),
                 height = 0, size=1.25)  +
  labs(x = "proportion", y = "cluster",
       color = "comm. type", shape = "comm. type", fill = "comm. type", size = "species pairs") + 
  
  scale_size_continuous(#trans = "log2",
    breaks = c(1, 50, 100),
    limits = c(1, 100),
    range = c(0.1, 5)
  ) +
  scale_color_manual(values = c("#999999", "#E69F00", "#56B4E9", "#009E73")) +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9", "#009E73")) +
  scale_shape_manual(values=c(24,21,22,25)) + 
  scale_y_continuous(breaks = 1:3, labels=levels(as.factor(dat_all$clust))) +
  coord_cartesian(xlim=c(0,1), ylim = c(0.7, 3.225)) +
  theme_bw()+
  theme(
    axis.text.y = element_text(size=10, family="Arial"),
    axis.text.x = element_text(size=10, family="Arial"),
    axis.title = element_text(size=12,color="black", family="Arial"),
    legend.title = element_text(size=11, color="black", family="Arial"),
    legend.text = element_text(size=9, family="Arial")) 

# panel B: proportion of clusters vs ecological group
plotB <- ggplot() +
  geom_point(data = dat_all,
             aes(y = as.numeric(factor(clust)) - 0.15 * (as.numeric(factor(king1)) - 2.5 + 0.5), 
                 x=pred, color=factor(king1), size = number_data),
             alpha = I(0.2)) +
  geom_point(data = predictions_kingdom,
             aes(y = as.numeric(factor(clust)) - 0.15 * (as.numeric(factor(king1)) - 2.5), x=pred,
                 color=factor(king1), shape = factor(king1), fill = factor(king1)),
             size = I(3)) +
  geom_errorbarh(data = predictions_kingdom,
                 aes(y = as.numeric(factor(clust)) - 0.15 * (as.numeric(factor(king1)) - 2.5), x=pred,
                     xmin=ci.lb, xmax=ci.ub, color=factor(king1)),
                 height = 0, size=1.25)  +
  labs(x = "proportion", y = "cluster",
       color = "ecol. group", shape = "ecol. group", fill = "ecol. group", size = "species pairs") + 
  
  scale_size_continuous(#trans = "log2",
    breaks = c(1, 50, 100),
    limits = c(1, 100),
    range = c(0.1, 5)
  ) +
  scale_color_manual(values = c("#323232", "#E69F00", "#56B4E9", "#009E73")) +
  scale_fill_manual(values = c("#323232", "#E69F00", "#56B4E9", "#009E73")) +
  scale_shape_manual(values=c(24,21,22,25)) + 
  scale_y_continuous(breaks = 1:3, labels=levels(as.factor(dat_all$clust))) +
  coord_cartesian(xlim=c(0,1), ylim = c(0.7, 3.225)) +
  theme_bw()+
  theme(
    axis.text.y = element_text(size=10, family="Arial"),
    axis.text.x = element_text(size=10, family="Arial"),
    axis.title = element_text(size=12,color="black", family="Arial"),
    legend.title = element_text(size=11, color="black", family="Arial"),
    legend.text = element_text(size=9, family="Arial")) 

# panel C: proportion of clusters vs sympatric
plotC <- ggplot() +
  geom_point(data = dat_all,
             aes(y = as.numeric(factor(clust)) - 0.15 * (2*(as.numeric(factor(sympatric)) - 1.5) + 0.5), 
                 x=pred, color=factor(sympatric), size = number_data),
             alpha = I(0.2)) +
  geom_point(data = predictions_sympa,
             aes(y = as.numeric(factor(clust)) - 0.15 * (2*(as.numeric(factor(sympatric)) - 1.5)), x=pred, #3
                 color=factor(sympatric), shape = factor(sympatric), fill = factor(sympatric)),
             size = I(3)) +
  geom_errorbarh(data = predictions_sympa,
                 aes(y = as.numeric(factor(clust)) - 0.15 * (2*(as.numeric(factor(sympatric)) - 1.5)), x=pred,
                     xmin=ci.lb, xmax=ci.ub, color=factor(sympatric)),
                 height = 0, size=1.25)  +
  labs(x = "proportion", y = "cluster",
       color = "sympatric", shape = "sympatric", fill = "sympatric", size = "species pairs") + 
  scale_size_continuous(#trans = "log2",
    breaks = c(1, 50, 100),
    limits = c(1, 100),
    range = c(0.1, 5)
  ) +
  scale_color_manual(values = c("#BB5566",  "#004488")) + 
  scale_fill_manual(values = c("#BB5566", "#004488")) +
  scale_shape_manual(values=c(24,21,22)) + #25
  scale_y_continuous(breaks = 1:3, labels=levels(as.factor(dat_all$clust))) +
  coord_cartesian(xlim=c(0,1),  ylim = c(0.7, 3.225)) +
  theme_bw()+
  theme(
    axis.text.y = element_text(size=10, family="Arial"),
    axis.text.x = element_text(size=10, family="Arial"),
    axis.title = element_text(size=12,color="black", family="Arial"),
    legend.title = element_text(size=11, color="black", family="Arial"),
    legend.text = element_text(size=9, family="Arial")) 



# panel D: proportion of clusters vs experimental setting
plotD <- ggplot() +
  geom_point(data = dat_all,
             aes(y = as.numeric(factor(clust)) - 0.15 * (1.5*(as.numeric(factor(setting)) - 2) + 0.5), 
                 x=pred, color=factor(setting), size = number_data),
             alpha = I(0.2)) +
  geom_point(data = predictions_setting,
             aes(y = as.numeric(factor(clust)) - 0.15 * (1.5*(as.numeric(factor(setting)) - 2)), x=pred,
                 color=factor(setting), shape = factor(setting), fill = factor(setting)),
             size = I(3)) +
  geom_errorbarh(data = predictions_setting,
                 aes(y = as.numeric(factor(clust)) - 0.15 * (1.5*(as.numeric(factor(setting)) - 2)), x=pred,
                     xmin=ci.lb, xmax=ci.ub, color=factor(setting)),
                 height = 0, size=1.25)  +
  labs(x = "proportion", y = "cluster",
       color = "exp. setting", shape = "exp. setting", fill = "exp. setting", size = "species pairs") + 
  scale_size_continuous(#trans = "log2",
    breaks = c(1, 50, 100),
    limits = c(1, 100),
    range = c(0.1, 5)
  ) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  scale_shape_manual(values=c(24,21,22)) + 
  scale_y_continuous(breaks = 1:3, labels=levels(as.factor(dat_all$clust))) +
  coord_cartesian(xlim=c(0,1), ylim = c(0.7, 3.225)) +
  theme_bw()+
  theme(
    axis.text.y = element_text(size=10, family="Arial"),
    axis.text.x = element_text(size=10, family="Arial"),
    axis.title = element_text(size=12,color="black", family="Arial"),
    legend.title = element_text(size=11, color="black", family="Arial"),
    legend.text = element_text(size=9, family="Arial")) 



# panel E: proportion of clusters vs population model
plotE <- ggplot() +
  geom_point(data = dat_all,
             aes(y = as.numeric(factor(clust)) - 0.15 * (1.5*(as.numeric(factor(model)) - 2) + 0.5), 
                 x=pred, color=factor(model), size = number_data),
             alpha = I(0.2)) +
  geom_point(data = predictions_model,
             aes(y = as.numeric(factor(clust)) - 0.15 * (1.5*(as.numeric(factor(model)) - 2)), x=pred,
                 color=factor(model), shape = factor(model), fill = factor(model)),
             size = I(3)) +
  geom_errorbarh(data = predictions_model,
                 aes(y = as.numeric(factor(clust)) - 0.15 * (1.5*(as.numeric(factor(model)) - 2)), x=pred,
                     xmin=ci.lb, xmax=ci.ub, color=factor(model)),
                 height = 0, size=1.25)  +
  labs(x = "proportion", y = "cluster",
       color = "pop. model", shape = "pop. model", fill = "pop. model", size = "species pairs") + 
  scale_size_continuous(#trans = "log2",
    breaks = c(1, 50, 100),
    limits = c(1, 100),
    range = c(0.1, 5)
  ) +
  scale_color_manual(values = c("#0072B2", "#D55E00", "#CC79A7")) + 
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#CC79A7")) +
  scale_shape_manual(values=c(24,21,22)) + #25
  scale_y_continuous(breaks = 1:3, labels=levels(as.factor(dat_all$clust))) +
  coord_cartesian(xlim=c(0,1), ylim = c(0.7, 3.225)) +
  theme_bw()+
  theme(
    axis.text.y = element_text(size=10, family="Arial"),
    axis.text.x = element_text(size=10, family="Arial"),
    axis.title = element_text(size=12,color="black", family="Arial"),
    legend.title = element_text(size=11, color="black", family="Arial"),
    legend.text = element_text(size=9, family="Arial")) 

# panel F: proportion of clusters vs growth method
plotF <- ggplot() +
  geom_point(data = dat_all,
             aes(y = as.numeric(factor(clust)) - 0.15 * (1.5*(as.numeric(factor(growthmethod)) - 2) + 0.5),
                 x=pred, color=factor(growthmethod), size = number_data),
             alpha = I(0.2)) +
  geom_point(data = predictions_growthmeth,
             aes(y = as.numeric(factor(clust)) - 0.15 * (1.5*(as.numeric(factor(growthmethod)) - 2)), x=pred,
                 color=factor(growthmethod), shape = factor(growthmethod), fill = factor(growthmethod)),
             size = I(3)) +
  geom_errorbarh(data = predictions_growthmeth,
                 aes(y = as.numeric(factor(clust)) - 0.15 * (1.5*(as.numeric(factor(growthmethod)) - 2)), x=pred,
                     xmin=ci.lb, xmax=ci.ub, color=factor(growthmethod)),
                 height = 0, size=1.25)  +
  labs(x = "proportion", y = "cluster",
       color = "growth meth.", shape = "growth meth.", fill = "growth meth.", size = "species pairs") + 
  scale_size_continuous(#trans = "log2",
    breaks = c(1, 50, 100),
    limits = c(1, 100),
    range = c(0.1, 5)
  ) +
  scale_color_manual(values = c("#DDAA33", "#BB5566", "#004488")) +  
  scale_fill_manual(values = c("#DDAA33", "#BB5566", "#004488")) +
  scale_shape_manual(values=c(24,21,22)) + 
  scale_y_continuous(breaks = 1:3, labels=levels(as.factor(dat_all$clust))) +
  coord_cartesian(xlim=c(0,1), ylim = c(0.7, 3.225)) +
  theme_bw()+
  theme(
    axis.text.y = element_text(size=10, family="Arial"),
    axis.text.x = element_text(size=10, family="Arial"),
    axis.title = element_text(size=12,color="black", family="Arial"),
    legend.title = element_text(size=11, color="black", family="Arial"),
    legend.text = element_text(size=9, family="Arial")) 






library(ggpubr)

ggarrange(plotA + theme(legend.box = "horizontal") +
            guides(size = guide_legend(order=1)
            ) +
            scale_y_continuous(limits = layer_scales(plotE)$y$range$range, 
                               breaks = 1:3, labels=levels(as.factor(dat_all$clust))),
          plotB +
            scale_size_continuous(breaks = NaN, 
                                  limits = c(1, 100), range = c(0.1, 5) ) +
            theme(legend.box = "horizontal") +
            guides(color = guide_legend(order=1),
                   shape = guide_legend(order=1),
                   fill = guide_legend(order=1)
            ), 
          plotC +
            scale_size_continuous(breaks = NaN, 
                                  limits = c(1, 100), range = c(0.1, 5) ) +
            theme(legend.box = "horizontal") +
            guides(color = guide_legend(order=1),
                   shape = guide_legend(order=1),
                   fill = guide_legend(order=1)
            ), 
          plotD +
            scale_size_continuous(breaks = NaN, 
                                  limits = c(1, 100), range = c(0.1, 5) ) +
            theme(legend.box = "horizontal") +
            guides(color = guide_legend(order=1),
                   shape = guide_legend(order=1),
                   fill = guide_legend(order=1)
            ),
          plotE +
            scale_size_continuous(breaks = NaN, 
                                  limits = c(1, 100), range = c(0.1, 5) ) +
            theme(legend.box = "horizontal") +
            guides(color = guide_legend(order=1),
                   shape = guide_legend(order=1),
                   fill = guide_legend(order=1)
            ),
          plotF +
            scale_size_continuous(breaks = NaN, 
                                  limits = c(1, 100), range = c(0.1, 5) ) +
            theme(legend.box = "horizontal") +
            guides(color = guide_legend(order=1),
                   shape = guide_legend(order=1),
                   fill = guide_legend(order=1)
            ),
          #common.legend = FALSE, legend = "right",
          labels = c("A", "B","C","D", "E","F"),
          ncol = 2, nrow = 3,align = "hv", widths = c(1, 1, 1),
          font.label = list(size = 16, family = "Arial", color = "black"),
          heights = c(1, 1, 1, 1, 1 ,1 ))

ggsave("meta_clusters_final.svg", width = 18, height = 16, units = "cm")



