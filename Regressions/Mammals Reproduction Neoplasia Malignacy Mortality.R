rm(list=ls()) #remove all objects into memmory
library(readxl)
library(ape)
library(nlme)
library(geiger)
library(ggeffects)
library(caper)
library(phylolm)
library(phytools)
library(MuMIn)
library(ggplot2)
library(rr2)
library(DHARMa)
library(dplyr)
library(tidyr)
library(ggpubr)
library(grid)
library(brms)
library(parallel)
library(tidybayes)
library(bayesplot)
library(Hmisc)
library(bayestestR)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Import phylogeny data ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#import the phylogenetic tree for mammals
setwd("C:/Users/antoi/My Drive/Postdoc/Zoo data/Biases/Phylogenies")
phylo_all <- read.tree("phylogeny_time_tree.tre") 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Import Neoplasia and malignancy data ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("C:/Users/antoi/My Drive/Postdoc/Zoo data/Biases/Scripts/Estimation of confounding effects/for github") #set working directory
data_neo_mal <- read_excel("Raw datasets.xlsx", sheet = "Mammal_reproduction_neo_mal")
rownames(data_neo_mal) <- data_neo_mal$Species

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import Mortality data----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

data_mortality <- read_excel("Raw datasets.xlsx", sheet = "Mammal_reproduction_mortality")
rownames(data_mortality) <-data_mortality$Species

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Calculate median for bilbiometric variabels to use in calculation or marinal effects ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#calculate the median of confounding variables to use in the calculation of marginal effects
median_logHindex <- median((c(data_neo_mal$logHIndex,data_mortality$logHIndex)))
median_logNInstitutions <- median((c(data_neo_mal$logNInstitutions, data_mortality$logNInstitutions)))
median_logNExSitu <- median((c(data_neo_mal$logNExSitu,  data_mortality$logNExSitu)))
median_ResidualslogNViews <- median((c(data_neo_mal$ResidualslogNViews,  data_mortality$ResidualslogNViews)))
median_logNExSitu <- median((c(data_neo_mal$logNExSitu,  data_mortality$logNExSitu)))
median_logMass <- median((c(data_neo_mal$logMass,  data_mortality$logMass)))
median_logMaxLongevity <- median((c(data_neo_mal$logMaxLongevity,  data_mortality$logMaxLongevity)))
median_logLactationLength  <- median((c(data_neo_mal$logLactationLength ,  data_mortality$logLactationLength)))
median_logGestationLength   <- median((c(data_neo_mal$logGestationLength  ,  data_mortality$logGestationLength )))
median_Littersize  <- median((c(data_neo_mal$Littersize  ,  data_mortality$Littersize )))
median_mean_GDP_5years  <- median((c(data_neo_mal$mean_GDP_5years  ,  data_mortality$mean_GDP_5years )))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Prepare var covar matrixes for brms models ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#check which species are both in the main dataframe and the phylogeny tree
obj_brms_neo_mal<- name.check(phylo_all, data_neo_mal)
obj_brms_neo_mal

#drop the species from the tree for which we have no cancer data
phylo_neo_mal <- drop.tip(phylo_all, obj_brms_neo_mal$tree_not_data)

#compute the var-covar matrix from the tree to use in brms
vcv_brms_OU_neo_mal <- vcv(corMartins(1, phylo_neo_mal, fixed=FALSE,form=~phylo_neo_mal$tip.label)) #OU model

#check which species are both in the main dataframe and the phylogeny tree
obj_brms_mortality<- name.check(phylo_all, data_mortality)
obj_brms_mortality

#drop the species from the tree for which we have no cancer data
phylo_mortality <- drop.tip(phylo_all, obj_brms_mortality$tree_not_data)

#compute the var-covar matrix from the tree to use in brms
vcv_brms_OU_mortality <- vcv(corMartins(1, phylo_mortality , fixed=FALSE,form=~phylo_mortality$tip.label)) #OU model

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Create a function to calculate the marginal effects for the crude and adjusted models
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

get_marginal_effects <- function(crude_model, adjusted_model, resp, effect, varext, NNecropsies, conditions_crude, conditions_adjusted){
  
  #calculate marginal effect for neoplasia for unadjusted model
  marginal_effects_crude <- conditional_effects(crude_model , conditions = conditions_crude, effect = effect, resp = resp)
  marginal_effects_crude <- data.frame(marginal_effects_crude[[1]])
  marginal_effects_crude$estimate__ <- marginal_effects_crude$estimate__/marginal_effects_crude$NNecropsies[1]
  marginal_effects_crude$lower__ <- marginal_effects_crude$lower__/marginal_effects_crude$NNecropsies[1]
  marginal_effects_crude$upper__ <- marginal_effects_crude$upper__/marginal_effects_crude$NNecropsies[1]
  marginal_effects_crude <- marginal_effects_crude[,c(varext, "estimate__", "lower__", "upper__")]
  marginal_effects_crude$Type <- "Crude"
  
  #calculate marginal effect for neoplasia for adjusted model
  conditions <- conditions_adjusted
  marginal_effects_adjusted <- conditional_effects(adjusted_model , conditions = conditions, effect = effect, resp = resp)
  marginal_effects_adjusted <- data.frame(marginal_effects_adjusted[[1]])
  marginal_effects_adjusted$estimate__ <- marginal_effects_adjusted$estimate__/marginal_effects_adjusted$NNecropsies[1]
  marginal_effects_adjusted$lower__ <- marginal_effects_adjusted$lower__/marginal_effects_adjusted$NNecropsies[1]
  marginal_effects_adjusted$upper__ <- marginal_effects_adjusted$upper__/marginal_effects_adjusted$NNecropsies[1]
  marginal_effects_adjusted <- marginal_effects_adjusted[,c(varext, "estimate__", "lower__", "upper__")]
  marginal_effects_adjusted$Type <- "Adjusted"
  
  #return the combined marginal effects
  marginal_effects_combined <- rbind(marginal_effects_crude, marginal_effects_adjusted)
  marginal_effects_combined$Type <- factor(marginal_effects_combined$Type, levels = c("Crude", "Adjusted"))
  return(marginal_effects_combined)
  
}  

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Create a function to help with model selection
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

diagnostic <- function(reference_model, adjusted_model){
  
  #extract fixed effect for the reference model
  fe_reference <- data.frame(fixef(reference_model)[,-2] )
  fe_reference$Significance <- as.numeric(sign(fe_reference[,2])+sign(fe_reference[,3]))
  fe_reference$Significance <- ifelse(fe_reference$Significance == 0, "", "**S**")
  
  #get the p-direction for the reference model
  pdir_ref <- data.frame(p_direction(reference_model))
  fe_reference$pdir <- pdir_ref[,2] 
  
  #get the fix effect for the adjusted model
  fe_adjusted <- data.frame(fixef(adjusted_model)[,-2]) 
  fe_adjusted$Significance <- as.numeric(sign(fe_adjusted[,2])+sign(fe_adjusted[,3]))
  fe_adjusted$Significance <- ifelse(fe_adjusted$Significance == 0, "", "**S**")
  fe_adjusted$PercentageChange <- rep("", nrow(fe_adjusted))
  
  #get the p-direction for the adjusted model
  pdir_adj <- data.frame(p_direction(adjusted_model))
  fe_adjusted$pdir <- pdir_adj[,2] 
  
  #calculate the percentage changes for the adjusted model when compared to the reference one
  for(i in 1:nrow(fe_adjusted)){
    sel <- which(rownames(fe_adjusted)[i] == rownames(fe_reference))
    if(length(sel) > 0){
      
      fe_adjusted$PercentageChange[i] <- round(((fe_adjusted[i,1] / fe_reference[sel,1]) - 1)  * 100, digits = 1) 
      
    }
  }
  
  list(Reference = fe_reference , Adjusted = fe_adjusted)
  
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Neoplasia ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Step 1: Univariate models ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#fit univariate model for logMass
fit_reproduction_neo_logMass <- brm(NNeoplasia | trials(NNecropsies)  ~ logMass + (1|gr(Species, cov = vcv_brms_OU)),
                                              family =binomial(),
                                              data = data_neo_mal,
                                              data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                              chains = 4,
                                              cores = 24,
                                              threads = threading(2),
                                              iter = 10000, warmup = 5000,
                                              thin = 10,
                                              seed = 1,
                                              control=list(max_treedepth = 12, adapt_delta = 0.99),
                                              save_pars = save_pars(all = TRUE))
fit_reproduction_neo_logMass
p_direction(fit_reproduction_neo_logMass)
#p_direction >0.90 -> retain for multivariate model -

#fit univariate model for logMaxLongevity
fit_reproduction_neo_logMaxLongevity <- brm(NNeoplasia | trials(NNecropsies)  ~ logMaxLongevity + (1|gr(Species, cov = vcv_brms_OU)),
                                    family =binomial(),
                                    data = data_neo_mal,
                                    data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                    chains = 4,
                                    cores = 24,
                                    threads = threading(2),
                                    iter = 10000, warmup = 5000,
                                    thin = 10,
                                    seed = 1,
                                    control=list(max_treedepth = 12, adapt_delta = 0.99),
                                    save_pars = save_pars(all = TRUE))
fit_reproduction_neo_logMaxLongevity
p_direction(fit_reproduction_neo_logMaxLongevity)
#p_direction <0.90 -> drop from the multivariate model -

#fit univariate model for LactationLength
fit_reproduction_neo_logLactationLength <- brm(NNeoplasia | trials(NNecropsies)  ~ logLactationLength + (1|gr(Species, cov = vcv_brms_OU)),
                                            family =binomial(),
                                            data = data_neo_mal,
                                            data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                            chains = 4,
                                            cores = 24,
                                            threads = threading(2),
                                            iter = 10000, warmup = 5000,
                                            thin = 10,
                                            seed = 1,
                                            control=list(max_treedepth = 12, adapt_delta = 0.99),
                                            save_pars = save_pars(all = TRUE))
fit_reproduction_neo_logLactationLength
p_direction(fit_reproduction_neo_logLactationLength)
#p_direction <0.90 -> exclude from the model

#fit univariate model for logGestationLength
fit_reproduction_neo_logGestationLength <- brm(NNeoplasia | trials(NNecropsies)  ~ logGestationLength + (1|gr(Species, cov = vcv_brms_OU)),
                                            family =binomial(),
                                            data = data_neo_mal,
                                            data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                            chains = 4,
                                            cores = 24,
                                            threads = threading(2),
                                            iter = 10000, warmup = 5000,
                                            thin = 10,
                                            seed = 1,
                                            control=list(max_treedepth = 12, adapt_delta = 0.99),
                                            save_pars = save_pars(all = TRUE))
fit_reproduction_neo_logGestationLength
p_direction(fit_reproduction_neo_logGestationLength)
#p_direction >0.90 -> retain for multivariate model -

#fit univariate model for placentation type
fit_reproduction_neo_Placentation <- brm(NNeoplasia | trials(NNecropsies)  ~ Placentation + (1|gr(Species, cov = vcv_brms_OU)),
                                               family =binomial(),
                                               data = data_neo_mal,
                                               data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                               chains = 4,
                                               cores = 24,
                                               threads = threading(2),
                                               iter = 10000, warmup = 5000,
                                               thin = 10,
                                               seed = 1,
                                               control=list(max_treedepth = 12, adapt_delta = 0.99),
                                               save_pars = save_pars(all = TRUE))
fit_reproduction_neo_Placentation
p_direction(fit_reproduction_neo_Placentation)
#p_direction >0.90 -> retain for multivariate model -

#fit univariate model for DietVertebrate 
fit_reproduction_neo_DietVertebrate <- brm(NNeoplasia | trials(NNecropsies)  ~ DietVertebrate + (1|gr(Species, cov = vcv_brms_OU)),
                                         family =binomial(),
                                         data = data_neo_mal,
                                         data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                         chains = 4,
                                         cores = 24,
                                         threads = threading(2),
                                         iter = 10000, warmup = 5000,
                                         thin = 10,
                                         seed = 1,
                                         control=list(max_treedepth = 12, adapt_delta = 0.99),
                                         save_pars = save_pars(all = TRUE))
fit_reproduction_neo_DietVertebrate
p_direction(fit_reproduction_neo_DietVertebrate)
#p_direction >0.90 ->  retain for multivariate model -

#fit univariate model for LitterSize 
fit_reproduction_neo_Littersize <- brm(NNeoplasia | trials(NNecropsies)  ~ Littersize + (1|gr(Species, cov = vcv_brms_OU)),
                                           family =binomial(),
                                           data = data_neo_mal,
                                           data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                           chains = 4,
                                           cores = 24,
                                           threads = threading(2),
                                           iter = 10000, warmup = 5000,
                                           thin = 10,
                                           seed = 1,
                                           control=list(max_treedepth = 12, adapt_delta = 0.99),
                                           save_pars = save_pars(all = TRUE))
fit_reproduction_neo_Littersize
p_direction(fit_reproduction_neo_Littersize)
#p_direction >0.90 ->  retain for multivariate model -

#fit univariate model for logHIndex
fit_reproduction_neo_logHIndex <- brm(NNeoplasia | trials(NNecropsies)  ~ logHIndex + (1|gr(Species, cov = vcv_brms_OU)),
                                       family =binomial(),
                                       data = data_neo_mal,
                                       data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                       chains = 4,
                                       cores = 24,
                                       threads = threading(2),
                                       iter = 10000, warmup = 5000,
                                       thin = 10,
                                       seed = 1,
                                       control=list(max_treedepth = 12, adapt_delta = 0.99),
                                       save_pars = save_pars(all = TRUE))
fit_reproduction_neo_logHIndex
p_direction(fit_reproduction_neo_logHIndex)
#p_direction >0.90 ->  retain for multivariate model -


#fit univariate model for mean_GDP_5years
fit_reproduction_neo_mean_GDP_5years <- brm(NNeoplasia | trials(NNecropsies)  ~ mean_GDP_5years + (1|gr(Species, cov = vcv_brms_OU)),
                                      family =binomial(),
                                      data = data_neo_mal,
                                      data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                      chains = 4,
                                      cores = 24,
                                      threads = threading(2),
                                      iter = 10000, warmup = 5000,
                                      thin = 10,
                                      seed = 1,
                                      control=list(max_treedepth = 12, adapt_delta = 0.99),
                                      save_pars = save_pars(all = TRUE))
fit_reproduction_neo_mean_GDP_5years
p_direction(fit_reproduction_neo_mean_GDP_5years)
#p_direction >0.90 ->  retain for multivariate model -

#fit univariate model for logNInstitutions
fit_reproduction_neo_logNInstitutions <- brm(NNeoplasia | trials(NNecropsies)  ~logNInstitutions + (1|gr(Species, cov = vcv_brms_OU)),
                                            family =binomial(),
                                            data = data_neo_mal,
                                            data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                            chains = 4,
                                            cores = 24,
                                            threads = threading(2),
                                            iter = 10000, warmup = 5000,
                                            thin = 10,
                                            seed = 1,
                                            control=list(max_treedepth = 12, adapt_delta = 0.99),
                                            save_pars = save_pars(all = TRUE))
fit_reproduction_neo_logNInstitutions
p_direction(fit_reproduction_neo_logNInstitutions)
#p_direction >0.90 ->  retain for multivariate model -

#fit univariate model for logNExSitu 
fit_reproduction_neo_logNExSitu  <- brm(NNeoplasia | trials(NNecropsies)  ~logNExSitu  + (1|gr(Species, cov = vcv_brms_OU)),
                                             family =binomial(),
                                             data = data_neo_mal,
                                             data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                             chains = 4,
                                             cores = 24,
                                             threads = threading(2),
                                             iter = 10000, warmup = 5000,
                                             thin = 10,
                                             seed = 1,
                                             control=list(max_treedepth = 12, adapt_delta = 0.99),
                                             save_pars = save_pars(all = TRUE))
fit_reproduction_neo_logNExSitu 
p_direction(fit_reproduction_neo_logNExSitu )
#p_direction >0.90 ->  retain for multivariate model -

#fit univariate model for logNExSitu 
fit_reproduction_neo_ResidualslogNViews  <- brm(NNeoplasia | trials(NNecropsies)  ~ResidualslogNViews  + (1|gr(Species, cov = vcv_brms_OU)),
                                        family =binomial(),
                                        data = data_neo_mal,
                                        data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                        chains = 4,
                                        cores = 24,
                                        threads = threading(2),
                                        iter = 10000, warmup = 5000,
                                        thin = 10,
                                        seed = 1,
                                        control=list(max_treedepth = 12, adapt_delta = 0.99),
                                        save_pars = save_pars(all = TRUE))
fit_reproduction_neo_ResidualslogNViews
p_direction(fit_reproduction_neo_ResidualslogNViews)
#p_direction >0.90 ->  retain for multivariate model -


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Step 2 : Drop variables part of the model selection process ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#fit a model with all the variables selected in step 1
fit_reproduction_neo_adjusted_m1 <- brm(NNeoplasia | trials(NNecropsies)  ~logMass + DietVertebrate + Littersize  + logGestationLength  + Placentation  + logNInstitutions + logHIndex + mean_GDP_5years + ResidualslogNViews + (1|gr(Species, cov = vcv_brms_OU)),
                                            family =binomial(),
                                            data = data_neo_mal,
                                            data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                            chains = 4,
                                            cores = 24,
                                            iter = 10000, warmup = 5000,
                                            thin = 10,
                                            seed = 1,
                                            control=list(max_treedepth = 12, adapt_delta = 0.99),
                                            save_pars = save_pars(all = TRUE))

fit_reproduction_neo_adjusted_m1 
summary(fit_reproduction_neo_adjusted_m1 , prob = 0.90)
p_direction(fit_reproduction_neo_adjusted_m1) 
#Littersize  has the lowest p-direction

#drop the Littersize variable
fit_reproduction_neo_adjusted_m2 <- brm(NNeoplasia | trials(NNecropsies)  ~logMass + DietVertebrate  + logGestationLength  + Placentation  + logNInstitutions + logHIndex + mean_GDP_5years + ResidualslogNViews + (1|gr(Species, cov = vcv_brms_OU)),
                                            family =binomial(),
                                            data = data_neo_mal,
                                            data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                            chains = 4,
                                            cores = 24,
                                            iter = 10000, warmup = 5000,
                                            thin = 10,
                                            seed = 1,
                                            control=list(max_treedepth = 12, adapt_delta = 0.99),
                                            save_pars = save_pars(all = TRUE))

fit_reproduction_neo_adjusted_m2
diagnostic(reference_model =fit_reproduction_neo_adjusted_m1, adjusted_model = fit_reproduction_neo_adjusted_m2)
p_direction(fit_reproduction_neo_adjusted_m2) 
#littersize doesn't change slopes of the other parameters by >25% -> drop

#drop logMass from the model
fit_reproduction_neo_adjusted_m3 <- brm(NNeoplasia | trials(NNecropsies)  ~DietVertebrate + logGestationLength  + Placentation  + logNInstitutions + logHIndex + mean_GDP_5years + ResidualslogNViews + (1|gr(Species, cov = vcv_brms_OU)),
                                        family =binomial(),
                                        data = data_neo_mal,
                                        data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                        chains = 4,
                                        cores = 24,
                                        iter = 10000, warmup = 5000,
                                        thin = 10,
                                        seed = 1,
                                        control=list(max_treedepth = 12, adapt_delta = 0.99),
                                        save_pars = save_pars(all = TRUE))

fit_reproduction_neo_adjusted_m3
diagnostic(reference_model =fit_reproduction_neo_adjusted_m2, adjusted_model = fit_reproduction_neo_adjusted_m3)
p_direction(fit_reproduction_neo_adjusted_m3) 
#logMass change the slope of the parameters >25% - keep 
#drop logHIndex next

#drop logHIndex
fit_reproduction_neo_adjusted_m4 <- brm(NNeoplasia | trials(NNecropsies)  ~ logMass + DietVertebrate + logGestationLength  + Placentation  + logNInstitutions + mean_GDP_5years + ResidualslogNViews +  (1|gr(Species, cov = vcv_brms_OU)),
                                        family =binomial(),
                                        data = data_neo_mal,
                                        data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                        chains = 4,
                                        cores = 24,
                                        iter = 10000, warmup = 5000,
                                        thin = 10,
                                        seed = 1,
                                        control=list(max_treedepth = 12, adapt_delta = 0.99),
                                        save_pars = save_pars(all = TRUE))

fit_reproduction_neo_adjusted_m4
diagnostic(reference_model =fit_reproduction_neo_adjusted_m2, adjusted_model = fit_reproduction_neo_adjusted_m4)
p_direction(fit_reproduction_neo_adjusted_m4) 
#dropping logmaxlongevity change the slpoe of logGestationLength >25% -> keep
#drop ResidualslogNViews next

#drop ResidualslogNViews
fit_reproduction_neo_adjusted_m5 <- brm(NNeoplasia | trials(NNecropsies)  ~logMass + DietVertebrate  + logGestationLength  + Placentation  + logNInstitutions + logHIndex + mean_GDP_5years + (1|gr(Species, cov = vcv_brms_OU)),
                                        family =binomial(),
                                        data = data_neo_mal,
                                        data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                        chains = 4,
                                        cores = 24,
                                        iter = 10000, warmup = 5000,
                                        thin = 10,
                                        seed = 1,
                                        control=list(max_treedepth = 12, adapt_delta = 0.99),
                                        save_pars = save_pars(all = TRUE))

fit_reproduction_neo_adjusted_m5
diagnostic(reference_model =fit_reproduction_neo_adjusted_m2, adjusted_model = fit_reproduction_neo_adjusted_m5)
p_direction(fit_reproduction_neo_adjusted_m5)
#ResidualslogNViews changes the slopes >25% -> keep

#drop the mean_GDP_5years variable
fit_reproduction_neo_adjusted_m6 <- brm(NNeoplasia | trials(NNecropsies)  ~logMass + DietVertebrate + logGestationLength  + Placentation  + logNInstitutions + logHIndex + ResidualslogNViews + (1|gr(Species, cov = vcv_brms_OU)),
                                        family =binomial(),
                                        data = data_neo_mal,
                                        data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                        chains = 4,
                                        cores = 24,
                                        iter = 10000, warmup = 5000,
                                        thin = 10,
                                        seed = 1,
                                        control=list(max_treedepth = 12, adapt_delta = 0.99),
                                        save_pars = save_pars(all = TRUE))

fit_reproduction_neo_adjusted_m6
diagnostic(reference_model =fit_reproduction_neo_adjusted_m2, adjusted_model = fit_reproduction_neo_adjusted_m6)
p_direction(fit_reproduction_neo_adjusted_m6) 
#mean_GDP_5yearschange the slope  >25% -> keep
#drop DietVertebrate next

#drop DietVertebrate
fit_reproduction_neo_adjusted_m7 <- brm(NNeoplasia | trials(NNecropsies)  ~logMass + logGestationLength  + Placentation  + logNInstitutions + logHIndex + mean_GDP_5years + ResidualslogNViews + (1|gr(Species, cov = vcv_brms_OU)),
                                        family =binomial(),
                                        data = data_neo_mal,
                                        data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                        chains = 4,
                                        cores = 24,
                                        iter = 10000, warmup = 5000,
                                        thin = 10,
                                        seed = 1,
                                        control=list(max_treedepth = 12, adapt_delta = 0.99),
                                        save_pars = save_pars(all = TRUE))

fit_reproduction_neo_adjusted_m7
diagnostic(reference_model =fit_reproduction_neo_adjusted_m2, adjusted_model = fit_reproduction_neo_adjusted_m7)
p_direction(fit_reproduction_neo_adjusted_m7) 
#dropping the DietVertebrate >25% -> keep
#drop logNInstitutions 

#drop logNInstitutions   
fit_reproduction_neo_adjusted_m8 <- brm(NNeoplasia | trials(NNecropsies)  ~logMass + DietVertebrate + logGestationLength  + Placentation  + logHIndex + mean_GDP_5years + ResidualslogNViews +(1|gr(Species, cov = vcv_brms_OU)),
                                        family =binomial(),
                                        data = data_neo_mal,
                                        data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                        chains = 4,
                                        cores = 24,
                                        iter = 10000, warmup = 5000,
                                        thin = 10,
                                        seed = 1,
                                        control=list(max_treedepth = 12, adapt_delta = 0.99),
                                        save_pars = save_pars(all = TRUE))

fit_reproduction_neo_adjusted_m8
diagnostic(reference_model =fit_reproduction_neo_adjusted_m2, adjusted_model = fit_reproduction_neo_adjusted_m8)
p_direction(fit_reproduction_neo_adjusted_m8) 
#drop logNInstitutions   changes the slopes  >25% -> keep

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Step 3 Final model ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#add back logLactationLength to the model
fit_reproduction_neo_adjusted_m9 <- brm(NNeoplasia | trials(NNecropsies)  ~logMass + DietVertebrate +logLactationLength +  logGestationLength  + Placentation  + logNInstitutions + logHIndex + mean_GDP_5years + ResidualslogNViews + (1|gr(Species, cov = vcv_brms_OU)),
                                        family =binomial(),
                                        data = data_neo_mal,
                                        data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                        chains = 4,
                                        cores = 24,
                                        iter = 10000, warmup = 5000,
                                        thin = 10,
                                        seed = 1,
                                        control=list(max_treedepth = 12, adapt_delta = 0.99),
                                        save_pars = save_pars(all = TRUE))

fit_reproduction_neo_adjusted_m9
diagnostic(reference_model =fit_reproduction_neo_adjusted_m2, adjusted_model = fit_reproduction_neo_adjusted_m9)
#adding back logLactationLength to the model change the slope of the parameters >25% 

#add back log max longevity
fit_reproduction_neo_adjusted_m10 <- brm(NNeoplasia | trials(NNecropsies)  ~ logMass + logMaxLongevity + DietVertebrate +logLactationLength +  logGestationLength  + Placentation  + logNInstitutions + logHIndex + mean_GDP_5years + ResidualslogNViews + (1|gr(Species, cov = vcv_brms_OU)),
                                        family =binomial(),
                                        data = data_neo_mal,
                                        data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                        chains = 4,
                                        cores = 24,
                                        iter = 10000, warmup = 5000,
                                        thin = 10,
                                        seed = 1,
                                        control=list(max_treedepth = 12, adapt_delta = 0.99),
                                        save_pars = save_pars(all = TRUE))

fit_reproduction_neo_adjusted_m10
diagnostic(reference_model =fit_reproduction_neo_adjusted_m9, adjusted_model = fit_reproduction_neo_adjusted_m10)
#adding back logMaxLongevity to the model change the slope of the parameters >25% 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Calculate and plot the marginal effects ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#crude model
crude_model_neo <- brm(NNeoplasia | trials(NNecropsies)  ~ logGestationLength  + Placentation + (1|gr(Species, cov = vcv_brms_OU)),
                       family =binomial(),
                       data = data_neo_mal,
                       data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                       chains = 4,
                       cores = 24,
                       threads = threading(2),
                       iter = 10000, warmup = 5000,
                       thin = 10,
                       seed = 1,
                       control=list(max_treedepth = 12, adapt_delta = 0.99),
                       save_pars = save_pars(all = TRUE))
crude_model_neo

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# NMalignancy ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Step 1: Univariate models ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#fit univariate model for logMass
fit_reproduction_mal_logMass <- brm(NMalignant | trials(NNecropsies)  ~ logMass + (1|gr(Species, cov = vcv_brms_OU)),
                                    family =binomial(),
                                    data = data_neo_mal,
                                    data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                    chains = 4,
                                    cores = 24,
                                    threads = threading(2),
                                    iter = 10000, warmup = 5000,
                                    thin = 10,
                                    seed = 1,
                                    control=list(max_treedepth = 12, adapt_delta = 0.99),
                                    save_pars = save_pars(all = TRUE))
fit_reproduction_mal_logMass
p_direction(fit_reproduction_mal_logMass)
#p_direction <0.90 -> exclude from the model

#fit univariate model for logMaxLongevity
fit_reproduction_mal_logMaxLongevity <- brm(NMalignant | trials(NNecropsies)  ~ logMaxLongevity + (1|gr(Species, cov = vcv_brms_OU)),
                                            family =binomial(),
                                            data = data_neo_mal,
                                            data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                            chains = 4,
                                            cores = 24,
                                            threads = threading(2),
                                            iter = 10000, warmup = 5000,
                                            thin = 10,
                                            seed = 1,
                                            control=list(max_treedepth = 12, adapt_delta = 0.99),
                                            save_pars = save_pars(all = TRUE))
fit_reproduction_mal_logMaxLongevity
p_direction(fit_reproduction_mal_logMaxLongevity)
#p_direction <0.90 -> exclude from the model -

#fit univariate model for LactationLength
fit_reproduction_mal_logLactationLength <- brm(NMalignant | trials(NNecropsies)  ~ logLactationLength + (1|gr(Species, cov = vcv_brms_OU)),
                                               family =binomial(),
                                               data = data_neo_mal,
                                               data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                               chains = 4,
                                               cores = 24,
                                               threads = threading(2),
                                               iter = 10000, warmup = 5000,
                                               thin = 10,
                                               seed = 1,
                                               control=list(max_treedepth = 12, adapt_delta = 0.99),
                                               save_pars = save_pars(all = TRUE))
fit_reproduction_mal_logLactationLength
p_direction(fit_reproduction_mal_logLactationLength)
#p_direction <0.90 -> exclude from the model

#fit univariate model for logGestationLength
fit_reproduction_mal_logGestationLength <- brm(NMalignant | trials(NNecropsies)  ~ logGestationLength + (1|gr(Species, cov = vcv_brms_OU)),
                                               family =binomial(),
                                               data = data_neo_mal,
                                               data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                               chains = 4,
                                               cores = 24,
                                               threads = threading(2),
                                               iter = 10000, warmup = 5000,
                                               thin = 10,
                                               seed = 1,
                                               control=list(max_treedepth = 12, adapt_delta = 0.99),
                                               save_pars = save_pars(all = TRUE))
fit_reproduction_mal_logGestationLength
p_direction(fit_reproduction_mal_logGestationLength)
#p_direction >0.90 -> retain for multivariate model -

#fit univariate model for placentation type
fit_reproduction_mal_Placentation <- brm(NMalignant | trials(NNecropsies)  ~ Placentation + (1|gr(Species, cov = vcv_brms_OU)),
                                         family =binomial(),
                                         data = data_neo_mal,
                                         data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                         chains = 4,
                                         cores = 24,
                                         threads = threading(2),
                                         iter = 10000, warmup = 5000,
                                         thin = 10,
                                         seed = 1,
                                         control=list(max_treedepth = 12, adapt_delta = 0.99),
                                         save_pars = save_pars(all = TRUE))
fit_reproduction_mal_Placentation
p_direction(fit_reproduction_mal_Placentation)
#p_direction >0.90 -> retain for multivariate model -

#fit univariate model for DietVertebrate 
fit_reproduction_mal_DietVertebrate <- brm(NMalignant | trials(NNecropsies)  ~ DietVertebrate + (1|gr(Species, cov = vcv_brms_OU)),
                                           family =binomial(),
                                           data = data_neo_mal,
                                           data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                           chains = 4,
                                           cores = 24,
                                           threads = threading(2),
                                           iter = 10000, warmup = 5000,
                                           thin = 10,
                                           seed = 1,
                                           control=list(max_treedepth = 12, adapt_delta = 0.99),
                                           save_pars = save_pars(all = TRUE))
fit_reproduction_mal_DietVertebrate
p_direction(fit_reproduction_mal_DietVertebrate)
#p_direction >0.90 ->  retain for multivariate model -

#fit univariate model for LitterSize 
fit_reproduction_mal_Littersize <- brm(NMalignant | trials(NNecropsies)  ~ Littersize + (1|gr(Species, cov = vcv_brms_OU)),
                                       family =binomial(),
                                       data = data_neo_mal,
                                       data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                       chains = 4,
                                       cores = 24,
                                       threads = threading(2),
                                       iter = 10000, warmup = 5000,
                                       thin = 10,
                                       seed = 1,
                                       control=list(max_treedepth = 12, adapt_delta = 0.99),
                                       save_pars = save_pars(all = TRUE))
fit_reproduction_mal_Littersize
p_direction(fit_reproduction_mal_Littersize)
#p_direction >0.90 ->  retain for multivariate model -

#fit univariate model for logHIndex
fit_reproduction_mal_logHIndex <- brm(NMalignant | trials(NNecropsies)  ~ logHIndex + (1|gr(Species, cov = vcv_brms_OU)),
                                      family =binomial(),
                                      data = data_neo_mal,
                                      data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                      chains = 4,
                                      cores = 24,
                                      threads = threading(2),
                                      iter = 10000, warmup = 5000,
                                      thin = 10,
                                      seed = 1,
                                      control=list(max_treedepth = 12, adapt_delta = 0.99),
                                      save_pars = save_pars(all = TRUE))
fit_reproduction_mal_logHIndex
p_direction(fit_reproduction_mal_logHIndex)
#p_direction >0.90 ->  retain for multivariate model -


#fit univariate model for mean_GDP_5years
fit_reproduction_mal_mean_GDP_5years <- brm(NMalignant | trials(NNecropsies)  ~ mean_GDP_5years + (1|gr(Species, cov = vcv_brms_OU)),
                                            family =binomial(),
                                            data = data_neo_mal,
                                            data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                            chains = 4,
                                            cores = 24,
                                            threads = threading(2),
                                            iter = 10000, warmup = 5000,
                                            thin = 10,
                                            seed = 1,
                                            control=list(max_treedepth = 12, adapt_delta = 0.99),
                                            save_pars = save_pars(all = TRUE))
fit_reproduction_mal_mean_GDP_5years
p_direction(fit_reproduction_mal_mean_GDP_5years)
#p_direction >0.90 ->  retain for multivariate model -

#fit univariate model for logNInstitutions
fit_reproduction_mal_logNInstitutions <- brm(NMalignant | trials(NNecropsies)  ~logNInstitutions + (1|gr(Species, cov = vcv_brms_OU)),
                                             family =binomial(),
                                             data = data_neo_mal,
                                             data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                             chains = 4,
                                             cores = 24,
                                             threads = threading(2),
                                             iter = 10000, warmup = 5000,
                                             thin = 10,
                                             seed = 1,
                                             control=list(max_treedepth = 12, adapt_delta = 0.99),
                                             save_pars = save_pars(all = TRUE))
fit_reproduction_mal_logNInstitutions
p_direction(fit_reproduction_mal_logNInstitutions)
#p_direction >0.90 ->  retain for multivariate model -

#fit univariate model for logNExSitu 
fit_reproduction_mal_ResidualslogNViews  <- brm(NMalignant | trials(NNecropsies)  ~logNExSitu  + (1|gr(Species, cov = vcv_brms_OU)),
                                                family =binomial(),
                                                data = data_neo_mal,
                                                data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                                chains = 4,
                                                cores = 24,
                                                threads = threading(2),
                                                iter = 10000, warmup = 5000,
                                                thin = 10,
                                                seed = 1,
                                                control=list(max_treedepth = 12, adapt_delta = 0.99),
                                                save_pars = save_pars(all = TRUE))
fit_reproduction_mal_ResidualslogNViews
p_direction(fit_reproduction_mal_ResidualslogNViews)
#p_direction >0.90 ->  retain for multivariate model -

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Step 2 : Drop variables part of the model selection process ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#fit a model with all the variables selected in step 1
fit_reproduction_mal_adjusted_m1 <- brm(NMalignant | trials(NNecropsies)  ~DietVertebrate + Littersize  + logGestationLength  + Placentation  + logNInstitutions + logHIndex + mean_GDP_5years + ResidualslogNViews + (1|gr(Species, cov = vcv_brms_OU)),
                                        family =binomial(),
                                        data = data_neo_mal,
                                        data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                        chains = 4,
                                        cores = 24,
                                        iter = 10000, warmup = 5000,
                                        thin = 10,
                                        seed = 1,
                                        control=list(max_treedepth = 12, adapt_delta = 0.99),
                                        save_pars = save_pars(all = TRUE))

fit_reproduction_mal_adjusted_m1 
summary(fit_reproduction_mal_adjusted_m1 , prob = 0.90)
p_direction(fit_reproduction_mal_adjusted_m1) 
#Littersize  has the lowest p-direction -> drop first

#model with Littersize  dropped
fit_reproduction_mal_adjusted_m2 <- brm(NMalignant | trials(NNecropsies)  ~ DietVertebrate   + logGestationLength  + Placentation  + logNInstitutions + logHIndex + mean_GDP_5years + ResidualslogNViews + (1|gr(Species, cov = vcv_brms_OU)),
                                        family =binomial(),
                                        data = data_neo_mal,
                                        data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                        chains = 4,
                                        cores = 24,
                                        iter = 10000, warmup = 5000,
                                        thin = 10,
                                        seed = 1,
                                        control=list(max_treedepth = 12, adapt_delta = 0.99),
                                        save_pars = save_pars(all = TRUE))

fit_reproduction_mal_adjusted_m2 
diagnostic(reference_model =fit_reproduction_mal_adjusted_m1, adjusted_model = fit_reproduction_mal_adjusted_m2)
#keep in the model (changes some slopes >25%)
#drop logHIndex  

#model with logHIndex  dropped
fit_reproduction_mal_adjusted_m3 <- brm(NMalignant | trials(NNecropsies)  ~ DietVertebrate + Littersize  + logGestationLength  + Placentation  + logNInstitutions + mean_GDP_5years + ResidualslogNViews + (1|gr(Species, cov = vcv_brms_OU)),
                                        family =binomial(),
                                        data = data_neo_mal,
                                        data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                        chains = 4,
                                        cores = 24,
                                        iter = 10000, warmup = 5000,
                                        thin = 10,
                                        seed = 1,
                                        control=list(max_treedepth = 12, adapt_delta = 0.99),
                                        save_pars = save_pars(all = TRUE))

fit_reproduction_mal_adjusted_m3 
diagnostic(reference_model =fit_reproduction_mal_adjusted_m1, adjusted_model = fit_reproduction_mal_adjusted_m3)
#keep in the model (changes some slopes >25%)
#drop mean_GDP_5years  next

#model with mean_GDP_5years dropped
fit_reproduction_mal_adjusted_m4 <- brm(NMalignant | trials(NNecropsies)  ~  DietVertebrate + Littersize  + logGestationLength  + Placentation  + logNInstitutions + logHIndex + ResidualslogNViews + (1|gr(Species, cov = vcv_brms_OU)),
                                        family =binomial(),
                                        data = data_neo_mal,
                                        data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                        chains = 4,
                                        cores = 24,
                                        iter = 10000, warmup = 5000,
                                        thin = 10,
                                        seed = 1,
                                        control=list(max_treedepth = 12, adapt_delta = 0.99),
                                        save_pars = save_pars(all = TRUE))

fit_reproduction_mal_adjusted_m4 
diagnostic(reference_model =fit_reproduction_mal_adjusted_m1, adjusted_model = fit_reproduction_mal_adjusted_m4)
#changes the slope >25% keep

#model with the diet style dropped
fit_reproduction_mal_adjusted_m5 <- brm(NMalignant | trials(NNecropsies)  ~ Littersize  + logGestationLength  + Placentation  + logNInstitutions + logHIndex + mean_GDP_5years + ResidualslogNViews + (1|gr(Species, cov = vcv_brms_OU)),
                                        family =binomial(),
                                        data = data_neo_mal,
                                        data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                        chains = 4,
                                        cores = 24,
                                        iter = 10000, warmup = 5000,
                                        thin = 10,
                                        seed = 1,
                                        control=list(max_treedepth = 12, adapt_delta = 0.99),
                                        save_pars = save_pars(all = TRUE))

fit_reproduction_mal_adjusted_m5 
diagnostic(reference_model =fit_reproduction_mal_adjusted_m1, adjusted_model = fit_reproduction_mal_adjusted_m5)
#changes the slope >25% keep
#drop logNInstitutions  next

#model with logNInstitutions dropped
fit_reproduction_mal_adjusted_m6 <- brm(NMalignant | trials(NNecropsies)  ~ DietVertebrate + Littersize  + logGestationLength  + Placentation  + logHIndex + mean_GDP_5years + ResidualslogNViews + (1|gr(Species, cov = vcv_brms_OU)),
                                        family =binomial(),
                                        data = data_neo_mal,
                                        data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                        chains = 4,
                                        cores = 24,
                                        iter = 10000, warmup = 5000,
                                        thin = 10,
                                        seed = 1,
                                        control=list(max_treedepth = 12, adapt_delta = 0.99),
                                        save_pars = save_pars(all = TRUE))

fit_reproduction_mal_adjusted_m6 
diagnostic(reference_model =fit_reproduction_mal_adjusted_m1, adjusted_model = fit_reproduction_mal_adjusted_m6)
#changes the slope >25% keep
#drop ResidualslogNViews next


#model with ResidualslogNViews  dropped
fit_reproduction_mal_adjusted_m7 <- brm(NMalignant | trials(NNecropsies)  ~ DietVertebrate + Littersize  + logGestationLength  + Placentation  + logNInstitutions + logHIndex + mean_GDP_5years  + (1|gr(Species, cov = vcv_brms_OU)),
                                        family =binomial(),
                                        data = data_neo_mal,
                                        data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                        chains = 4,
                                        cores = 24,
                                        iter = 10000, warmup = 5000,
                                        thin = 10,
                                        seed = 1,
                                        control=list(max_treedepth = 12, adapt_delta = 0.99),
                                        save_pars = save_pars(all = TRUE))

fit_reproduction_mal_adjusted_m7 
diagnostic(reference_model =fit_reproduction_mal_adjusted_m1, adjusted_model = fit_reproduction_mal_adjusted_m7)
#changes the slope >25% keep

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Step 3 Final model ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#add back logMass to the model
fit_reproduction_mal_adjusted_m8 <- brm(NMalignant | trials(NNecropsies)  ~ logMass + DietVertebrate + Littersize  + logGestationLength  + Placentation  + logNInstitutions + logHIndex + mean_GDP_5years + ResidualslogNViews + (1|gr(Species, cov = vcv_brms_OU)),
                                        family =binomial(),
                                        data = data_neo_mal,
                                        data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                        chains = 4,
                                        cores = 24,
                                        iter = 10000, warmup = 5000,
                                        thin = 10,
                                        seed = 1,
                                        control=list(max_treedepth = 12, adapt_delta = 0.99),
                                        save_pars = save_pars(all = TRUE))

fit_reproduction_mal_adjusted_m8
diagnostic(reference_model =fit_reproduction_mal_adjusted_m1, adjusted_model = fit_reproduction_mal_adjusted_m8)
#changes the slope >25% keep

#add back logLactationLength to the model
fit_reproduction_mal_adjusted_m9 <- brm(NMalignant | trials(NNecropsies)  ~ logMass + DietVertebrate + Littersize + logLactationLength   + logGestationLength  + Placentation  + logNInstitutions + logHIndex + mean_GDP_5years + ResidualslogNViews + (1|gr(Species, cov = vcv_brms_OU)),
                                         family =binomial(),
                                         data = data_neo_mal,
                                         data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                         chains = 4,
                                         cores = 24,
                                         iter = 10000, warmup = 5000,
                                         thin = 10,
                                         seed = 1,
                                         control=list(max_treedepth = 12, adapt_delta = 0.99),
                                         save_pars = save_pars(all = TRUE))

fit_reproduction_mal_adjusted_m9 
diagnostic(reference_model =fit_reproduction_mal_adjusted_m8, adjusted_model = fit_reproduction_mal_adjusted_m9)
#changes slope > 25% keep

#add back logMaxLongevity to the model
fit_reproduction_mal_adjusted_m10 <- brm(NMalignant | trials(NNecropsies)  ~ logMass+ logMaxLongevity + DietVertebrate + Littersize + logLactationLength   + logGestationLength  + Placentation  + logNInstitutions + logHIndex + mean_GDP_5years + ResidualslogNViews + (1|gr(Species, cov = vcv_brms_OU)),
                                        family =binomial(),
                                        data = data_neo_mal,
                                        data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                        chains = 4,
                                        cores = 24,
                                        iter = 10000, warmup = 5000,
                                        thin = 10,
                                        seed = 1,
                                        control=list(max_treedepth = 12, adapt_delta = 0.99),
                                        save_pars = save_pars(all = TRUE))

fit_reproduction_mal_adjusted_m10 
diagnostic(reference_model =fit_reproduction_mal_adjusted_m9, adjusted_model = fit_reproduction_mal_adjusted_m10)
#logMaxLongevity changes slopes <25% -drop

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Calculate and plot the marginal effects ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#crude model
crude_model_mal <- brm(NMalignant | trials(NNecropsies)  ~ logGestationLength  + Placentation + (1|gr(Species, cov = vcv_brms_OU)),
                                    family =binomial(),
                                    data = data_neo_mal,
                                    data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                    chains = 4,
                                    cores = 24,
                                    threads = threading(2),
                                    iter = 10000, warmup = 5000,
                                    thin = 10,
                                    seed = 1,
                                    control=list(max_treedepth = 12, adapt_delta = 0.99),
                                    save_pars = save_pars(all = TRUE))
crude_model_mal

#create the conditions for the crude and adjusted models
cdt_crude_mal <-  list(NNecropsies = 100,
                       logLactationLength = median_logLactationLength,
                       logGestationLength  = median_logGestationLength,  
                       Littersize = median_Littersize,
                       DietVertebrate = 1)

cdt_adj_mal <- list( logNInstitutions = median_logNInstitutions,
                     logMass = median_logMass,
                     logMaxLongevity = median_logMaxLongevity,
                     ResidualslogNViews = median_ResidualslogNViews,
                     mean_GDP_5years= median_mean_GDP_5years,
                     logLactationLength = median_logLactationLength,
                     logGestationLength  = median_logGestationLength,  
                     logNExSitu   = median_logNExSitu,
                     NNecropsies = 100,
                     DietVertebrate = 1,
                     Littersize = median_Littersize)


conditional_effects(fit_reproduction_mal_adjusted_m11 , conditions = cdt_adj_mal, effect = "Placentation", resp = "Placentation")

marginal_effects_mal <- get_marginal_effects(crude_model = crude_model_mal , 
                                             adjusted_model = fit_reproduction_mal_adjusted_m9,
                                             resp = "NMalignant",
                                             effect = "Placentation",
                                             varext = "Placentation",
                                             NNecropsies = 100, 
                                             conditions_crude = cdt_crude_mal, 
                                             conditions_adjusted = cdt_adj_mal)

#plot the marginal effects for mutation rate
df_malignant <- marginal_effects_mal
p_reproduction_mal <- ggplot(data = df_malignant , aes(x = Placentation, y = estimate__, group = Type, color = Type, fill = Type))+
  geom_point(size = 4, position = position_dodge(width = 0.5))+
  geom_errorbar(mapping = aes(x = Placentation, ymin = lower__, ymax = upper__), alpha = 1, position = position_dodge(width = 0.5),  width = 0, size = 1.5)+
  theme_classic()+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=20,face="bold"))+
  theme(legend.position = c(0.2, 0.90))+
  theme(plot.title = element_text(size = 18, face = "bold"))+
  theme(legend.text=element_text(size=12),
        legend.title = element_text(size = 12))+
  scale_color_manual(values = c("Crude" = "#E69F00", "Adjusted" = "#009E73")) +
  scale_fill_manual(values = c("Crude" = "#E69F00", "Adjusted" = "#009E73")) +
  #theme(axis.title.y = element_blank())+  # Remove axis titles
  xlab("Placentation type")+
  ylab("Malignancy prevalence")+
  labs(color = "Model type:", fill = "Model type:")+
  ylim(0,0.3)+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
#xlim(-3.5,3.5)
p_reproduction_mal

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# NMortality ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Step 1: Univariate models ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#fit univariate model for logMass
fit_reproduction_mor_logMass <- brm(NMortality | trials(NNecropsies)  ~ logMass + (1|gr(Species, cov = vcv_brms_OU)),
                                    family =binomial(),
                                    data = data_mortality,
                                    data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                                    chains = 4,
                                    cores = 24,
                                    threads = threading(2),
                                    iter = 10000, warmup = 5000,
                                    thin = 10,
                                    seed = 1,
                                    control=list(max_treedepth = 12, adapt_delta = 0.99),
                                    save_pars = save_pars(all = TRUE))
fit_reproduction_mor_logMass
p_direction(fit_reproduction_mor_logMass)
#p_direction <0.90 -> keep in the model

#fit univariate model for logMaxLongevity
fit_reproduction_mor_logMaxLongevity <- brm(NMortality | trials(NNecropsies)  ~ logMaxLongevity + (1|gr(Species, cov = vcv_brms_OU)),
                                            family =binomial(),
                                            data = data_mortality,
                                            data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                                            chains = 4,
                                            cores = 24,
                                            threads = threading(2),
                                            iter = 10000, warmup = 5000,
                                            thin = 10,
                                            seed = 1,
                                            control=list(max_treedepth = 12, adapt_delta = 0.99),
                                            save_pars = save_pars(all = TRUE))
fit_reproduction_mor_logMaxLongevity
p_direction(fit_reproduction_mor_logMaxLongevity)
#p_direction <0.90 -> exclude from the model -

#fit univariate model for logLactationLength
fit_reproduction_mor_logLactationLength <- brm(NMortality | trials(NNecropsies)  ~ logLactationLength + (1|gr(Species, cov = vcv_brms_OU)),
                                               family =binomial(),
                                               data = data_mortality,
                                               data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                                               chains = 4,
                                               cores = 24,
                                               threads = threading(2),
                                               iter = 10000, warmup = 5000,
                                               thin = 10,
                                               seed = 1,
                                               control=list(max_treedepth = 12, adapt_delta = 0.99),
                                               save_pars = save_pars(all = TRUE))
fit_reproduction_mor_logLactationLength
p_direction(fit_reproduction_mor_logLactationLength)
#p_direction <0.90 -> exclude from the model

#fit univariate model for logGestationLength
fit_reproduction_mor_logGestationLength <- brm(NMortality | trials(NNecropsies)  ~ logGestationLength  + (1|gr(Species, cov = vcv_brms_OU)),
                                               family =binomial(),
                                               data = data_mortality,
                                               data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                                               chains = 4,
                                               cores = 24,
                                               threads = threading(2),
                                               iter = 10000, warmup = 5000,
                                               thin = 10,
                                               seed = 1,
                                               control=list(max_treedepth = 12, adapt_delta = 0.99),
                                               save_pars = save_pars(all = TRUE))
fit_reproduction_mor_logGestationLength
p_direction(fit_reproduction_mor_logGestationLength)
#p_direction >0.90 -> retain for multivariate model -

#fit univariate model for placentation type
fit_reproduction_mor_Placentation <- brm(NMortality | trials(NNecropsies)  ~ Placentation  + (1|gr(Species, cov = vcv_brms_OU)),
                                         family =binomial(),
                                         data = data_mortality,
                                         data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                                         chains = 4,
                                         cores = 24,
                                         threads = threading(2),
                                         iter = 10000, warmup = 5000,
                                         thin = 10,
                                         seed = 1,
                                         control=list(max_treedepth = 12, adapt_delta = 0.99),
                                         save_pars = save_pars(all = TRUE))
fit_reproduction_mor_Placentation
p_direction(fit_reproduction_mor_Placentation)
#p_direction >0.90 -> retain for multivariate model -

#fit univariate model for DietVertebrate 
fit_reproduction_mor_DietVertebrate <-  brm(NMortality | trials(NNecropsies)  ~ DietVertebrate   + (1|gr(Species, cov = vcv_brms_OU)),
                                            family =binomial(),
                                            data = data_mortality,
                                            data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                                            chains = 4,
                                            cores = 24,
                                            threads = threading(2),
                                            iter = 10000, warmup = 5000,
                                            thin = 10,
                                            seed = 1,
                                            control=list(max_treedepth = 12, adapt_delta = 0.99),
                                            save_pars = save_pars(all = TRUE))
fit_reproduction_mor_DietVertebrate
p_direction(fit_reproduction_mor_DietVertebrate)
#p_direction >0.90 ->  retain for multivariate model -

#fit univariate model for LitterSize 
fit_reproduction_mor_Littersize <- brm(NMortality | trials(NNecropsies)  ~ Littersize   + (1|gr(Species, cov = vcv_brms_OU)),
                                       family =binomial(),
                                       data = data_mortality,
                                       data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                                       chains = 4,
                                       cores = 24,
                                       threads = threading(2),
                                       iter = 10000, warmup = 5000,
                                       thin = 10,
                                       seed = 1,
                                       control=list(max_treedepth = 12, adapt_delta = 0.99),
                                       save_pars = save_pars(all = TRUE))
fit_reproduction_mor_Littersize
p_direction(fit_reproduction_mor_Littersize)
#p_direction >0.90 ->  retain for multivariate model -

#fit univariate model for logHIndex
fit_reproduction_mor_logHIndex <- brm(NMortality | trials(NNecropsies)  ~  logHIndex   + (1|gr(Species, cov = vcv_brms_OU)),
                                      family =binomial(),
                                      data = data_mortality,
                                      data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                                      chains = 4,
                                      cores = 24,
                                      threads = threading(2),
                                      iter = 10000, warmup = 5000,
                                      thin = 10,
                                      seed = 1,
                                      control=list(max_treedepth = 12, adapt_delta = 0.99),
                                      save_pars = save_pars(all = TRUE))
fit_reproduction_mor_logHIndex
p_direction(fit_reproduction_mor_logHIndex)
#p_direction >0.90 ->  retain for multivariate model -


#fit univariate model for mean_GDP_5years
fit_reproduction_mor_mean_GDP_5years <- brm(NMortality | trials(NNecropsies)  ~  mean_GDP_5years   + (1|gr(Species, cov = vcv_brms_OU)),
                                            family =binomial(),
                                            data = data_mortality,
                                            data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                                            chains = 4,
                                            cores = 24,
                                            threads = threading(2),
                                            iter = 10000, warmup = 5000,
                                            thin = 10,
                                            seed = 1,
                                            control=list(max_treedepth = 12, adapt_delta = 0.99),
                                            save_pars = save_pars(all = TRUE))
fit_reproduction_mor_mean_GDP_5years
p_direction(fit_reproduction_mor_mean_GDP_5years)
#p_direction <0.90 ->  remove from model -

#fit univariate model for logNInstitutions
fit_reproduction_mal_logNInstitutions <- brm(NMortality | trials(NNecropsies)  ~  logNInstitutions   + (1|gr(Species, cov = vcv_brms_OU)),
                                             family =binomial(),
                                             data = data_mortality,
                                             data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                                             chains = 4,
                                             cores = 24,
                                             threads = threading(2),
                                             iter = 10000, warmup = 5000,
                                             thin = 10,
                                             seed = 1,
                                             control=list(max_treedepth = 12, adapt_delta = 0.99),
                                             save_pars = save_pars(all = TRUE))
fit_reproduction_mal_logNInstitutions
p_direction(fit_reproduction_mal_logNInstitutions)
#p_direction <0.90 ->  remove drom the model -

#fit univariate model for ResidualslogNViews
fit_reproduction_mor_ResidualslogNViews  <- brm(NMortality | trials(NNecropsies)  ~  ResidualslogNViews   + (1|gr(Species, cov = vcv_brms_OU)),
                                                family =binomial(),
                                                data = data_mortality,
                                                data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                                                chains = 4,
                                                cores = 24,
                                                threads = threading(2),
                                                iter = 10000, warmup = 5000,
                                                thin = 10,
                                                seed = 1,
                                                control=list(max_treedepth = 12, adapt_delta = 0.99),
                                                save_pars = save_pars(all = TRUE))
fit_reproduction_mor_ResidualslogNViews
p_direction(fit_reproduction_mor_ResidualslogNViews)
#p_direction >0.90 ->  retain for multivariate model -

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Step 2 : Drop variables part of the model selection process ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#fit a model with all the variables selected in step 1
fit_reproduction_mor_adjusted_m1 <- brm(NMortality | trials(NNecropsies)  ~  logMass + Placentation + logGestationLength + DietVertebrate  + Littersize + logHIndex + ResidualslogNViews + (1|gr(Species, cov = vcv_brms_OU)),
                                        family =binomial(),
                                        data = data_mortality,
                                        data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                                        chains = 4,
                                        cores = 24,
                                        threads = threading(2),
                                        iter = 10000, warmup = 5000,
                                        thin = 10,
                                        seed = 1,
                                        control=list(max_treedepth = 12, adapt_delta = 0.99),
                                        save_pars = save_pars(all = TRUE))

fit_reproduction_mor_adjusted_m1 
summary(fit_reproduction_mor_adjusted_m1 , prob = 0.90)
p_direction(fit_reproduction_mor_adjusted_m1) 
#Littersize   has the lowest p-direction -> drop first

#drop Littersize  first
fit_reproduction_mor_adjusted_m2 <- brm(NMortality | trials(NNecropsies)  ~  logMass + Placentation + logGestationLength + DietVertebrate  + logHIndex + ResidualslogNViews + (1|gr(Species, cov = vcv_brms_OU)),
                                        family =binomial(),
                                        data = data_mortality,
                                        data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                                        chains = 4,
                                        cores = 24,
                                        threads = threading(2),
                                        iter = 10000, warmup = 5000,
                                        thin = 10,
                                        seed = 1,
                                        control=list(max_treedepth = 12, adapt_delta = 0.99),
                                        save_pars = save_pars(all = TRUE))

fit_reproduction_mor_adjusted_m2 
diagnostic(reference_model =fit_reproduction_mor_adjusted_m1, adjusted_model = fit_reproduction_mor_adjusted_m2)
p_direction(fit_reproduction_mor_adjusted_m2) 
# Littersize change the slopes <25% -> drop from the model
# drop logMass next

#drop logMass
fit_reproduction_mor_adjusted_m3 <- brm(NMortality | trials(NNecropsies)  ~   Placentation + logGestationLength + DietVertebrate  + logHIndex + ResidualslogNViews + (1|gr(Species, cov = vcv_brms_OU)),
                                        family =binomial(),
                                        data = data_mortality,
                                        data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                                        chains = 4,
                                        cores = 24,
                                        threads = threading(2),
                                        iter = 10000, warmup = 5000,
                                        thin = 10,
                                        seed = 1,
                                        control=list(max_treedepth = 12, adapt_delta = 0.99),
                                        save_pars = save_pars(all = TRUE))

fit_reproduction_mor_adjusted_m3
diagnostic(reference_model =fit_reproduction_mor_adjusted_m2, adjusted_model = fit_reproduction_mor_adjusted_m3)
p_direction(fit_reproduction_mor_adjusted_m3) 
# logMass changes the slopes <25% (PlacentationHemochorial) -> keep in the model
#drop diet next


#model with placentation droped
fit_reproduction_mor_adjusted_m4 <- brm(NMortality | trials(NNecropsies)  ~  logMass + Placentation + logGestationLength  + logHIndex + ResidualslogNViews + (1|gr(Species, cov = vcv_brms_OU)),
                                        family =binomial(),
                                        data = data_mortality,
                                        data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                                        chains = 4,
                                        cores = 24,
                                        threads = threading(2),
                                        iter = 10000, warmup = 5000,
                                        thin = 10,
                                        seed = 1,
                                        control=list(max_treedepth = 12, adapt_delta = 0.99),
                                        save_pars = save_pars(all = TRUE))

fit_reproduction_mor_adjusted_m4
diagnostic(reference_model =fit_reproduction_mor_adjusted_m2, adjusted_model = fit_reproduction_mor_adjusted_m4)
#changes the slope >25% -> keep in the model
#drop logHIndex next  next

#drop logHIndex 
fit_reproduction_mor_adjusted_m5 <- brm(NMortality | trials(NNecropsies)  ~  logMass + Placentation + logGestationLength + DietVertebrate  + ResidualslogNViews + (1|gr(Species, cov = vcv_brms_OU)),
                                        family =binomial(),
                                        data = data_mortality,
                                        data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                                        chains = 4,
                                        cores = 24,
                                        threads = threading(2),
                                        iter = 10000, warmup = 5000,
                                        thin = 10,
                                        seed = 1,
                                        control=list(max_treedepth = 12, adapt_delta = 0.99),
                                        save_pars = save_pars(all = TRUE))

fit_reproduction_mor_adjusted_m5
diagnostic(reference_model =fit_reproduction_mor_adjusted_m2, adjusted_model = fit_reproduction_mor_adjusted_m5)
#logHIndex  changes the intercept and logHIndex >25% -> keep
#drop logHIndex next

#model with ResidualslogNViews droped
fit_reproduction_mor_adjusted_m6 <- brm(NMortality | trials(NNecropsies)  ~  logMass + Placentation + logGestationLength + DietVertebrate  + logHIndex + (1|gr(Species, cov = vcv_brms_OU)),
                                        family =binomial(),
                                        data = data_mortality,
                                        data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                                        chains = 4,
                                        cores = 24,
                                        threads = threading(2),
                                        iter = 10000, warmup = 5000,
                                        thin = 10,
                                        seed = 1,
                                        control=list(max_treedepth = 12, adapt_delta = 0.99),
                                        save_pars = save_pars(all = TRUE))

fit_reproduction_mor_adjusted_m6
diagnostic(reference_model =fit_reproduction_mor_adjusted_m2, adjusted_model = fit_reproduction_mor_adjusted_m6)
p_direction(fit_reproduction_mor_adjusted_m6) 
#changes slopes >25% keep in the model

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Step 3 Final model ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#add the logLactationLegnt back to the model
fit_reproduction_mor_adjusted_m7 <- brm(NMortality | trials(NNecropsies)  ~  logMass + Placentation + logGestationLength + logLactationLength + DietVertebrate  + logHIndex + ResidualslogNViews + (1|gr(Species, cov = vcv_brms_OU)),
                                        family =binomial(),
                                        data = data_mortality,
                                        data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                                        chains = 4,
                                        cores = 24,
                                        threads = threading(2),
                                        iter = 10000, warmup = 5000,
                                        thin = 10,
                                        seed = 1,
                                        control=list(max_treedepth = 12, adapt_delta = 0.99),
                                        save_pars = save_pars(all = TRUE))

fit_reproduction_mor_adjusted_m7
diagnostic(reference_model =fit_reproduction_mor_adjusted_m2, adjusted_model = fit_reproduction_mor_adjusted_m7)
p_direction(fit_reproduction_mor_adjusted_m7) 
#changes slopes >25% -> keep

#add the logNInstutions back to the model
fit_reproduction_mor_adjusted_m8 <- brm(NMortality | trials(NNecropsies)  ~  logMass + Placentation + logGestationLength + logLactationLength + DietVertebrate  +	logNInstitutions +  logHIndex + ResidualslogNViews + (1|gr(Species, cov = vcv_brms_OU)),
                                        family =binomial(),
                                        data = data_mortality,
                                        data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                                        chains = 4,
                                        cores = 24,
                                        threads = threading(2),
                                        iter = 10000, warmup = 5000,
                                        thin = 10,
                                        seed = 1,
                                        control=list(max_treedepth = 12, adapt_delta = 0.99),
                                        save_pars = save_pars(all = TRUE))

fit_reproduction_mor_adjusted_m8
diagnostic(reference_model =fit_reproduction_mor_adjusted_m7, adjusted_model = fit_reproduction_mor_adjusted_m8)
p_direction(fit_reproduction_mor_adjusted_m8) 
#keep in the model

#add back logMaxLongevity
fit_reproduction_mor_adjusted_m9 <- brm(NMortality | trials(NNecropsies)  ~  logMass + logMaxLongevity + Placentation + logGestationLength + logLactationLength + DietVertebrate  +	logNInstitutions +  logHIndex + ResidualslogNViews + (1|gr(Species, cov = vcv_brms_OU)),
                                         family =binomial(),
                                         data = data_mortality,
                                         data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                                         chains = 4,
                                         cores = 24,
                                         threads = threading(2),
                                         iter = 10000, warmup = 5000,
                                         thin = 10,
                                         seed = 1,
                                         control=list(max_treedepth = 12, adapt_delta = 0.99),
                                         save_pars = save_pars(all = TRUE))

fit_reproduction_mor_adjusted_m9
diagnostic(reference_model =fit_reproduction_mor_adjusted_m8, adjusted_model = fit_reproduction_mor_adjusted_m9)
p_direction(fit_reproduction_mor_adjusted_m9) 
#drop logMaxLongevity


#add back mean_GDP_5years
fit_reproduction_mor_adjusted_m10 <- brm(NMortality | trials(NNecropsies)  ~  logMass  + Placentation + logGestationLength + logLactationLength + DietVertebrate +	logNInstitutions +  logHIndex + ResidualslogNViews + mean_GDP_5years + (1|gr(Species, cov = vcv_brms_OU)),
                                         family =binomial(),
                                         data = data_mortality,
                                         data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                                         chains = 4,
                                         cores = 24,
                                         threads = threading(2),
                                         iter = 10000, warmup = 5000,
                                         thin = 10,
                                         seed = 1,
                                         control=list(max_treedepth = 12, adapt_delta = 0.99),
                                         save_pars = save_pars(all = TRUE))

fit_reproduction_mor_adjusted_m10
diagnostic(reference_model =fit_reproduction_mor_adjusted_m8, adjusted_model = fit_reproduction_mor_adjusted_m10)
p_direction(fit_reproduction_mor_adjusted_m10) 
#keep this one as a final model

