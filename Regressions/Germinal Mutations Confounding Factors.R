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
phylo_all <- read.tree("phylogeny germinal mutations.nwk") 

  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import Neoplasia and malignancy data ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("C:/Users/antoi/My Drive/Postdoc/Zoo data/Biases/Scripts/Estimation of confounding effects/for github")

data_mortality <- read_excel("Raw datasets.xlsx", sheet = "Germinal_mutations")
rownames(data_mortality) <- data_mortality$Species
  
#calculate the median of confounding variables to use in the calculation of marginal effects
median_logHIndex <- median(data_mortality$logHIndex)
median_logNInstitutions <- median(data_mortality$logNInstitutions)
median_logNExSitu <- median(data_mortality$logNExSitu)
median_ResidualslogNViews <- median(data_mortality$ResidualslogNViews)
median_mean_GDP_5years <- median(data_mortality$mean_GDP_5years)
median_logNExSitu <- median(data_mortality$logNExSitu)
median_logMass <- median(data_mortality$logMass)
median_logMaxLongevity <- median(data_mortality$logMaxLongevity)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Prepare var covar matrixes for brms models ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
# Step 1: Univariate models ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Alpha = 0.1

#univariate model with log mass
fit_mortality_mass <- brm( NMortality | trials(NNecropsies)  ~  logMass + (1|gr(Species, cov = vcv_brms_OU)),
                               family =binomial(),
                               data = data_mortality,
                               data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                               chains = 4,
                               cores = 24,
                               iter = 10000, warmup = 5000,
                               thin = 10,
                               seed = 1,
                               control=list(max_treedepth = 12, adapt_delta = 0.99),
                               save_pars = save_pars(all = TRUE))

summary(fit_mortality_mass, prob = 0.90)
p_direction(fit_mortality_mass) 
#p_direction >0.90 -> retain for multivariate model -


#univariate model with log longevity
fit_mortality_longevity <- brm( NMortality | trials(NNecropsies)  ~  logMaxLongevity + (1|gr(Species, cov = vcv_brms_OU)),
                           family =binomial(),
                           data = data_mortality,
                           data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                           chains = 4,
                           cores = 24,
                           iter = 10000, warmup = 5000,
                           thin = 10,
                           seed = 1,
                           control=list(max_treedepth = 12, adapt_delta = 0.99),
                           save_pars = save_pars(all = TRUE))

summary(fit_mortality_longevity, prob = 0.90)
p_direction(fit_mortality_longevity) 
#p_direction >0.90 -> retain for multivariate model


#univariate model with trophic level
fit_mortality_trophic <- brm( NMortality | trials(NNecropsies)  ~  Tophic.Level + (1|gr(Species, cov = vcv_brms_OU)),
                                family =binomial(),
                                data = data_mortality,
                                data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                                chains = 4,
                                cores = 24,
                                iter = 10000, warmup = 5000,
                                thin = 10,
                                seed = 1,
                                control=list(max_treedepth = 12, adapt_delta = 0.99),
                                save_pars = save_pars(all = TRUE))

summary(fit_mortality_trophic, prob = 0.90)
p_direction(fit_mortality_trophic) 
#p_direction <0.90 -> exclude from multivariate model
 
#univariate model with mutation rate
fit_mortality_mutation <- brm( NMortality | trials(NNecropsies)  ~  logMutationRate + (1|gr(Species, cov = vcv_brms_OU)),
                                      family =binomial(),
                                      data = data_mortality,
                                      data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                                      chains = 4,
                                      cores = 24,
                                      iter = 10000, warmup = 5000,
                                      thin = 10,
                                      seed = 1,
                                      control=list(max_treedepth = 12, adapt_delta = 0.99),
                                      save_pars = save_pars(all = TRUE))

summary(fit_mortality_mutation, prob = 0.90)
p_direction(fit_mortality_mutation) 
#p_direction >0.90 -> retain for multivariate model

#model with logHIndex
fit_mortality_logHIndex <- brm( NMortality | trials(NNecropsies)  ~ logHIndex + (1|gr(Species, cov = vcv_brms_OU)),
                      family =binomial(),
                      data = data_mortality,
                      data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                      chains = 4,
                      cores = 24,
                      iter = 10000, warmup = 5000,
                      thin = 10,
                      seed = 1,
                      control=list(max_treedepth = 12, adapt_delta = 0.99),
                      save_pars = save_pars(all = TRUE))

summary(fit_mortality_logHIndex, prob = 0.90)
p_direction(fit_mortality_logHIndex) 
#p_direction >0.90 -> retain for multivariate model

#model with logNInstitutions
fit_mortality_logNInstitutions <- brm( NMortality | trials(NNecropsies)  ~ logNInstitutions + (1|gr(Species, cov = vcv_brms_OU)),
                                family =binomial(),
                                data = data_mortality,
                                data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                                chains = 4,
                                cores = 24,
                                iter = 10000, warmup = 5000,
                                thin = 10,
                                seed = 1,
                                control=list(max_treedepth = 12, adapt_delta = 0.99),
                                save_pars = save_pars(all = TRUE))

summary(fit_mortality_logNInstitutions, prob = 0.90)
p_direction(fit_mortality_logNInstitutions) 
#p_direction <0.90 -> exclude from multivariate model

#model with mean_GDP_5years
fit_mortality_mean_GDP_5years <- brm( NMortality | trials(NNecropsies)  ~ mean_GDP_5years + (1|gr(Species, cov = vcv_brms_OU)),
                                 family =binomial(),
                                 data = data_mortality,
                                 data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                                 chains = 4,
                                 cores = 24,
                                 iter = 10000, warmup = 5000,
                                 thin = 10,
                                 seed = 1,
                                 control=list(max_treedepth = 12, adapt_delta = 0.99),
                                 save_pars = save_pars(all = TRUE))

summary(fit_mortality_mean_GDP_5years, prob = 0.90)
p_direction(fit_mortality_mean_GDP_5years) 
#p_direction >0.90 -> retain for multivariate model

#model with ResidualslogNViews
fit_mortality_ResidualslogNViews <- brm( NMortality | trials(NNecropsies)  ~ ResidualslogNViews + (1|gr(Species, cov = vcv_brms_OU)),
                                      family =binomial(),
                                      data = data_mortality,
                                      data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                                      chains = 4,
                                      cores = 24,
                                      iter = 10000, warmup = 5000,
                                      thin = 10,
                                      seed = 1,
                                      control=list(max_treedepth = 12, adapt_delta = 0.99),
                                      save_pars = save_pars(all = TRUE))

summary(fit_mortality_ResidualslogNViews, prob = 0.90)
p_direction(fit_mortality_ResidualslogNViews) 
#p_direction <0.90 -> exclude from multivariate model


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Step 2: Drop part of the model selection ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#fit a model with all variables retained in step 1 
fit_adjusted_mortality_m1 <- brm(NMortality | trials(NNecropsies)  ~ logMass + logMaxLongevity + logMutationRate + logHIndex + mean_GDP_5years + (1|gr(Species, cov = vcv_brms_OU)),
                                 family =binomial(),
                                 data = data_mortality,
                                 data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                                 chains = 4,
                                 cores = 24,
                                 iter = 10000, warmup = 5000,
                                 thin = 10,
                                 seed = 1,
                                 control=list(max_treedepth = 12, adapt_delta = 0.99),
                                 save_pars = save_pars(all = TRUE))

fit_adjusted_mortality_m1
p_direction(fit_adjusted_mortality_m1) 
#logMass has the lowest p-direction (outside logMutationRate)

#fit a model with logMass dropped
fit_adjusted_mortality_m2 <- brm(NMortality | trials(NNecropsies)  ~ logMaxLongevity + logMutationRate + logHIndex + mean_GDP_5years + (1|gr(Species, cov = vcv_brms_OU)),
                                 family =binomial(),
                                 data = data_mortality,
                                 data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                                 chains = 4,
                                 cores = 24,
                                 iter = 10000, warmup = 5000,
                                 thin = 10,
                                 seed = 1,
                                 control=list(max_treedepth = 12, adapt_delta = 0.99),
                                 save_pars = save_pars(all = TRUE))

fit_adjusted_mortality_m2

diagnostic(reference_model =fit_adjusted_mortality_m1, adjusted_model = fit_adjusted_mortality_m2)
#some slopes changed >25% (including logMutationRate) ->  keep logMass in the model
#mean_GDP_5years has the lowest p-direction in model m1 -> drop next

#fit a model with mean_GDP_5years dropped
fit_adjusted_mortality_m3 <- brm(NMortality | trials(NNecropsies)  ~  logMass + logMaxLongevity + logMutationRate + logHIndex  + (1|gr(Species, cov = vcv_brms_OU)),
                                 family =binomial(),
                                 data = data_mortality,
                                 data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                                 chains = 4,
                                 cores = 24,
                                 iter = 10000, warmup = 5000,
                                 thin = 10,
                                 seed = 1,
                                 control=list(max_treedepth = 12, adapt_delta = 0.99),
                                 save_pars = save_pars(all = TRUE))

fit_adjusted_mortality_m3
diagnostic(reference_model =fit_adjusted_mortality_m1, adjusted_model = fit_adjusted_mortality_m3)
#dropping mean_GDP_5years changes the slope of logMutation rate >25% -> retain in the model
#keep model m2 and then next drop logMaxLongevity to test its effect

#fit a model with logMaxLongevity dropped
fit_adjusted_mortality_m4 <- brm(NMortality | trials(NNecropsies)  ~ logMass + logMutationRate + logHIndex + mean_GDP_5years + (1|gr(Species, cov = vcv_brms_OU)),
                                 family =binomial(),
                                 data = data_mortality,
                                 data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                                 chains = 4,
                                 cores = 24,
                                 iter = 10000, warmup = 5000,
                                 thin = 10,
                                 seed = 1,
                                 control=list(max_treedepth = 12, adapt_delta = 0.99),
                                 save_pars = save_pars(all = TRUE))

fit_adjusted_mortality_m4
diagnostic(reference_model =fit_adjusted_mortality_m2, adjusted_model = fit_adjusted_mortality_m4)
#dropping logMaxLongevity changes the slope of logMutation rate >25% -> retain in the model
#keep model m2 and then next drop logHindex to test its effect

#fit a model with logHindex dropped
fit_adjusted_mortality_m5 <- brm(NMortality | trials(NNecropsies)  ~ logMass + logMaxLongevity + logMutationRate + mean_GDP_5years + (1|gr(Species, cov = vcv_brms_OU)),
                                 family =binomial(),
                                 data = data_mortality,
                                 data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                                 chains = 4,
                                 cores = 24,
                                 iter = 10000, warmup = 5000,
                                 thin = 10,
                                 seed = 1,
                                 control=list(max_treedepth = 12, adapt_delta = 0.99),
                                 save_pars = save_pars(all = TRUE))

fit_adjusted_mortality_m5
diagnostic(reference_model =fit_adjusted_mortality_m2, adjusted_model = fit_adjusted_mortality_m5)
#dropping logHindex changes the slope of logMutationRate rate >25% -> retain in the model
#model m2 selected for the next step

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 3: Reintroducing variables dropped in Step 1 in the model ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#refit the model selected in step 2
#logMaxLongevity + logMutationRate + logHIndex + mean_GDP_5years 
fit_adjusted_mortality_m6 <- brm(NMortality | trials(NNecropsies)  ~  logMass + logMaxLongevity + logMutationRate + logHIndex + mean_GDP_5years + (1|gr(Species, cov = vcv_brms_OU)),
                                 family =binomial(),
                                 data = data_mortality,
                                 data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                                 chains = 4,
                                 cores = 24,
                                 iter = 10000, warmup = 5000,
                                 thin = 10,
                                 seed = 1,
                                 control=list(max_treedepth = 12, adapt_delta = 0.99),
                                 save_pars = save_pars(all = TRUE))

fit_adjusted_mortality_m6 
p_direction(fit_adjusted_mortality_m6) 

#reintroducing logNInstitutions
fit_adjusted_mortality_m7 <- brm(NMortality | trials(NNecropsies)  ~  logMass + logMaxLongevity + logMutationRate + logHIndex + mean_GDP_5years + logNInstitutions + (1|gr(Species, cov = vcv_brms_OU)),
                                 family =binomial(),
                                 data = data_mortality,
                                 data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                                 chains = 4,
                                 cores = 24,
                                 iter = 10000, warmup = 5000,
                                 thin = 10,
                                 seed = 1,
                                 control=list(max_treedepth = 12, adapt_delta = 0.99),
                                 save_pars = save_pars(all = TRUE))

fit_adjusted_mortality_m7 
diagnostic(reference_model =fit_adjusted_mortality_m6, adjusted_model = fit_adjusted_mortality_m7)
#logNInstitutions has an effect on the slope of logMass <25% -> drop from the model


#Add ResidualslogNViews to the model selected selected in step 2
fit_adjusted_mortality_m8 <- brm(NMortality | trials(NNecropsies)  ~  logMass + logMaxLongevity + logMutationRate + logHIndex + mean_GDP_5years  + ResidualslogNViews + (1|gr(Species, cov = vcv_brms_OU)),
                                 family =binomial(),
                                 data = data_mortality,
                                 data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                                 chains = 4,
                                 cores = 24,
                                 iter = 10000, warmup = 5000,
                                 thin = 10,
                                 seed = 1,
                                 control=list(max_treedepth = 12, adapt_delta = 0.99),
                                 save_pars = save_pars(all = TRUE))

fit_adjusted_mortality_m8 

diagnostic(reference_model =fit_adjusted_mortality_m6, adjusted_model = fit_adjusted_mortality_m8)
#ResidualslogNViews has an effect on the slope of logMass    >25% -> retain in the model

fit_adjusted_mortality_final <- fit_adjusted_mortality_m8 
pp_check(fit_adjusted_mortality_final, ndraws = 100)

