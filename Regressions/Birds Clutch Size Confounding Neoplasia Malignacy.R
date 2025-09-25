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

import_dir <- "C:/Users/antoi/My Drive/Postdoc/Zoo data/Biases/Data"
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import Neoplasia and malignancy data ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("C:/Users/antoi/My Drive/Postdoc/Zoo data/Biases/Scripts/Estimation of confounding effects/for github")
data_neo_mal <- read_excel("Raw datasets.xlsx", sheet = "Clutch_Size")
rownames(data_neo_mal) <- data_neo_mal$Species

#calculate the median of confounding variables to use in the calculation of marginal effects
median_logHindex <- median(data_neo_mal$logHIndex)
median_logNInstitutions <- median(data_neo_mal$logNInstitutions)
median_logNExSitu <- median(data_neo_mal$logNExSitu)
median_ResidualslogNViews <- median(data_neo_mal$ResidualslogNViews)
median_logNExSitu <- median(data_neo_mal$logNExSitu)
median_logMass <- median(data_neo_mal$logMass)
median_logMaxLongevity <- median(data_neo_mal$logMaxLongevity)
median_mean_GDP_5years <- median(data_neo_mal$mean_GDP_5years)
median_IncubationLength <- median(data_neo_mal$IncubationLength)

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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Neoplasia  ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Step 1 univariate models  ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#fit univariate model with logMass
fit_reproduction_neo_logMass <- brm(NNeoplasia | trials(NNecropsies)  ~ logMass + (1|gr(Species, cov = vcv_brms_OU)),
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
summary(fit_reproduction_neo_logMass, prob = 0.90)
p_direction(fit_reproduction_neo_logMass) 
#p_direction <0.90 -> exclude from the model - 

#fit univariate model with logMaxLongevity
fit_reproduction_neo_logMaxLongevity <- brm(NNeoplasia | trials(NNecropsies)  ~ logMaxLongevity + (1|gr(Species, cov = vcv_brms_OU)),
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
summary(fit_reproduction_neo_logMaxLongevity, prob = 0.90)
p_direction(fit_reproduction_neo_logMaxLongevity) 
#p_direction <0.90 -> exclude from the model -

#fit univariate model with Incubation length
fit_reproduction_neo_IncubationLength <- brm(NNeoplasia | trials(NNecropsies)  ~ IncubationLength + (1|gr(Species, cov = vcv_brms_OU)),
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
summary(fit_reproduction_neo_IncubationLength, prob = 0.90)
p_direction(fit_reproduction_neo_IncubationLength) 
#p_direction <0.90 -> exclude from the model - 

#fit univariate model with ClutchSize
fit_reproduction_neo_ClutchSize <- brm(NNeoplasia | trials(NNecropsies)  ~ ClutchSize + (1|gr(Species, cov = vcv_brms_OU)),
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
summary(fit_reproduction_neo_ClutchSize, prob = 0.90)
p_direction(fit_reproduction_neo_ClutchSize) 
#p_direction >0.90 -> retain for multivariate model -

#fit univariate model with logHIndex
fit_reproduction_neo_logHIndex <- brm(NNeoplasia | trials(NNecropsies)  ~ logHIndex + (1|gr(Species, cov = vcv_brms_OU)),
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
summary(fit_reproduction_neo_logHIndex, prob = 0.90)
p_direction(fit_reproduction_neo_logHIndex) 
#p_direction <90 -> exclude from multivariate model -

#fit univariate model with logNInstitutions
fit_reproduction_neo_logNInstitutions <- brm(NNeoplasia | trials(NNecropsies)  ~ logNInstitutions + (1|gr(Species, cov = vcv_brms_OU)),
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
summary(fit_reproduction_neo_logNInstitutions, prob = 0.90)
p_direction(fit_reproduction_neo_logNInstitutions) 
#p_direction >0.90 -> retain for multivariate model -

#fit univariate model with mean_GDP_5years
fit_reproduction_neo_mean_GDP_5years <- brm(NNeoplasia | trials(NNecropsies)  ~ mean_GDP_5years + (1|gr(Species, cov = vcv_brms_OU)),
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
summary(fit_reproduction_neo_mean_GDP_5years, prob = 0.90)
p_direction(fit_reproduction_neo_mean_GDP_5years) 
#p_direction <0.90 -> exclude from the model -


#fit univariate model with ResidualslogNViews
fit_reproduction_neo_ResidualslogNViews <- brm(NNeoplasia | trials(NNecropsies)  ~ ResidualslogNViews + (1|gr(Species, cov = vcv_brms_OU)),
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
summary(fit_reproduction_neo_ResidualslogNViews, prob = 0.90)
p_direction(fit_reproduction_neo_ResidualslogNViews) 
#p_direction <0.90 -> exclude from the model -

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Step 2 drop variables from the model  ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#fit a model with all variables retained in Step 1
fit_reproduction_neo_m1 <- brm(NNeoplasia | trials(NNecropsies)  ~ ClutchSize  + logNInstitutions + (1|gr(Species, cov = vcv_brms_OU)),
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
fit_reproduction_neo_m1
p_direction(fit_reproduction_neo_m1) 
#logNInstitutions has the lowest p-direction - drop it first

#drop logNInstitutions from the model
fit_reproduction_neo_m2 <- brm(NNeoplasia | trials(NNecropsies)  ~ ClutchSize + (1|gr(Species, cov = vcv_brms_OU)),
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
fit_reproduction_neo_m2

diagnostic(reference_model =fit_reproduction_neo_m1, adjusted_model = fit_reproduction_neo_m2)
#dropping logNInstitutions change the slopes of the intercep  >25%  -> keep in the model

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Step 3 building the final model ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#check if reintroducing logMaxLongevity changes the slopes of the final model
fit_reproduction_neo_m3 <- brm(NNeoplasia | trials(NNecropsies)  ~ logMaxLongevity + ClutchSize  + logNInstitutions + (1|gr(Species, cov = vcv_brms_OU)),
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

fit_reproduction_neo_m3 

diagnostic(reference_model =fit_reproduction_neo_m2, adjusted_model = fit_reproduction_neo_m3)
#adding logMaxLongevity does change the intercept by >25% ->  keep in the model


#check if reintroducing logMass changes the slopes of the final model
fit_reproduction_neo_m4 <- brm(NNeoplasia | trials(NNecropsies)  ~ logMaxLongevity + logMass + ClutchSize  + logNInstitutions + (1|gr(Species, cov = vcv_brms_OU)),
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

fit_reproduction_neo_m4 

diagnostic(reference_model =fit_reproduction_neo_m3, adjusted_model = fit_reproduction_neo_m4)
#adding back the logMass to the model has an effect on the logMaxLongevity  and ClutchSize by  >25% -> keep
  

#check if reintroducing IncubationLength changes the slopes of the final model
fit_reproduction_neo_m5 <- brm(NNeoplasia | trials(NNecropsies)  ~ logMaxLongevity + logMass + ClutchSize  + logNInstitutions + IncubationLength + (1|gr(Species, cov = vcv_brms_OU)),
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

fit_reproduction_neo_m5 

diagnostic(reference_model =fit_reproduction_neo_m4 , adjusted_model = fit_reproduction_neo_m5)
#adding IncubationLength  does not change the slopes of the other parameters by >25% -> drop from the model

#check if reintroducing logHIndex changes the slopes of the final model
fit_reproduction_neo_m6 <- brm(NNeoplasia | trials(NNecropsies)  ~ logMaxLongevity + logHIndex + logMass + ClutchSize  + logNInstitutions + (1|gr(Species, cov = vcv_brms_OU)),
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

fit_reproduction_neo_m6 

diagnostic(reference_model =fit_reproduction_neo_m4, adjusted_model = fit_reproduction_neo_m6)
#adding logHIndex does changes the slope of logMaxLongevity by >25% -> keep in the model

#check if reintroducing  ResidualslogNViews changes the slopes of the final model
fit_reproduction_neo_m7 <- brm(NNeoplasia | trials(NNecropsies)  ~ ResidualslogNViews + logMaxLongevity + logHIndex + logMass + ClutchSize  + logNInstitutions + (1|gr(Species, cov = vcv_brms_OU)),
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

fit_reproduction_neo_m7 

diagnostic(reference_model =fit_reproduction_neo_m6, adjusted_model = fit_reproduction_neo_m7)
#adding ResidualslogNViews  does not changes the slope of log h-index >25% -> keep in the model

#check if reintroducing   mean_GDP_5years changes the slopes of the final model
fit_reproduction_neo_m8 <- brm(NNeoplasia | trials(NNecropsies)  ~  mean_GDP_5years + ResidualslogNViews + logMaxLongevity + logHIndex + logMass + ClutchSize  + logNInstitutions + (1|gr(Species, cov = vcv_brms_OU)),
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

fit_reproduction_neo_m8 

diagnostic(reference_model =fit_reproduction_neo_m7, adjusted_model = fit_reproduction_neo_m8)
#adding mean_GDP_5years change the slopes of the other parameters by >25% -> keep in the model

final_model <- fit_reproduction_neo_m8 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Calculate and plot the marginal effects ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#create the conditions for the crude and adjusted models
cdt_crude_neo <-  list(NNecropsies = 100)

cdt_adj_neo <- list( logNInstitutions = median_logNInstitutions,
                     logMass = median_logMass,
                     logMaxLongevity  = median_logMaxLongevity ,
                     logHIndex = median_logHindex,
                     mean_GDP_5years = median_mean_GDP_5years,
                     ResidualslogNViews = median_ResidualslogNViews,
                      NNecropsies = 100)

marginal_effects_neo <- get_marginal_effects(crude_model = fit_reproduction_neo_ClutchSize , 
                                                                 adjusted_model = final_model ,
                                                                 resp = "NNeoplasia",
                                                                 effect = "ClutchSize",
                                                                 varext = "ClutchSize",
                                                                 NNecropsies = 100, 
                                                                 conditions_crude = cdt_crude_neo, 
                                                                 conditions_adjusted = cdt_adj_neo)

#plot the marginal effects for Clutch Size
df_neoplasia <- subset(marginal_effects_neo, ClutchSize >= min(na.omit(data_neo_mal$ClutchSize)) & ClutchSize <= max(na.omit(data_neo_mal$ClutchSize)))
p_ClutchSize_neo <- ggplot(data = df_neoplasia  , aes(x = ClutchSize, y = estimate__, group = Type, color = Type, fill = Type))+
  geom_line(size = 1.5)+
  geom_ribbon(mapping = aes(x = ClutchSize, ymin = lower__, ymax = upper__), color = NA, alpha = 0.25)+
  theme_classic()+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=20,face="bold"))+
  theme(legend.position = c(0.2, 0.5))+
  theme(plot.title = element_text(size = 18, face = "bold"))+
  theme(legend.text=element_text(size=12),
        legend.title = element_text(size = 12))+
  scale_color_manual(values = c("Crude" = "#E69F00", "Adjusted" = "#009E73")) +
  scale_fill_manual(values = c("Crude" = "#E69F00", "Adjusted" = "#009E73")) +
  #theme(axis.title.y = element_blank())+  # Remove axis titles
  xlab("Clutch Size")+
  ylab("Neoplasia prevalence")+
  labs(color = "Model type:", fill = "Model type:")+
  ylim(0,0.9)
#xlim(-3.5,3.5)
p_ClutchSize_neo 

####Add marginal effect for body mass

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Malignancy  ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Step 1 univariate models  ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#fit univariate model with logMass
fit_reproduction_mal_logMass <- brm(NMalignant | trials(NNecropsies)  ~ logMass + (1|gr(Species, cov = vcv_brms_OU)),
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
summary(fit_reproduction_mal_logMass, prob = 0.90)
p_direction(fit_reproduction_mal_logMass) 
#p_direction <0.90 -> exclude from the model - 

#fit univariate model with logMaxLongevity
fit_reproduction_mal_logMaxLongevity <- brm(NMalignant | trials(NNecropsies)  ~ logMaxLongevity + (1|gr(Species, cov = vcv_brms_OU)),
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
summary(fit_reproduction_mal_logMaxLongevity, prob = 0.90)
p_direction(fit_reproduction_mal_logMaxLongevity) 
#p_direction <0.90 -> exclude from the model - 

#fit univariate model with Incubation length
fit_reproduction_mal_IncubationLength <- brm(NMalignant | trials(NNecropsies)  ~ IncubationLength + (1|gr(Species, cov = vcv_brms_OU)),
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
summary(fit_reproduction_mal_IncubationLength, prob = 0.90)
p_direction(fit_reproduction_mal_IncubationLength) 
#p_direction <0.90 -> exclude from the model - 

#fit univariate model with ClutchSize
fit_reproduction_mal_ClutchSize <- brm(NMalignant | trials(NNecropsies)  ~ ClutchSize + (1|gr(Species, cov = vcv_brms_OU)),
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
summary(fit_reproduction_mal_ClutchSize, prob = 0.90)
p_direction(fit_reproduction_mal_ClutchSize) 
#p_direction >0.90 -> include in the model - 

#fit univariate model with logHIndex
fit_reproduction_mal_logHIndex <- brm(NMalignant | trials(NNecropsies)  ~ logHIndex + (1|gr(Species, cov = vcv_brms_OU)),
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
summary(fit_reproduction_mal_logHIndex, prob = 0.90)
p_direction(fit_reproduction_mal_logHIndex) 
#p_direction <0.90 -> exclude from the model - 


#fit univariate model with logNInstitutions
fit_reproduction_mal_logNInstitutions <- brm(NMalignant | trials(NNecropsies)  ~ logNInstitutions + (1|gr(Species, cov = vcv_brms_OU)),
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
summary(fit_reproduction_mal_logNInstitutions, prob = 0.90)
p_direction(fit_reproduction_mal_logNInstitutions) 
#p_direction >0.90 -> include in the model - 

#fit univariate model with mean_GDP_5years
fit_reproduction_mal_mean_GDP_5years <- brm(NMalignant | trials(NNecropsies)  ~ mean_GDP_5years + (1|gr(Species, cov = vcv_brms_OU)),
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
summary(fit_reproduction_mal_mean_GDP_5years, prob = 0.90)
p_direction(fit_reproduction_mal_mean_GDP_5years) 
#mean_GDP_5years is not significant


#fit univariate model with ResidualslogNViews
fit_reproduction_mal_ResidualslogNViews <- brm(NMalignant | trials(NNecropsies)  ~ ResidualslogNViews + (1|gr(Species, cov = vcv_brms_OU)),
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
summary(fit_reproduction_mal_ResidualslogNViews, prob = 0.90)
p_direction(fit_reproduction_mal_ResidualslogNViews) 
#ResidualslogNViews is significant

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Step 2 drop variables from the model  ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#fit a model with all variables retained in Step 1
fit_reproduction_mal_m1 <- brm(NMalignant | trials(NNecropsies)  ~ ClutchSize  + logNInstitutions + (1|gr(Species, cov = vcv_brms_OU)),
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
fit_reproduction_mal_m1
#logNExSitu has the lowest p-direction - drop it first

#drop logNInstitutions from the model
fit_reproduction_mal_m2 <- brm(NMalignant | trials(NNecropsies)  ~ ClutchSize + (1|gr(Species, cov = vcv_brms_OU)),
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
fit_reproduction_mal_m2

diagnostic(reference_model =fit_reproduction_mal_m1, adjusted_model = fit_reproduction_mal_m2)
#dropping logNInstitutions changes the intercept  >25%  -> keep in the model

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Step 3 building the final model ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#check if reintroducing logMaxLongevity changes the slopes of the final model
fit_reproduction_mal_m3 <- brm(NMalignant | trials(NNecropsies)  ~ logMaxLongevity + ClutchSize  + logNInstitutions + (1|gr(Species, cov = vcv_brms_OU)),
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

fit_reproduction_mal_m3 

diagnostic(reference_model =fit_reproduction_mal_m1, adjusted_model = fit_reproduction_mal_m3)
#adding logMaxLongevity change the slopes of the other parameters by >25% -> add back to the model


#check if reintroducing logMass changes the slopes of the final model
fit_reproduction_mal_m4 <- brm(NMalignant | trials(NNecropsies)  ~ logMass + logMaxLongevity + ClutchSize  + logNInstitutions + (1|gr(Species, cov = vcv_brms_OU)),
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

fit_reproduction_mal_m4 

diagnostic(reference_model =fit_reproduction_mal_m3, adjusted_model = fit_reproduction_mal_m4)
#adding back the logMass to the model has an effect on the ClutchSize  slope by  >25% -> keep


#check if reintroducing IncubationLength changes the slopes of the final model
fit_reproduction_mal_m5 <- brm(NMalignant | trials(NNecropsies)  ~ IncubationLength + logMass + logMaxLongevity + ClutchSize  + logNInstitutions + (1|gr(Species, cov = vcv_brms_OU)),
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

fit_reproduction_mal_m5 

diagnostic(reference_model =fit_reproduction_mal_m4 , adjusted_model = fit_reproduction_mal_m5)
#adding IncubationLength  does not change the slopes of the other parameters by >25% -> drop from the model

#check if reintroducing logHIndex changes the slopes of the final model
fit_reproduction_mal_m6 <- brm(NMalignant | trials(NNecropsies)  ~ logHIndex +  logMass + logMaxLongevity + ClutchSize  + logNInstitutions + (1|gr(Species, cov = vcv_brms_OU)),
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

fit_reproduction_mal_m6 

diagnostic(reference_model =fit_reproduction_mal_m4, adjusted_model = fit_reproduction_mal_m6)
#adding logHIndex doesn't change the slopes of the other parameters by >25% -> drop from model

#check if reintroducing  ResidualslogNViews changes the slopes of the final model
fit_reproduction_mal_m7 <- brm(NMalignant | trials(NNecropsies)  ~ ResidualslogNViews +  logMass + logMaxLongevity + ClutchSize  + logNInstitutions + (1|gr(Species, cov = vcv_brms_OU)),
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

fit_reproduction_mal_m7 

diagnostic(reference_model =fit_reproduction_mal_m4, adjusted_model = fit_reproduction_mal_m7)
#adding ResidualslogNViews does change the slopes of logMaxLongevity  by >25% -> drop from the model

#check if reintroducing   mean_GDP_5years changes the slopes of the final model
fit_reproduction_mal_m8 <- brm(NMalignant | trials(NNecropsies)  ~  mean_GDP_5years + ResidualslogNViews +  logMass + logMaxLongevity + ClutchSize  + logNInstitutions +   (1|gr(Species, cov = vcv_brms_OU)),
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

fit_reproduction_mal_m8

diagnostic(reference_model =fit_reproduction_mal_m7, adjusted_model = fit_reproduction_mal_m8)
#adding mean_GDP_5years change the slopes of the other parameters by >25% -> keep in the model

fit_reproduction_mal_final <- fit_reproduction_mal_m8


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Calculate and plot the marginal effects ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#create the conditions for the crude and adjusted models
cdt_crude_mal <-  list(NNecropsies = 100)

cdt_adj_mal <- list( logNInstitutions = median_logNInstitutions,
                     logMass = median_logMass,
                     logMaxLongevity = median_logMaxLongevity,
                     ResidualslogNViews = median_ResidualslogNViews,
                     mean_GDP_5years= median_mean_GDP_5years,
                     NNecropsies = 100)

marginal_effects_mal <- get_marginal_effects(crude_model = fit_reproduction_mal_ClutchSize , 
                                                       adjusted_model = fit_reproduction_mal_final,
                                                       resp = "NMalignant",
                                                       effect = "ClutchSize",
                                                       varext = "ClutchSize",
                                                       NNecropsies = 100, 
                                                       conditions_crude = cdt_crude_mal, 
                                                       conditions_adjusted = cdt_adj_mal)

#plot the marginal effects for mutation rate
df_mortality <- subset(marginal_effects_mal, ClutchSize >= min(na.omit(data_neo_mal$ClutchSize)) & ClutchSize <= max(na.omit(data_neo_mal$ClutchSize)))
p_ClutchSize_mal <- ggplot(data = df_mortality , aes(x = ClutchSize, y = estimate__, group = Type, color = Type, fill = Type))+
  geom_line(size = 1.5)+
  geom_ribbon(mapping = aes(x = ClutchSize, ymin = lower__, ymax = upper__), color = NA, alpha = 0.25)+
  theme_classic()+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=20,face="bold"))+
  theme(legend.position = c(0.2, 0.5))+
  theme(plot.title = element_text(size = 18, face = "bold"))+
  theme(legend.text=element_text(size=12),
        legend.title = element_text(size = 12))+
  scale_color_manual(values = c("Crude" = "#E69F00", "Adjusted" = "#009E73")) +
  scale_fill_manual(values = c("Crude" = "#E69F00", "Adjusted" = "#009E73")) +
  #theme(axis.title.y = element_blank())+  # Remove axis titles
  xlab("Clutch Size")+
  ylab("Malignancy prevalence")+
  labs(color = "Model type:", fill = "Model type:")+
  ylim(0,1)
#xlim(-3.5,3.5)
p_ClutchSize_mal 

