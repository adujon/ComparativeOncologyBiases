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
setwd("C:/Users/antoi/My Drive/Postdoc/Zoo data/Biases/Scripts/Estimation of confounding effects/for github")

data_neo_mal_mass <- data_neo_mal<- read_excel("Raw datasets.xlsx", sheet = "Peto_neoplasia_malignancy")
rownames(data_neo_mal_mass) <- data_neo_mal_mass$Species

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import Mortality data----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

data_mortality <- read_excel("Raw datasets.xlsx", sheet = "Peto_mortality")
rownames(data_mortality) <- data_mortality$Species

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Calculate median for bilbiometric variabels to use in calculation or marinal effects ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#calculate the median of confounding variables to use in the calculation of marginal effects
median_logHindex <- median(na.omit(c(data_neo_mal$logHIndex,data_mortality$logHIndex)))
median_logNInstitutions <- median(na.omit(c(data_neo_mal$logNInstitutions, data_mortality$logNInstitutions)))
median_logNExSitu <- median(na.omit(c(data_neo_mal$logNExSitu,  data_mortality$logNExSitu)))
median_ResidualslogNViews <- median(na.omit(c(data_neo_mal$ResidualslogNViews,  data_mortality$ResidualslogNViews)))
median_logNExSitu <- median(na.omit(c(data_neo_mal$logNExSitu,  data_mortality$logNExSitu)))
median_logMass <- median(na.omit(c(data_neo_mal$logMass,  data_mortality$logMass)))
median_logMaxLongevity <- median(na.omit(c(data_neo_mal$logMaxLongevity,  data_mortality$logMaxLongevity)))
median_mean_GDP_5years <- median(na.omit(c(data_neo_mal$mean_GDP_5years,  data_mortality$mean_GDP_5years)))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Prepare var covar matrixes for brms models ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#check which species are both in the main dataframe and the phylogeny tree
obj_brms_neo_mal<- name.check(phylo_all, data_neo_mal_mass)
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

#univariate model with log mass
fit_peto_neo_logMass <-  brm(NNeoplasia | trials(NNecropsies)  ~ logMass + (1|gr(Species, cov =  vcv_brms_OU)),
                               family =binomial(),
                               data = data_neo_mal_mass,
                               data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                               chains = 4,
                               cores = 24,
                               iter = 10000, warmup = 5000,
                               thin = 5,
                               seed = 1,
                               control=list(max_treedepth = 12, adapt_delta = 0.99),
                               save_pars = save_pars(all = TRUE))
fit_peto_neo_logMass
summary(fit_peto_neo_logMass,prob = 0.90)
p_direction(fit_peto_neo_logMass)
#p-direction >0.90 -> keep in the model


#univariate model with log max longevity
fit_peto_neo_logMaxLongevity <-  brm(NNeoplasia | trials(NNecropsies)  ~ logMaxLongevity + (1|gr(Species, cov =  vcv_brms_OU)),
                             family =binomial(),
                             data = data_neo_mal_mass,
                             data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                             chains = 4,
                             cores = 24,
                             iter = 10000, warmup = 5000,
                             thin = 5,
                             seed = 1,
                             control=list(max_treedepth = 12, adapt_delta = 0.99),
                             save_pars = save_pars(all = TRUE))
fit_peto_neo_logMaxLongevity
summary(fit_peto_neo_logMaxLongevity,prob = 0.90)
p_direction(fit_peto_neo_logMaxLongevity)
#p-direction >0.90 -> keep in the model

#univariate model with class
fit_peto_neo_Class <-  brm(NNeoplasia | trials(NNecropsies)  ~ Class + (1|gr(Species, cov =  vcv_brms_OU)),
                                     family =binomial(),
                                     data = data_neo_mal_mass,
                                     data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                     chains = 4,
                                     cores = 24,
                                     iter = 10000, warmup = 5000,
                                     thin = 5,
                                     seed = 1,
                                     control=list(max_treedepth = 12, adapt_delta = 0.99),
                                     save_pars = save_pars(all = TRUE))
fit_peto_neo_Class
summary(fit_peto_neo_Class,prob = 0.90)
p_direction(fit_peto_neo_Class)
#p-direction >0.90 -> keep in the model

#univariate model with logHIndex
fit_peto_neo_logHIndex <-  brm(NNeoplasia | trials(NNecropsies)  ~ logHIndex + (1|gr(Species, cov =  vcv_brms_OU)),
                           family =binomial(),
                           data = data_neo_mal_mass,
                           data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                           chains = 4,
                           cores = 24,
                           iter = 10000, warmup = 5000,
                           thin = 5,
                           seed = 1,
                           control=list(max_treedepth = 12, adapt_delta = 0.99),
                           save_pars = save_pars(all = TRUE))
fit_peto_neo_logHIndex
summary(fit_peto_neo_logHIndex,prob = 0.90)
p_direction(fit_peto_neo_logHIndex)
#p-direction >0.90 -> keep in the model

#univariate model with logNInstitutions
fit_peto_neo_logNInstitutions<-  brm(NNeoplasia | trials(NNecropsies)  ~ logHIndex + (1|gr(Species, cov =  vcv_brms_OU)),
                               family =binomial(),
                               data = data_neo_mal_mass,
                               data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                               chains = 4,
                               cores = 24,
                               iter = 10000, warmup = 5000,
                               thin = 5,
                               seed = 1,
                               control=list(max_treedepth = 12, adapt_delta = 0.99),
                               save_pars = save_pars(all = TRUE))
fit_peto_neo_logNInstitutions
summary(fit_peto_neo_logNInstitutions,prob = 0.90)
p_direction(fit_peto_neo_logNInstitutions)
#p-direction >0.90 -> keep in the model

#univariate model with mean_GDP_5years
fit_peto_neo_mean_GDP_5years<-  brm(NNeoplasia | trials(NNecropsies)  ~ mean_GDP_5years + (1|gr(Species, cov =  vcv_brms_OU)),
                                     family =binomial(),
                                     data = data_neo_mal_mass,
                                     data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                     chains = 4,
                                     cores = 24,
                                     iter = 10000, warmup = 5000,
                                     thin = 5,
                                     seed = 1,
                                     control=list(max_treedepth = 12, adapt_delta = 0.99),
                                     save_pars = save_pars(all = TRUE))
fit_peto_neo_mean_GDP_5years
summary(fit_peto_neo_mean_GDP_5years,prob = 0.90)
p_direction(fit_peto_neo_mean_GDP_5years)


#univariate model with ResidualslogNViews 
fit_peto_neo_ResidualslogNViews <-  brm(NNeoplasia | trials(NNecropsies)  ~ ResidualslogNViews  + (1|gr(Species, cov =  vcv_brms_OU)),
                                     family =binomial(),
                                     data = data_neo_mal_mass,
                                     data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                     chains = 4,
                                     cores = 24,
                                     iter = 10000, warmup = 5000,
                                     thin = 5,
                                     seed = 1,
                                     control=list(max_treedepth = 12, adapt_delta = 0.99),
                                     save_pars = save_pars(all = TRUE))
fit_peto_neo_ResidualslogNViews 
summary(fit_peto_neo_ResidualslogNViews ,prob = 0.90)
p_direction(fit_peto_neo_ResidualslogNViews )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Step 2: Drop part of the model selection ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#fit a model with all the variables retained in step 1
fit_peto_neo_m1 <- brm(NNeoplasia | trials(NNecropsies)  ~logMaxLongevity + logMass + Class + logHIndex + logNInstitutions  + mean_GDP_5years + ResidualslogNViews + (1|gr(Species, cov =  vcv_brms_OU)),
                       family =binomial(),
                       data = data_neo_mal_mass,
                       data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                       chains = 4,
                       cores = 24,
                       iter = 10000, warmup = 5000,
                       thin = 5,
                       seed = 1,
                       control=list(max_treedepth = 12, adapt_delta = 0.99),
                       save_pars = save_pars(all = TRUE))
fit_peto_neo_m1
p_direction(fit_peto_neo_m1) 
#drop mean_GDP_5years  first

#drop mean_GDP_5years  from the model
fit_peto_neo_m2 <- brm(NNeoplasia | trials(NNecropsies)  ~logMaxLongevity + logMass + Class + logHIndex + logNInstitutions + ResidualslogNViews + (1|gr(Species, cov =  vcv_brms_OU)),
                       family =binomial(),
                       data = data_neo_mal_mass,
                       data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                       chains = 4,
                       cores = 24,
                       iter = 10000, warmup = 5000,
                       thin = 5,
                       seed = 1,
                       control=list(max_treedepth = 12, adapt_delta = 0.99),
                       save_pars = save_pars(all = TRUE))
fit_peto_neo_m2
diagnostic(reference_model =fit_peto_neo_m1, adjusted_model = fit_peto_neo_m2)
p_direction(fit_peto_neo_m2) 
#mean_GDP_5years doesn't change the slope of other paramater by>25 % -> drop from the model
#drop class nest

#drop class  from the model
fit_peto_neo_m3 <- brm(NNeoplasia | trials(NNecropsies)  ~logMaxLongevity + logMass + logHIndex + logNInstitutions + ResidualslogNViews + (1|gr(Species, cov =  vcv_brms_OU)),
                       family =binomial(),
                       data = data_neo_mal_mass,
                       data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                       chains = 4,
                       cores = 24,
                       iter = 10000, warmup = 5000,
                       thin = 5,
                       seed = 1,
                       control=list(max_treedepth = 12, adapt_delta = 0.99),
                       save_pars = save_pars(all = TRUE))
fit_peto_neo_m3
diagnostic(reference_model =fit_peto_neo_m2, adjusted_model = fit_peto_neo_m3)
p_direction(fit_peto_neo_m3) 
#dropping class change the slope of logMass and logHIndex > 25% -> keep in the model
#drop longevity next from the model

#drop logMaxLongevity  from the model
fit_peto_neo_m4 <- brm(NNeoplasia | trials(NNecropsies)  ~logMass + Class + logHIndex + logNInstitutions + ResidualslogNViews + (1|gr(Species, cov =  vcv_brms_OU)),
                       family =binomial(),
                       data = data_neo_mal_mass,
                       data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                       chains = 4,
                       cores = 24,
                       iter = 10000, warmup = 5000,
                       thin = 5,
                       seed = 1,
                       control=list(max_treedepth = 12, adapt_delta = 0.99),
                       save_pars = save_pars(all = TRUE))
fit_peto_neo_m4
diagnostic(reference_model =fit_peto_neo_m2, adjusted_model = fit_peto_neo_m4)
p_direction(fit_peto_neo_m4) 
#droping longevity change the slope of logMass >25%
#drop logHIndex next

#drop logHIndex from the model
fit_peto_neo_m5 <- brm(NNeoplasia | trials(NNecropsies)  ~logMaxLongevity + logMass + Class  + logNInstitutions + ResidualslogNViews + (1|gr(Species, cov =  vcv_brms_OU)),
                       family =binomial(),
                       data = data_neo_mal_mass,
                       data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                       chains = 4,
                       cores = 24,
                       iter = 10000, warmup = 5000,
                       thin = 5,
                       seed = 1,
                       control=list(max_treedepth = 12, adapt_delta = 0.99),
                       save_pars = save_pars(all = TRUE))
fit_peto_neo_m5
diagnostic(reference_model =fit_peto_neo_m4, adjusted_model = fit_peto_neo_m5)
p_direction(fit_peto_neo_m5) 
#log logHIndex  change the slope of ClassAves  and logMaxLongevity >25%

#drop logNInstitutions from the model
fit_peto_neo_m6 <- brm(NNeoplasia | trials(NNecropsies)  ~logMaxLongevity + logMass + Class + logHIndex  + ResidualslogNViews + (1|gr(Species, cov =  vcv_brms_OU)),
                       family =binomial(),
                       data = data_neo_mal_mass,
                       data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                       chains = 4,
                       cores = 24,
                       iter = 10000, warmup = 5000,
                       thin = 5,
                       seed = 1,
                       control=list(max_treedepth = 12, adapt_delta = 0.99),
                       save_pars = save_pars(all = TRUE))
fit_peto_neo_m6
diagnostic(reference_model =fit_peto_neo_m4, adjusted_model = fit_peto_neo_m6)
p_direction(fit_peto_neo_m6) 
#dropping logNInstutition change the slope of logMaxLongevity and ClassAves  >25% -> keep in the model

#drop ResidualslogNViews from the model
fit_peto_neo_m7 <- brm(NNeoplasia | trials(NNecropsies)  ~logMaxLongevity + logMass + Class + logHIndex + logNInstitutions + (1|gr(Species, cov =  vcv_brms_OU)),
                       family =binomial(),
                       data = data_neo_mal_mass,
                       data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                       chains = 4,
                       cores = 24,
                       iter = 10000, warmup = 5000,
                       thin = 5,
                       seed = 1,
                       control=list(max_treedepth = 12, adapt_delta = 0.99),
                       save_pars = save_pars(all = TRUE))
fit_peto_neo_m7
diagnostic(reference_model =fit_peto_neo_m4, adjusted_model = fit_peto_neo_m7)
p_direction(fit_peto_neo_m7) 
#dropping ResidualslogNViews change the slope of logMaxLongevity and ClassAves  >25% -> keep in the model



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Step 3: Final model ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


fit_peto_neo_m8 <- brm(NNeoplasia | trials(NNecropsies)  ~logMaxLongevity * logMass + Class + logHIndex + logNInstitutions + ResidualslogNViews + (1|gr(Species, cov =  vcv_brms_OU)),
                       family =binomial(),
                       data = data_neo_mal_mass,
                       data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                       chains = 4,
                       cores = 24,
                       iter = 10000, warmup = 5000,
                       thin = 5,
                       seed = 1,
                       control=list(max_treedepth = 12, adapt_delta = 0.99),
                       save_pars = save_pars(all = TRUE))
fit_peto_neo_m8
diagnostic(reference_model =fit_peto_neo_m2, adjusted_model = fit_peto_neo_m8)
#adding an interaction between logmass and longevity changes the slope of logmass >25%


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Calculate and plot the marginal effects ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

crude_model_neo <- brm(NNeoplasia | trials(NNecropsies)  ~ logMaxLongevity * logMass + (1|gr(Species, cov =  vcv_brms_OU)),
                   family =binomial(),
                   data = data_neo_mal_mass,
                   data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                   chains = 4,
                   cores = 24,
                   iter = 10000, warmup = 5000,
                   thin = 5,
                   seed = 1,
                   control=list(max_treedepth = 12, adapt_delta = 0.99),
                   save_pars = save_pars(all = TRUE))
  
#create the conditions for the crude and adjusted models
cdt_crude_peto_neo <-  list(NNecropsies = 100,
                            logMaxLongevity = median_logMaxLongevity)

cdt_adj_peto_neo <- list(logMass = median_logMass,
                             logHIndex = median_logHindex, 
                             logMaxLongevity = median_logMaxLongevity,
                             logNInstitutions = median_logNInstitutions,
                             ResidualslogNViews = median_ResidualslogNViews,
                             Class = "Aves",
                             NNecropsies = 100)

marginal_effects_peto_neo <- get_marginal_effects(crude_model = crude_model_neo, 
                                                       adjusted_model = fit_peto_neo_m8,
                                                       resp = "NNeoplasia",
                                                       effect = "logMass",
                                                       varext = "logMass",
                                                       NNecropsies = 100, 
                                                       conditions_crude = cdt_crude_peto_neo , 
                                                       conditions_adjusted = cdt_adj_peto_neo)

#plot the marginal effects for mutation rate
df_peto_neo <- subset(marginal_effects_peto_neo, logMass >= min(data_neo_mal_mass$logMass) & logMass <= max(data_neo_mal_mass$logMass))
p_peto_neo <- ggplot(data = df_peto_neo , aes(x = logMass, y = estimate__, group = Type, color = Type, fill = Type))+
  geom_line(size = 1.5)+
  geom_ribbon(mapping = aes(x = logMass, ymin = lower__, ymax = upper__), color = NA, alpha = 0.25)+
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
  xlab("log10 body mass (g)")+
  ylab("Neoplasia prevalence")+
  labs(color = "Model type:", fill = "Model type:")+
  ylim(0,0.5)
#xlim(-3.5,3.5)
p_peto_neo

###add plot for longevity marginal effects

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Malignancy ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Step 1: Univariate models ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#univariate model with log mass
fit_peto_mal_logMass <-  brm(NMalignant | trials(NNecropsies)  ~ logMass + (1|gr(Species, cov =  vcv_brms_OU)),
                             family =binomial(),
                             data = data_neo_mal_mass,
                             data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                             chains = 4,
                             cores = 24,
                             iter = 10000, warmup = 5000,
                             thin = 5,
                             seed = 1,
                             control=list(max_treedepth = 12, adapt_delta = 0.99),
                             save_pars = save_pars(all = TRUE))
fit_peto_mal_logMass
summary(fit_peto_mal_logMass,prob = 0.90)
p_direction(fit_peto_mal_logMass)
#p-direction >0.90 -> keep in the model


#univariate model with log max longevity
fit_peto_mal_logMaxLongevity <-  brm(NMalignant | trials(NNecropsies)  ~ logMaxLongevity + (1|gr(Species, cov =  vcv_brms_OU)),
                                     family =binomial(),
                                     data = data_neo_mal_mass,
                                     data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                     chains = 4,
                                     cores = 24,
                                     iter = 10000, warmup = 5000,
                                     thin = 5,
                                     seed = 1,
                                     control=list(max_treedepth = 12, adapt_delta = 0.99),
                                     save_pars = save_pars(all = TRUE))
fit_peto_mal_logMaxLongevity
summary(fit_peto_mal_logMaxLongevity,prob = 0.90)
p_direction(fit_peto_mal_logMaxLongevity)
#p-direction >0.90 -> keep in the model

#univariate model with class
fit_peto_mal_Class <-  brm(NMalignant | trials(NNecropsies)  ~ Class + (1|gr(Species, cov =  vcv_brms_OU)),
                           family =binomial(),
                           data = data_neo_mal_mass,
                           data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                           chains = 4,
                           cores = 24,
                           iter = 10000, warmup = 5000,
                           thin = 5,
                           seed = 1,
                           control=list(max_treedepth = 12, adapt_delta = 0.99),
                           save_pars = save_pars(all = TRUE))
fit_peto_mal_Class
summary(fit_peto_mal_Class,prob = 0.90)
p_direction(fit_peto_mal_Class)
#p-direction >0.90 -> keep in the model

#univariate model with logHIndex
fit_peto_mal_logHIndex <-  brm(NMalignant | trials(NNecropsies)  ~ logHIndex + (1|gr(Species, cov =  vcv_brms_OU)),
                               family =binomial(),
                               data = data_neo_mal_mass,
                               data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                               chains = 4,
                               cores = 24,
                               iter = 10000, warmup = 5000,
                               thin = 5,
                               seed = 1,
                               control=list(max_treedepth = 12, adapt_delta = 0.99),
                               save_pars = save_pars(all = TRUE))
fit_peto_mal_logHIndex
summary(fit_peto_mal_logHIndex,prob = 0.90)
p_direction(fit_peto_mal_logHIndex)
#p-direction >0.90 -> keep in the model

#univariate model with logNInstitutions
fit_peto_mal_logNInstitutions<-  brm(NMalignant | trials(NNecropsies)  ~ logHIndex + (1|gr(Species, cov =  vcv_brms_OU)),
                                     family =binomial(),
                                     data = data_neo_mal_mass,
                                     data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                     chains = 4,
                                     cores = 24,
                                     iter = 10000, warmup = 5000,
                                     thin = 5,
                                     seed = 1,
                                     control=list(max_treedepth = 12, adapt_delta = 0.99),
                                     save_pars = save_pars(all = TRUE))
fit_peto_mal_logNInstitutions
summary(fit_peto_mal_logNInstitutions,prob = 0.90)
p_direction(fit_peto_mal_logNInstitutions)
#p-direction >0.90 -> keep in the model

#univariate model with mean_GDP_5years
fit_peto_mal_mean_GDP_5years<-  brm(NMalignant | trials(NNecropsies)  ~ mean_GDP_5years + (1|gr(Species, cov =  vcv_brms_OU)),
                                    family =binomial(),
                                    data = data_neo_mal_mass,
                                    data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                    chains = 4,
                                    cores = 24,
                                    iter = 10000, warmup = 5000,
                                    thin = 5,
                                    seed = 1,
                                    control=list(max_treedepth = 12, adapt_delta = 0.99),
                                    save_pars = save_pars(all = TRUE))
fit_peto_mal_mean_GDP_5years
summary(fit_peto_mal_mean_GDP_5years,prob = 0.90)
p_direction(fit_peto_mal_mean_GDP_5years)
#p-direction >0.90 -> keep in the model

#univariate model with ResidualslogNViews 
fit_peto_mal_ResidualslogNViews <-  brm(NMalignant | trials(NNecropsies)  ~ ResidualslogNViews  + (1|gr(Species, cov =  vcv_brms_OU)),
                                        family =binomial(),
                                        data = data_neo_mal_mass,
                                        data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                                        chains = 4,
                                        cores = 24,
                                        iter = 10000, warmup = 5000,
                                        thin = 5,
                                        seed = 1,
                                        control=list(max_treedepth = 12, adapt_delta = 0.99),
                                        save_pars = save_pars(all = TRUE))
fit_peto_mal_ResidualslogNViews 
summary(fit_peto_mal_ResidualslogNViews ,prob = 0.90)
p_direction(fit_peto_mal_ResidualslogNViews )
#p-direction >0.90 -> keep in the model


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Step 2: Drop part of the model selection ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#fit a model with all the variables retained in step 1
fit_peto_mal_m1 <- brm(NMalignant | trials(NNecropsies)  ~ logMaxLongevity + logMass + Class + logHIndex + logNInstitutions  + mean_GDP_5years + ResidualslogNViews + (1|gr(Species, cov =  vcv_brms_OU)),
                       family =binomial(),
                       data = data_neo_mal_mass,
                       data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                       chains = 4,
                       cores = 24,
                       iter = 10000, warmup = 5000,
                       thin = 5,
                       seed = 1,
                       control=list(max_treedepth = 12, adapt_delta = 0.99),
                       save_pars = save_pars(all = TRUE))
fit_peto_mal_m1
p_direction(fit_peto_mal_m1) 
#drop mean_GDP_5years  first


#model with mean_GDP_5years dropped
fit_peto_mal_m2 <- brm(NMalignant | trials(NNecropsies)  ~ logMaxLongevity + logMass + Class + logHIndex + logNInstitutions  + ResidualslogNViews + (1|gr(Species, cov =  vcv_brms_OU)),
                       family =binomial(),
                       data = data_neo_mal_mass,
                       data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                       chains = 4,
                       cores = 24,
                       iter = 10000, warmup = 5000,
                       thin = 5,
                       seed = 1,
                       control=list(max_treedepth = 12, adapt_delta = 0.99),
                       save_pars = save_pars(all = TRUE))
fit_peto_mal_m2
diagnostic(reference_model =fit_peto_mal_m1, adjusted_model = fit_peto_mal_m2)
p_direction(fit_peto_mal_m2)
#droping mean_GDP_5years doesn't change the slopes >25% - drop from the model
#drop class next

#drop class from the model
fit_peto_mal_m3 <- brm(NMalignant | trials(NNecropsies)  ~ logMaxLongevity + logMass + logHIndex + logNInstitutions  + ResidualslogNViews + (1|gr(Species, cov =  vcv_brms_OU)),
                       family =binomial(),
                       data = data_neo_mal_mass,
                       data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                       chains = 4,
                       cores = 24,
                       iter = 10000, warmup = 5000,
                       thin = 5,
                       seed = 1,
                       control=list(max_treedepth = 12, adapt_delta = 0.99),
                       save_pars = save_pars(all = TRUE))
fit_peto_mal_m3
diagnostic(reference_model =fit_peto_mal_m2, adjusted_model = fit_peto_mal_m3)
p_direction(fit_peto_mal_m3)
#excluding class change the slope of logMass >25% -> keep in the model
#drop logHIndex next

#model with logHIndex dropped
fit_peto_mal_m4 <- brm(NMalignant | trials(NNecropsies)  ~ logMaxLongevity + logMass + Class + logNInstitutions  + ResidualslogNViews + (1|gr(Species, cov =  vcv_brms_OU)),
                       family =binomial(),
                       data = data_neo_mal_mass,
                       data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                       chains = 4,
                       cores = 24,
                       iter = 10000, warmup = 5000,
                       thin = 5,
                       seed = 1,
                       control=list(max_treedepth = 12, adapt_delta = 0.99),
                       save_pars = save_pars(all = TRUE))
fit_peto_mal_m4
diagnostic(reference_model =fit_peto_mal_m2, adjusted_model = fit_peto_mal_m4)
p_direction(fit_peto_mal_m4)
#droping logHIndex change the slopes of Class and  logMaxLongevity >25% - keep in the model
#drop logNInstitutions next

#drop logNInstitutions from the model
fit_peto_mal_m5 <- brm(NMalignant | trials(NNecropsies)  ~ logMaxLongevity + logMass + Class + logHIndex  + ResidualslogNViews + (1|gr(Species, cov =  vcv_brms_OU)),
                       family =binomial(),
                       data = data_neo_mal_mass,
                       data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                       chains = 4,
                       cores = 24,
                       iter = 10000, warmup = 5000,
                       thin = 5,
                       seed = 1,
                       control=list(max_treedepth = 12, adapt_delta = 0.99),
                       save_pars = save_pars(all = TRUE))
fit_peto_mal_m5
diagnostic(reference_model =fit_peto_mal_m2, adjusted_model = fit_peto_mal_m5)
p_direction(fit_peto_mal_m5)
#drop ResidualslogNViews changes the slope of logHIndex, class and logMaxLongevity >2% -> keep in the model

#model with ResidualslogNViews dropped
fit_peto_mal_m6 <- brm(NMalignant | trials(NNecropsies)  ~ logMaxLongevity + logMass + Class + logHIndex + logNInstitutions + (1|gr(Species, cov =  vcv_brms_OU)),
                       family =binomial(),
                       data = data_neo_mal_mass,
                       data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                       chains = 4,
                       cores = 24,
                       iter = 10000, warmup = 5000,
                       thin = 5,
                       seed = 1,
                       control=list(max_treedepth = 12, adapt_delta = 0.99),
                       save_pars = save_pars(all = TRUE))
fit_peto_mal_m6
diagnostic(reference_model =fit_peto_mal_m2, adjusted_model = fit_peto_mal_m6)
p_direction(fit_peto_mal_m6)
#droping ResidualslogNViews change the slopes of class and lofMass >25% - keep in the model

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Step 3: Final model ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#check for an interaction between mass and longevity
fit_peto_mal_m7 <- brm(NMalignant | trials(NNecropsies)  ~ logMaxLongevity * logMass + Class + logHIndex + logNInstitutions  + ResidualslogNViews + (1|gr(Species, cov =  vcv_brms_OU)),
                       family =binomial(),
                       data = data_neo_mal_mass,
                       data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                       chains = 4,
                       cores = 24,
                       iter = 10000, warmup = 5000,
                       thin = 5,
                       seed = 1,
                       control=list(max_treedepth = 12, adapt_delta = 0.99),
                       save_pars = save_pars(all = TRUE))

diagnostic(reference_model =fit_peto_mal_m2, adjusted_model = fit_peto_mal_m7)
p_direction(fit_peto_mal_m7)
#adding an interaction term changes the slope of logmass >25% - keep in the model

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Calculate and plot the marginal effects ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

crude_model_mal <- brm(NMalignant | trials(NNecropsies)  ~ logMaxLongevity * logMass + (1|gr(Species, cov =  vcv_brms_OU)),
                       family =binomial(),
                       data = data_neo_mal_mass,
                       data2 = list(vcv_brms_OU = vcv_brms_OU_neo_mal), 
                       chains = 4,
                       cores = 24,
                       iter = 10000, warmup = 5000,
                       thin = 5,
                       seed = 1,
                       control=list(max_treedepth = 12, adapt_delta = 0.99),
                       save_pars = save_pars(all = TRUE))

#create the conditions for the crude and adjusted models
cdt_crude_peto_mal <-  list(NNecropsies = 100,
                            logMaxLongevity = median_logMaxLongevity)

cdt_adj_peto_mal <- list(logMass = median_logMass,
                         logHIndex = median_logHindex, 
                         logMaxLongevity = median_logMaxLongevity,
                         logNInstitutions = median_logNInstitutions,
                         ResidualslogNViews = median_ResidualslogNViews,
                         Class = "Mammalia",
                         NNecropsies = 100)

marginal_effects_peto_mal <- get_marginal_effects(crude_model = crude_model_mal, 
                                                  adjusted_model = fit_peto_mal_m7,
                                                  resp = "NNeoplasia",
                                                  effect = "logMass",
                                                  varext = "logMass",
                                                  NNecropsies = 100, 
                                                  conditions_crude = cdt_crude_peto_mal , 
                                                  conditions_adjusted = cdt_adj_peto_mal)

#plot the marginal effects for mutation rate
df_peto_mal <- subset(marginal_effects_peto_mal, logMass >= min(data_neo_mal_mass$logMass) & logMass <= max(data_neo_mal_mass$logMass))
p_peto_mal <- ggplot(data = df_peto_mal , aes(x = logMass, y = estimate__, group = Type, color = Type, fill = Type))+
  geom_line(size = 1.5)+
  geom_ribbon(mapping = aes(x = logMass, ymin = lower__, ymax = upper__), color = NA, alpha = 0.25)+
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
  xlab("log10 body mass (g)")+
  ylab("Malignancy prevalence")+
  labs(color = "Model type:", fill = "Model type:")+
  ylim(0,0.5)
#xlim(-3.5,3.5)
p_peto_mal

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Malignancy ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Step 1: Univariate models ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#univariate model with log mass
fit_peto_mor_logMass <-  brm(NMortality | trials(NNecropsies)  ~ logMass + (1|gr(Species, cov =  vcv_brms_OU)),
                             family =binomial(),
                             data = data_mortality,
                             data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                             chains = 4,
                             cores = 24,
                             iter = 10000, warmup = 5000,
                             thin = 5,
                             seed = 1,
                             control=list(max_treedepth = 12, adapt_delta = 0.99),
                             save_pars = save_pars(all = TRUE))
fit_peto_mor_logMass
summary(fit_peto_mor_logMass,prob = 0.90)
p_direction(fit_peto_mor_logMass)
#p-direction >0.90 -> keep in the model


#univariate model with log max longevity
fit_peto_mor_logMaxLongevity <-  brm(NMortality | trials(NNecropsies)  ~ logMaxLongevity + (1|gr(Species, cov =  vcv_brms_OU)),
                                     family =binomial(),
                                     data = data_mortality,
                                     data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                                     chains = 4,
                                     cores = 24,
                                     iter = 10000, warmup = 5000,
                                     thin = 5,
                                     seed = 1,
                                     control=list(max_treedepth = 12, adapt_delta = 0.99),
                                     save_pars = save_pars(all = TRUE))
fit_peto_mor_logMaxLongevity
summary(fit_peto_mor_logMaxLongevity,prob = 0.90)
p_direction(fit_peto_mor_logMaxLongevity)
#p-direction <0.90 -> keep in the model because Peto

#univariate model with logHIndex
fit_peto_mor_logHIndex <-  brm(NMortality | trials(NNecropsies)  ~ logHIndex + (1|gr(Species, cov =  vcv_brms_OU)),
                               family =binomial(),
                               data = data_mortality,
                               data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                               chains = 4,
                               cores = 24,
                               iter = 10000, warmup = 5000,
                               thin = 5,
                               seed = 1,
                               control=list(max_treedepth = 12, adapt_delta = 0.99),
                               save_pars = save_pars(all = TRUE))
fit_peto_mor_logHIndex
summary(fit_peto_mor_logHIndex,prob = 0.90)
p_direction(fit_peto_mor_logHIndex)
#p-direction <0.90 -> drop from the model

#univariate model with logNInstitutions
fit_peto_mor_logNInstitutions<-  brm(NMortality | trials(NNecropsies)  ~ logHIndex + (1|gr(Species, cov =  vcv_brms_OU)),
                                     family =binomial(),
                                     data = data_mortality,
                                     data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                                     chains = 4,
                                     cores = 24,
                                     iter = 10000, warmup = 5000,
                                     thin = 5,
                                     seed = 1,
                                     control=list(max_treedepth = 12, adapt_delta = 0.99),
                                     save_pars = save_pars(all = TRUE))
fit_peto_mor_logNInstitutions
summary(fit_peto_mor_logNInstitutions,prob = 0.90)
p_direction(fit_peto_mor_logNInstitutions)
#p-direction <0.90 -> drop from the model

#univariate model with mean_GDP_5years
fit_peto_mor_mean_GDP_5years<-  brm(NMortality | trials(NNecropsies)  ~ mean_GDP_5years + (1|gr(Species, cov =  vcv_brms_OU)),
                                    family =binomial(),
                                    data = data_mortality,
                                    data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                                    chains = 4,
                                    cores = 24,
                                    iter = 10000, warmup = 5000,
                                    thin = 5,
                                    seed = 1,
                                    control=list(max_treedepth = 12, adapt_delta = 0.99),
                                    save_pars = save_pars(all = TRUE))
fit_peto_mor_mean_GDP_5years
summary(fit_peto_mor_mean_GDP_5years,prob = 0.90)
p_direction(fit_peto_mor_mean_GDP_5years)
#p-direction >0.90 -> keep in the model

#univariate model with ResidualslogNViews 
fit_peto_mor_ResidualslogNViews <-  brm(NMortality | trials(NNecropsies)  ~ ResidualslogNViews  + (1|gr(Species, cov =  vcv_brms_OU)),
                                        family =binomial(),
                                        data = data_mortality,
                                        data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                                        chains = 4,
                                        cores = 24,
                                        iter = 10000, warmup = 5000,
                                        thin = 5,
                                        seed = 1,
                                        control=list(max_treedepth = 12, adapt_delta = 0.99),
                                        save_pars = save_pars(all = TRUE))
fit_peto_mor_ResidualslogNViews 
summary(fit_peto_mor_ResidualslogNViews ,prob = 0.90)
p_direction(fit_peto_mor_ResidualslogNViews )
#p-direction >0.90 -> keep in the model

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Step 2: Drop part of the model selection ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#fir a model with all the variables selected in step 1
fit_peto_mor_m1 <-  brm(NMortality | trials(NNecropsies)  ~ logMass + logMaxLongevity + mean_GDP_5years + ResidualslogNViews  + (1|gr(Species, cov =  vcv_brms_OU)),
                                        family =binomial(),
                                        data = data_mortality,
                                        data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                                        chains = 4,
                                        cores = 24,
                                        iter = 10000, warmup = 5000,
                                        thin = 5,
                                        seed = 1,
                                        control=list(max_treedepth = 12, adapt_delta = 0.99),
                                        save_pars = save_pars(all = TRUE))
fit_peto_mor_m1
p_direction(fit_peto_mor_m1)
#mean_GDP_5years has the lowest p-direction, drop first

#model with GDP dropped
fit_peto_mor_m2 <-  brm(NMortality | trials(NNecropsies)  ~ logMass + logMaxLongevity + ResidualslogNViews  + (1|gr(Species, cov =  vcv_brms_OU)),
                        family =binomial(),
                        data = data_mortality,
                        data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                        chains = 4,
                        cores = 24,
                        iter = 10000, warmup = 5000,
                        thin = 5,
                        seed = 1,
                        control=list(max_treedepth = 12, adapt_delta = 0.99),
                        save_pars = save_pars(all = TRUE))
diagnostic(reference_model =fit_peto_mor_m1, adjusted_model = fit_peto_mor_m2)
fit_peto_mor_m2
p_direction(fit_peto_mor_m2)
#dropping the GDP from the model change the slope of other parameters <25% -> drop from the model
#drop ResidualslogNViews  next

#model with ResidualslogNViews dropped
fit_peto_mor_m3 <-  brm(NMortality | trials(NNecropsies)  ~ logMass + logMaxLongevity + (1|gr(Species, cov =  vcv_brms_OU)),
                        family =binomial(),
                        data = data_mortality,
                        data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                        chains = 4,
                        cores = 24,
                        iter = 10000, warmup = 5000,
                        thin = 5,
                        seed = 1,
                        control=list(max_treedepth = 12, adapt_delta = 0.99),
                        save_pars = save_pars(all = TRUE))
diagnostic(reference_model =fit_peto_mor_m2, adjusted_model = fit_peto_mor_m3)
fit_peto_mor_m3
p_direction(fit_peto_mor_m3)
#dropping ResidualslogNViews changes the slope logMass  and logMaxLongevity  >25% -> keep in the model


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Step 3: Final model ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#add back parameters dropped from step 1
#add back logNInstitutions to the model
fit_peto_mor_m3 <-  brm(NMortality | trials(NNecropsies)  ~ logMass + logMaxLongevity + logNInstitutions + ResidualslogNViews  + (1|gr(Species, cov =  vcv_brms_OU)),
                        family =binomial(),
                        data = data_mortality,
                        data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                        chains = 4,
                        cores = 24,
                        iter = 10000, warmup = 5000,
                        thin = 5,
                        seed = 1,
                        control=list(max_treedepth = 12, adapt_delta = 0.99),
                        save_pars = save_pars(all = TRUE))
diagnostic(reference_model =fit_peto_mor_m2, adjusted_model = fit_peto_mor_m3)
fit_peto_mor_m3
p_direction(fit_peto_mor_m3)
#logNInstitutions doesn't change the slopes by >25% -> exclude from model

#add back logHIndex to the model
fit_peto_mor_m4 <-  brm(NMortality | trials(NNecropsies)  ~ logMass + logMaxLongevity + logHIndex + ResidualslogNViews  + (1|gr(Species, cov =  vcv_brms_OU)),
                        family =binomial(),
                        data = data_mortality,
                        data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                        chains = 4,
                        cores = 24,
                        iter = 10000, warmup = 5000,
                        thin = 5,
                        seed = 1,
                        control=list(max_treedepth = 12, adapt_delta = 0.99),
                        save_pars = save_pars(all = TRUE))
diagnostic(reference_model =fit_peto_mor_m2, adjusted_model = fit_peto_mor_m4)
fit_peto_mor_m4
p_direction(fit_peto_mor_m4)
#logHIndex change the slope of  logMaxLongevity by >25% -> keep in the model

#check interaction effect between mass and lognevity
fit_peto_mor_m5 <-  brm(NMortality | trials(NNecropsies)  ~ logMass * logMaxLongevity + logHIndex + ResidualslogNViews  + (1|gr(Species, cov =  vcv_brms_OU)),
                        family =binomial(),
                        data = data_mortality,
                        data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                        chains = 4,
                        cores = 24,
                        iter = 10000, warmup = 5000,
                        thin = 5,
                        seed = 1,
                        control=list(max_treedepth = 12, adapt_delta = 0.99),
                        save_pars = save_pars(all = TRUE))
diagnostic(reference_model =fit_peto_mor_m4, adjusted_model = fit_peto_mor_m5)
fit_peto_mor_m5
p_direction(fit_peto_mor_m5)
#adding the interaction effect change the slope of logMass and logMaxLongevity -> keep it

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Calculate and plot the marginal effects ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

crude_model_mor <- brm(NMortality | trials(NNecropsies)  ~ logMaxLongevity * logMass + (1|gr(Species, cov =  vcv_brms_OU)),
                       family =binomial(),
                       data = data_mortality,
                       data2 = list(vcv_brms_OU = vcv_brms_OU_mortality), 
                       chains = 4,
                       cores = 24,
                       iter = 10000, warmup = 5000,
                       thin = 5,
                       seed = 1,
                       control=list(max_treedepth = 12, adapt_delta = 0.99),
                       save_pars = save_pars(all = TRUE))

#create the conditions for the crude and adjusted models
cdt_crude_peto_mor <-  list(NNecropsies = 100,
                            logMaxLongevity = median_logMaxLongevity)

cdt_adj_peto_mor <- list(logHIndex = median_logHindex, 
                         ResidualslogNViews = median_ResidualslogNViews,
                         logMaxLongevity = median_logMaxLongevity,  
                         NNecropsies = 100)

marginal_effects_peto_mal <- get_marginal_effects(crude_model = crude_model_mor, 
                                                  adjusted_model = fit_peto_mor_m5,
                                                  resp = "NMortality",
                                                  effect = "logMass",
                                                  varext = "logMass",
                                                  NNecropsies = 100, 
                                                  conditions_crude = cdt_crude_peto_mor, 
                                                  conditions_adjusted = cdt_adj_peto_mor)

#plot the marginal effects for mutation rate
df_peto_mal <- subset(marginal_effects_peto_mal, logMass >= min(data_mortality$logMass) & logMass <= max(data_mortality$logMass))
p_peto_mal <- ggplot(data = df_peto_mal , aes(x = logMass, y = estimate__, group = Type, color = Type, fill = Type))+
  geom_line(size = 1.5)+
  geom_ribbon(mapping = aes(x = logMass, ymin = lower__, ymax = upper__), color = NA, alpha = 0.25)+
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
  xlab("log10 body mass (g)")+
  ylab("Malignancy prevalence")+
  labs(color = "Model type:", fill = "Model type:")+
  ylim(0,0.5)
#xlim(-3.5,3.5)
p_peto_mal

