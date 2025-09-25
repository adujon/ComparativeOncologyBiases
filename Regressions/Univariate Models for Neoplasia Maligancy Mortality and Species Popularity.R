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

#setup import and export folders
export_dir <- "C:/Users/antoi/My Drive/Postdoc/Zoo data/Peto/models"
import_dir <- "C:/Users/antoi/My Drive/Postdoc/Zoo data/Biases/Data"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Import phylogeny data ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Neoplasia and malignancy  ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#import the phylogenetic tree for mammals
setwd("C:/Users/antoi/My Drive/Postdoc/Zoo data/Biases/Phylogenies")
phylo_all <- read.tree("Phylogeny_all_species.tre") 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import Neoplasia and malignancy data ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("C:/Users/antoi/My Drive/Postdoc/Zoo data/Biases/Scripts/Data exploration/For Github")

data_neo_mal <- read_xlsx("Raw dataset confounding variables.xlsx", sheet = "Neoplasia_malignancy", col_types = c("text", "text", rep("numeric", 21)))
data_neo_mal <- data.frame(data_neo_mal)
rownames(data_neo_mal) <- data_neo_mal$Species
data_neo_mal$NeoplasiaPrevalence <- data_neo_mal$NNeoplasia/data_neo_mal$NNecropsies
data_neo_mal$MalignancyPrevalence <- data_neo_mal$NMalignant/data_neo_mal$NNecropsies

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import Mortality data----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

data_mortality <- read_xlsx("Raw dataset confounding variables.xlsx", sheet = "Mortality", col_types = c("text", rep("numeric", 20)))
data_mortality <- data.frame(data_mortality)
rownames(data_mortality) <- data_mortality$Species
data_mortality$MortalityPrevalence <- data_mortality$NMortality / data_mortality$NNecropsies

#define colors for the plots
my_colors <- c("#E69F00", "#56B4E9", "limegreen", "purple")

#create an empty dataframe to store the effect size for each confounding variables
effect_sizes <- data.frame(OR  = NULL, Q2.5 = NULL, Q97.5 = NULL, Variable = NULL, Tumour = NULL)

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

#compute the var-covar matrix from all species
vcv_brms_OU_all <- vcv(corMartins(1, phylo_all , fixed=FALSE,form=~phylo_all$tip.label)) #OU model

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Tumor prevalence ~ mean_GDP_5years ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Neoplasia and malignancy ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#fit multivariate model predicting neoplsia and malignangy prevalence in function of mean_GDP_5years
fit_mean_GDP_5years_neo_mal <- brm(mvbind(NNeoplasia, NMalignant) | trials(NNecropsies)  ~ mean_GDP_5years + (1|gr(Species, cov = vcv_brms_OU)),
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

fit_mean_GDP_5years_neo_mal

#calculate odd ratios  and ORs for neoplasia
fe <- format(round(exp(fixed.effects(fit_mean_GDP_5years_neo_mal)), digits = 1), nsmall = 1)
OR <- paste("OR = ", fe[3,1], " (95%CI:", fe[3,3], "-", fe[3,4], ")", sep = "")
rsq_lbl <- as.expression(bquote(~.(OR)))
es_neoplasia <-  data.frame(OR  = fe[3,1], Q2.5 = fe[3,3], Q97.5 = fe[3,4], Variable = "GDP per capita", Tumour = "Neoplasia")
effect_sizes <- rbind(effect_sizes, es_neoplasia)


#calculate marginal effect for neoplasia
cond_effects <- conditional_effects(fit_mean_GDP_5years_neo_mal , effect = "mean_GDP_5years", resp = "NNeoplasia")
cond_effects_plot <- data.frame(cond_effects[[1]])
cond_effects_plot$estimate__ <- cond_effects_plot$estimate__/cond_effects_plot$NNecropsies[1]
cond_effects_plot$lower__ <- cond_effects_plot$lower__/cond_effects_plot$NNecropsies[1]
cond_effects_plot$upper__ <- cond_effects_plot$upper__/cond_effects_plot$NNecropsies[1]

#plot marginal effect for neoplasia (24 missing values)
p_mean_GDP_5years_neoplasia <- ggplot()+
  geom_ribbon(data = cond_effects_plot, mapping = aes(x = mean_GDP_5years, ymin = lower__, ymax = upper__), alpha = 0.25)+
  geom_line(data = cond_effects_plot, mapping = aes(x = mean_GDP_5years, y = estimate__), size = 1.5)+
  geom_point(data = data_neo_mal, mapping = aes(x = mean_GDP_5years, y = NeoplasiaPrevalence, color = Class), size = 3.5, alpha = 0.75)+
  theme_classic()+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=20,face="bold"))+
  theme(legend.position = c(0.2, 0.5))+
  theme(plot.title = element_text(size = 20, face = "bold"))+
  theme(legend.text=element_text(size=25),
        legend.title = element_text(size = 25))+
  guides(fill=guide_legend(title="Taxa: "))+
  scale_color_manual(values = my_colors) +
  ggtitle("Weighted average GDP per capita")+
  xlab("Weighted average GDP per capita\n (2018-2023, 10000's of USD)")+
  ylab("Neoplasia prevalence")+
  ylim(0,0.70)+
  xlim(1.5,7.5)+
  annotate(geom="text", x=4, y=0.68, label=rsq_lbl, size = 8)
p_mean_GDP_5years_neoplasia

#calculate odd ratios for malignancy
fe <- format(round(exp(fixed.effects(fit_mean_GDP_5years_neo_mal)), digits = 1), nsmall = 1)
OR <- paste("OR = ", fe[4,1], " (95%CI:", fe[4,3], "-", fe[4,4], ")", sep = "")
rsq_lbl <- as.expression(bquote(~.(OR)))
es_malignant <-  data.frame(OR  = fe[4,1], Q2.5 = fe[4,3], Q97.5 = fe[4,4], Variable = "GDP per capita", Tumour = "Malignancy")
effect_sizes <- rbind(effect_sizes, es_malignant)

#calculate marginal effect for malignancy
cond_effects <- conditional_effects(fit_mean_GDP_5years_neo_mal , "mean_GDP_5years", rep = "NMalignant")
cond_effects_plot <- data.frame(cond_effects[[1]])
cond_effects_plot$estimate__ <- cond_effects_plot$estimate__/cond_effects_plot$NNecropsies[1]
cond_effects_plot$lower__ <- cond_effects_plot$lower__/cond_effects_plot$NNecropsies[1]
cond_effects_plot$upper__ <- cond_effects_plot$upper__/cond_effects_plot$NNecropsies[1]

#plot the marginal effect for malignancy (24 missing values)
p_mean_GDP_5years_malignant <- ggplot()+
  geom_ribbon(data = cond_effects_plot, mapping = aes(x = mean_GDP_5years, ymin = lower__, ymax = upper__), alpha = 0.25)+
  geom_line(data = cond_effects_plot, mapping = aes(x = mean_GDP_5years, y = estimate__), size = 1.5)+
  geom_point(data = data_neo_mal, mapping = aes(x = mean_GDP_5years, y = MalignancyPrevalence, color = Class), size = 2)+
  theme_classic()+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=20,face="bold"))+
  theme(legend.position = c(0.2, 0.5))+
  theme(plot.title = element_text(size = 20, face = "bold"))+
  theme(legend.text=element_text(size=25),
        legend.title = element_text(size = 25))+
  guides(fill=guide_legend(title="Taxa: "))+
  scale_color_manual(values = my_colors) +
  ggtitle("Weighted average GDP per capita")+
  xlab("Weighted average GDP per capita\n (2018-2023, 10000's of USD)")+
  ylab("Malignancy prevalence")+
  ylim(0,0.70)+
  xlim(1.5,7.5)+
  annotate(geom="text", x=4, y=0.68, label=rsq_lbl, size = 8)
p_mean_GDP_5years_malignant

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Mortality ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#fit multivariate model predicting mortality prevalence in function of mean_GDP_5years
fit_mean_GDP_5years_mortality <- brm(NMortality | trials(NNecropsies)  ~ mean_GDP_5years + (1|gr(Species, cov = vcv_brms_OU)),
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

fit_mean_GDP_5years_mortality

#calculate odd ratios  for mortality
fe <- format(round(exp(fixed.effects(fit_mean_GDP_5years_mortality)), digits = 1), nsmall = 1)
OR <- paste("OR = ", fe[2,1], " (95%CI:", fe[2,3], "-", fe[2,4], ")", sep = "")
rsq_lbl <- as.expression(bquote(~.(OR)))
es_mortality <-  data.frame(OR  = fe[2,1], Q2.5 = fe[2,3], Q97.5 = fe[2,4], Variable = "GDP per capita", Tumour = "Mortality")
effect_sizes <- rbind(effect_sizes, es_mortality)

#calculate marginal effect for mortality
cond_effects <- conditional_effects(fit_mean_GDP_5years_mortality , effect = "mean_GDP_5years", resp = "NMalignant")
cond_effects_plot <- data.frame(cond_effects[[1]])
cond_effects_plot$estimate__ <- cond_effects_plot$estimate__/cond_effects_plot$NNecropsies[1]
cond_effects_plot$lower__ <- cond_effects_plot$lower__/cond_effects_plot$NNecropsies[1]
cond_effects_plot$upper__ <- cond_effects_plot$upper__/cond_effects_plot$NNecropsies[1]

#plot marginal effect for mortality (3 missing values)
p_mean_GDP_5years_mortality <- ggplot()+
  geom_ribbon(data = cond_effects_plot, mapping = aes(x = mean_GDP_5years, ymin = lower__, ymax = upper__), alpha = 0.25)+
  geom_line(data = cond_effects_plot, mapping = aes(x = mean_GDP_5years, y = estimate__), size = 1.5)+
  geom_point(data =data_mortality, mapping = aes(x = mean_GDP_5years, y = MortalityPrevalence), size = 3.5, alpha = 0.75, color = "limegreen")+
  theme_classic()+
  xlab("log10 number of Wikipedia views")+
  ylab("Neoplasia prevalence")+
  theme_classic()+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=20,face="bold"))+
  theme(legend.position = c(0.2, 0.5))+
  theme(plot.title = element_text(size = 20, face = "bold"))+
  theme(legend.text=element_text(size=12),
        legend.title = element_text(size = 12))+
  ggtitle("Weighted average GDP per capita")+
  xlab("Weighted average GDP per capita\n (2018-2023, 10000's of USD)")+
  ylab("Mortality prevalence")+
  ylim(0,0.70)+
  xlim(1.5,7.5)+
  annotate(geom="text", x=4, y=0.68, label=rsq_lbl, size = 8)
p_mean_GDP_5years_mortality

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Tumor prevalence ~ logNViews ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Neoplasia and malignancy ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#fit multivariate model predicting neoplsia and malignangy prevalence in function of logNViews
fit_logNViews_neo_mal <- brm(mvbind(NNeoplasia, NMalignant) | trials(NNecropsies)  ~ logNViews + (1|gr(Species, cov = vcv_brms_OU)),
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

fit_logNViews_neo_mal

#calculate odd ratios for neoplasia
fe <- format(round(exp(fixed.effects(fit_logNViews_neo_mal)), digits = 1), nsmall = 1)
OR <- paste("OR = ", fe[3,1], " (95%CI:", fe[3,3], "-", fe[3,4], ")", sep = "")
rsq_lbl <- as.expression(bquote(~.(OR)))
es_neoplasia <-  data.frame(OR  = fe[3,1], Q2.5 = fe[3,3], Q97.5 = fe[3,4], Variable = "log10 number of Wikipedia views", Tumour = "Neoplasia")
effect_sizes <- rbind(effect_sizes, es_neoplasia)

#calculate marginal effect for neoplasia
cond_effects <- conditional_effects(fit_logNViews_neo_mal , effect = "logNViews", resp = "NNeoplasia")
cond_effects_plot <- data.frame(cond_effects[[1]])
cond_effects_plot$estimate__ <- cond_effects_plot$estimate__/cond_effects_plot$NNecropsies[1]
cond_effects_plot$lower__ <- cond_effects_plot$lower__/cond_effects_plot$NNecropsies[1]
cond_effects_plot$upper__ <- cond_effects_plot$upper__/cond_effects_plot$NNecropsies[1]

#plot marginal effect for neoplasia
p_logNViews_neoplasia <- ggplot()+
  geom_ribbon(data = cond_effects_plot, mapping = aes(x = logNViews, ymin = lower__, ymax = upper__), alpha = 0.25)+
  geom_line(data = cond_effects_plot, mapping = aes(x = logNViews, y = estimate__), size = 1.5)+
  geom_point(data = data_neo_mal, mapping = aes(x = logNViews, y = NeoplasiaPrevalence, color = Class), size = 3.5, alpha = 0.75)+
  theme_classic()+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=20,face="bold"))+
  theme(legend.position = c(0.2, 0.5))+
  theme(legend.text=element_text(size=25),
        legend.title = element_text(size = 25))+
  theme(plot.title = element_text(size = 20, face = "bold"))+
  ggtitle("Wikipedia views")+
  xlab("log10 number of Wikipedia views")+
  ylab("Neoplasia prevalence")+
  scale_color_manual(values = my_colors) +
  ylim(0,0.70)+
  xlim(3,8)+
  annotate(geom="text", x=5, y=0.68, label=rsq_lbl, size = 8)
p_logNViews_neoplasia

#calculate odd ratios for malignancy
rsq <- rstantools::bayes_R2(fit_logNViews_neo_mal, re.form = NA)
fe <- format(round(exp(fixed.effects(fit_logNViews_neo_mal)), digits = 1), nsmall = 1)
OR <- paste("OR = ", fe[4,1], " (95%CI:", fe[4,3], "-", fe[4,4], ")", sep = "")
rsq_lbl <- as.expression(bquote(~.(OR)))
es_malignant<-  data.frame(OR  = fe[4,1], Q2.5 = fe[4,3], Q97.5 = fe[4,4], Variable = "log10 Number of views", Tumour = "Malignancy")
effect_sizes <- rbind(effect_sizes, es_malignant)

#calculate marginal effect for malignancy
cond_effects <- conditional_effects(fit_logNViews_neo_mal , "logNViews", rep = "NMalignant")
cond_effects_plot <- data.frame(cond_effects[[1]])
cond_effects_plot$estimate__ <- cond_effects_plot$estimate__/cond_effects_plot$NNecropsies[1]
cond_effects_plot$lower__ <- cond_effects_plot$lower__/cond_effects_plot$NNecropsies[1]
cond_effects_plot$upper__ <- cond_effects_plot$upper__/cond_effects_plot$NNecropsies[1]

#plot the marginal effect for malignancy
p_logNViews_malignant <- ggplot()+
  geom_ribbon(data = cond_effects_plot, mapping = aes(x = logNViews, ymin = lower__, ymax = upper__), alpha = 0.25)+
  geom_line(data = cond_effects_plot, mapping = aes(x = logNViews, y = estimate__), size = 1.5)+
  geom_point(data = data_neo_mal, mapping = aes(x = logNViews, y = MalignancyPrevalence, color = Class), size = 3, alpha = 0.75)+
  theme_classic()+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=20,face="bold"))+
  theme(legend.position = c(0.2, 0.5))+
  theme(plot.title = element_text(size = 20, face = "bold"))+
  theme(legend.text=element_text(size=12),
        legend.title = element_text(size = 12))+
  ggtitle("Wikipedia views")+
  xlab("log10 number of Wikipedia views")+
  ylab("Malignancy prevalence")+
  scale_color_manual(values = my_colors) +
  ylim(0,0.70)+
  xlim(3,8)+
  annotate(geom="text", x=5, y=0.68, label=rsq_lbl, size = 8)
p_logNViews_malignant

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Mortality ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#fit multivariate model predicting mortality prevalence prevalence in function of logNViews
fit_logNViews_mortality <- brm(NMortality | trials(NNecropsies)  ~ logNViews + (1|gr(Species, cov = vcv_brms_OU)),
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

fit_logNViews_mortality

#calculate odd ratios for mortality
rsq <- rstantools::bayes_R2(fit_logNViews_mortality, re.form = NA)
fe <- format(round(exp(fixed.effects(fit_logNViews_mortality)), digits = 1), nsmall = 1)
OR <- paste("OR = ", fe[2,1], " (95%CI:", fe[2,3], "-", fe[2,4], ")", sep = "")
rsq_lbl <- as.expression(bquote(~.(OR)))
es_mortality<-  data.frame(OR  = fe[2,1], Q2.5 = fe[2,3], Q97.5 = fe[2,4], Variable = "log10 Number of views", Tumour = "Mortality")
effect_sizes <- rbind(effect_sizes, es_mortality)

#calculate marginal effect for mortality
cond_effects <- conditional_effects(fit_logNViews_mortality , effect = "logNViews", resp = "NMalignant")
cond_effects_plot <- data.frame(cond_effects[[1]])
cond_effects_plot$estimate__ <- cond_effects_plot$estimate__/cond_effects_plot$NNecropsies[1]
cond_effects_plot$lower__ <- cond_effects_plot$lower__/cond_effects_plot$NNecropsies[1]
cond_effects_plot$upper__ <- cond_effects_plot$upper__/cond_effects_plot$NNecropsies[1]

#plot marginal effect for mortality
p_logNViews_mortality <- ggplot()+
  geom_ribbon(data = cond_effects_plot, mapping = aes(x = logNViews, ymin = lower__, ymax = upper__), alpha = 0.25)+
  geom_line(data = cond_effects_plot, mapping = aes(x = logNViews, y = estimate__), size = 1.5)+
  geom_point(data =data_mortality, mapping = aes(x = logNViews, y = MortalityPrevalence), size = 3.5, alpha = 0.75, color = "limegreen")+
  theme_classic()+
  xlab("log10 Number of views")+
  ylab("Neoplasia prevalence")+
  theme_classic()+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=20,face="bold"))+
  theme(legend.position = c(0.2, 0.5))+
  theme(plot.title = element_text(size = 20, face = "bold"))+
  theme(legend.text=element_text(size=12),
        legend.title = element_text(size = 12))+
  ggtitle("Wikipedia views")+
  xlab("log10 number of Wikipedia views")+
  ylab("Mortality prevalence")+
  ylim(0,0.70)+
  xlim(3,7)+
  annotate(geom="text", x=5, y=0.68, label=rsq_lbl, size = 8)
p_logNViews_mortality

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Tumor prevalence ~ ResidualslogNViews ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Neoplasia and malignancy ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#fit multivariate model predicting neoplsia and malignangy prevalence in function of ResidualslogNViews
fit_ResidualslogNViews_neo_mal <- brm(mvbind(NNeoplasia, NMalignant) | trials(NNecropsies)  ~ ResidualslogNViews + (1|gr(Species, cov = vcv_brms_OU)),
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

fit_ResidualslogNViews_neo_mal

#calculate odd ratios for neoplasia
fe <- format(round(exp(fixed.effects(fit_ResidualslogNViews_neo_mal)), digits = 1), nsmall = 1)
OR <- paste("OR = ", fe[3,1], " (95%CI:", fe[3,3], "-", fe[3,4], ")", sep = "")
rsq_lbl <- as.expression(bquote(~.(OR)))
es_neoplasia<-  data.frame(OR  = fe[3,1], Q2.5 = fe[3,3], Q97.5 = fe[3,4], Variable = "log10 Residuals number of views", Tumour = "Neoplasia")
effect_sizes <- rbind(effect_sizes, es_neoplasia)

#calculate marginal effect for neoplasia
cond_effects <- conditional_effects(fit_ResidualslogNViews_neo_mal , effect = "ResidualslogNViews", resp = "NNeoplasia")
cond_effects_plot <- data.frame(cond_effects[[1]])
cond_effects_plot$estimate__ <- cond_effects_plot$estimate__/cond_effects_plot$NNecropsies[1]
cond_effects_plot$lower__ <- cond_effects_plot$lower__/cond_effects_plot$NNecropsies[1]
cond_effects_plot$upper__ <- cond_effects_plot$upper__/cond_effects_plot$NNecropsies[1]

#plot marginal effect for neoplasia
p_ResidualslogNViews_neoplasia <- ggplot()+
  geom_ribbon(data = cond_effects_plot, mapping = aes(x = ResidualslogNViews, ymin = lower__, ymax = upper__), alpha = 0.25)+
  geom_line(data = cond_effects_plot, mapping = aes(x = ResidualslogNViews, y = estimate__), size = 1.5)+
  geom_point(data = data_neo_mal, mapping = aes(x = ResidualslogNViews, y = NeoplasiaPrevalence, color = Class), size = 3.5, aplha = 0.75)+
  theme_classic()+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=20,face="bold"))+
  theme(legend.position = c(0.2, 0.5))+
  theme(plot.title = element_text(size = 20, face = "bold"))+
  theme(legend.text=element_text(size=12),
        legend.title = element_text(size = 12))+
  ggtitle("Wikipedia views (residuals)")+
  xlab("log10 number of Wikipedia views (Residuals)")+
  ylab("Neoplasia prevalence")+
  scale_color_manual(values = my_colors) +
  ylim(0,0.70)+
  xlim(-2.5,2)+
  annotate(geom="text", x= 0, y=0.68, label=rsq_lbl, size = 8)
p_ResidualslogNViews_neoplasia

#calculate odd ratios  for malignancy
fe <- format(round(exp(fixed.effects(fit_ResidualslogNViews_neo_mal)), digits = 1), nsmall = 1)
OR <- paste("OR = ", fe[4,1], " (95%CI:", fe[4,3], "-", fe[4,4], ")", sep = "")
rsq_lbl <- as.expression(bquote(~.(OR)))
es_malignant<-  data.frame(OR  = fe[4,1], Q2.5 = fe[4,3], Q97.5 = fe[4,4], Variable = "log10 Residuals number of views", Tumour = "Malignancy")
effect_sizes <- rbind(effect_sizes, es_malignant)

#calculate marginal effect for malignancy
cond_effects <- conditional_effects(fit_ResidualslogNViews_neo_mal , "ResidualslogNViews", rep = "NMalignant")
cond_effects_plot <- data.frame(cond_effects[[1]])
cond_effects_plot$estimate__ <- cond_effects_plot$estimate__/cond_effects_plot$NNecropsies[1]
cond_effects_plot$lower__ <- cond_effects_plot$lower__/cond_effects_plot$NNecropsies[1]
cond_effects_plot$upper__ <- cond_effects_plot$upper__/cond_effects_plot$NNecropsies[1]

#plot the marginal effect for malignancy
p_ResidualslogNViews_malignant <- ggplot()+
  geom_ribbon(data = cond_effects_plot, mapping = aes(x = ResidualslogNViews, ymin = lower__, ymax = upper__), alpha = 0.25)+
  geom_line(data = cond_effects_plot, mapping = aes(x = ResidualslogNViews, y = estimate__), size = 1.5)+
  geom_point(data = data_neo_mal, mapping = aes(x = ResidualslogNViews, y = MalignancyPrevalence, color = Class), size = 2)+
  theme_classic()+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=20,face="bold"))+
  theme(legend.position = c(0.2, 0.5))+
  theme(plot.title = element_text(size = 20, face = "bold"))+
  theme(legend.text=element_text(size=12),
        legend.title = element_text(size = 12))+
  ggtitle("Wikipedia views (residuals)")+
  xlab("log10 number of Wikipedia views (Residuals)")+
  ylab("Malignancy prevalence")+
  scale_color_manual(values = my_colors) +
  ylim(0,0.70)+
  xlim(-2.5,2)+
  annotate(geom="text", x= 0, y=0.68, label=rsq_lbl, size = 8)
p_ResidualslogNViews_malignant

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Mortality ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#fit multivariate model predicting mortality prevalence in function of ResidualslogNViews
fit_ResidualslogNViews_mortality <- brm(NMortality | trials(NNecropsies)  ~ ResidualslogNViews + (1|gr(Species, cov = vcv_brms_OU)),
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

fit_ResidualslogNViews_mortality

#calculate odd ratios  for mortality
fe <- format(round(exp(fixed.effects(fit_ResidualslogNViews_mortality)), digits = 1), nsmall = 1)
OR <- paste("OR = ", fe[2,1], " (95%CI:", fe[2,3], "-", fe[2,4], ")", sep = "")
rsq_lbl <- as.expression(bquote(~.(OR)))
es_mortality<-  data.frame(OR  = fe[2,1], Q2.5 = fe[2,3], Q97.5 = fe[2,4], Variable = "log10 Residuals number of views", Tumour = "Mortality")
effect_sizes <- rbind(effect_sizes, es_mortality)

#calculate marginal effect for mortality
cond_effects <- conditional_effects(fit_ResidualslogNViews_mortality , effect = "ResidualslogNViews", resp = "NMalignant")
cond_effects_plot <- data.frame(cond_effects[[1]])
cond_effects_plot$estimate__ <- cond_effects_plot$estimate__/cond_effects_plot$NNecropsies[1]
cond_effects_plot$lower__ <- cond_effects_plot$lower__/cond_effects_plot$NNecropsies[1]
cond_effects_plot$upper__ <- cond_effects_plot$upper__/cond_effects_plot$NNecropsies[1]

#plot marginal effect for mortality
p_ResidualslogNViews_mortality <- ggplot()+
  geom_ribbon(data = cond_effects_plot, mapping = aes(x = ResidualslogNViews, ymin = lower__, ymax = upper__), alpha = 0.25)+
  geom_line(data = cond_effects_plot, mapping = aes(x = ResidualslogNViews, y = estimate__), size = 1.5)+
  geom_point(data =data_mortality, mapping = aes(x = ResidualslogNViews, y = MortalityPrevalence), size = 3.5, alpha = 0.75, color = "limegreen")+
  theme_classic()+
  xlab("log10 N Views")+
  ylab("Neoplasia prevalence")+
  theme_classic()+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=20,face="bold"))+
  theme(legend.position = c(0.2, 0.5))+
  theme(plot.title = element_text(size = 20, face = "bold"))+
  theme(legend.text=element_text(size=12),
        legend.title = element_text(size = 12))+
  ggtitle("Wikipedia views (residuals)")+
  xlab("log10 number of Wikipedia views (Residuals)")+
  ylab("Mortality prevalence")+
  xlim(-2.5,2)+
  annotate(geom="text", x= 0, y=0.68, label=rsq_lbl, size = 8)
p_ResidualslogNViews_mortality

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Tumor prevalence ~ logHIndex ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Neoplasia and malignancy ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#fit multivariate model predicting neoplsia and malignangy prevalence in function of logHIndex
fit_logHIndex_neo_mal <- brm(mvbind(NNeoplasia, NMalignant) | trials(NNecropsies)  ~ logHIndex + (1|gr(Species, cov = vcv_brms_OU)),
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

fit_logHIndex_neo_mal

#calculate odd ratios  for neoplasia
fe <- format(round(exp(fixed.effects(fit_logHIndex_neo_mal)), digits = 1), nsmall = 1)
OR <- paste("OR = ", fe[3,1], " (95%CI:", fe[3,3], "-", fe[3,4], ")", sep = "")
rsq_lbl <- as.expression(bquote(~.(OR)))
es_neoplasia<-  data.frame(OR  = fe[3,1], Q2.5 = fe[3,3], Q97.5 = fe[3,4], Variable = "log10 H-Index + 1", Tumour = "Neoplasia")
effect_sizes <- rbind(effect_sizes, es_neoplasia)

#calculate marginal effect for neoplasia
cond_effects <- conditional_effects(fit_logHIndex_neo_mal , effect = "logHIndex", resp = "NNeoplasia")
cond_effects_plot <- data.frame(cond_effects[[1]])
cond_effects_plot$estimate__ <- cond_effects_plot$estimate__/cond_effects_plot$NNecropsies[1]
cond_effects_plot$lower__ <- cond_effects_plot$lower__/cond_effects_plot$NNecropsies[1]
cond_effects_plot$upper__ <- cond_effects_plot$upper__/cond_effects_plot$NNecropsies[1]

#plot marginal effect for neoplasia
p_logHIndex_neoplasia <- ggplot()+
  geom_ribbon(data = cond_effects_plot, mapping = aes(x = logHIndex, ymin = lower__, ymax = upper__), alpha = 0.25)+
  geom_line(data = cond_effects_plot, mapping = aes(x = logHIndex, y = estimate__), size = 1.5)+
  geom_point(data = data_neo_mal, mapping = aes(x = logHIndex, y = NeoplasiaPrevalence, color = Class), size = 3.5, alpha = 0.75)+
  theme_classic()+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=20,face="bold"))+
  theme(legend.position = c(0.2, 0.5))+
  theme(plot.title = element_text(size = 20, face = "bold"))+
  theme(legend.text=element_text(size=12),
        legend.title = element_text(size = 12))+
  ggtitle("H-index")+
  xlab("log10 H-index+1")+
  ylab("Neoplasia prevalence")+
  scale_color_manual(values = my_colors)+
  ylim(0,0.70)+
  xlim(0,3)+
  annotate(geom="text", x=1.4, y=0.68, label=rsq_lbl, size = 8)
p_logHIndex_neoplasia

#calculate odd ratios  for malignancy
fe <- format(round(exp(fixed.effects(fit_logHIndex_neo_mal)), digits = 1), nsmall = 1)
OR <- paste("OR = ", fe[4,1], " (95%CI:", fe[4,3], "-", fe[4,4], ")", sep = "")
rsq_lbl <- as.expression(bquote(~.(OR)))
es_malignant<-  data.frame(OR  = fe[4,1], Q2.5 = fe[4,3], Q97.5 = fe[4,4], Variable = "log10 H-Index + 1", Tumour = "Malignancy")
effect_sizes <- rbind(effect_sizes, es_malignant)

#calculate marginal effect for malignancy
cond_effects <- conditional_effects(fit_logHIndex_neo_mal , "logHIndex", rep = "NMalignant")
cond_effects_plot <- data.frame(cond_effects[[1]])
cond_effects_plot$estimate__ <- cond_effects_plot$estimate__/cond_effects_plot$NNecropsies[1]
cond_effects_plot$lower__ <- cond_effects_plot$lower__/cond_effects_plot$NNecropsies[1]
cond_effects_plot$upper__ <- cond_effects_plot$upper__/cond_effects_plot$NNecropsies[1]

#plot the marginal effect for malignancy
p_logHIndex_malignant <- ggplot()+
  geom_ribbon(data = cond_effects_plot, mapping = aes(x = logHIndex, ymin = lower__, ymax = upper__), alpha = 0.25)+
  geom_line(data = cond_effects_plot, mapping = aes(x = logHIndex, y = estimate__), size = 1.5)+
  geom_point(data = data_neo_mal, mapping = aes(x = logHIndex, y = MalignancyPrevalence, color = Class), size = 3.5, alpha = 0.75)+
  theme_classic()+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=20,face="bold"))+
  theme(legend.position = c(0.2, 0.5))+
  theme(plot.title = element_text(size = 20, face = "bold"))+
  theme(legend.text=element_text(size=12),
        legend.title = element_text(size = 12))+
  ggtitle("H-Index")+
  xlab("log10 H-index+1")+
  ylab("Malignancy prevalence")+
  scale_color_manual(values = my_colors)+
  ylim(0,0.70)+
  xlim(0,3)+
  annotate(geom="text", x=1.4, y=0.68, label=rsq_lbl, size = 8)
p_logHIndex_malignant

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Mortality ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#fit multivariate model predicting mortality prevalence prevalence in function of logHIndex
fit_logHIndex_mortality <- brm(NMortality | trials(NNecropsies)  ~ logHIndex + (1|gr(Species, cov = vcv_brms_OU)),
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

fit_logHIndex_mortality

#compute the R2 for mortality
fe <- format(round(exp(fixed.effects(fit_logHIndex_mortality)), digits = 1), nsmall = 1)
OR <- paste("OR = ", fe[2,1], " (95%CI:", fe[2,3], "-", fe[2,4], ")", sep = "")
rsq_lbl <- as.expression(bquote(~.(OR)))
es_mortality<-  data.frame(OR  = fe[2,1], Q2.5 = fe[2,3], Q97.5 = fe[2,4], Variable = "log10 H-Index + 1", Tumour = "Mortality")
effect_sizes <- rbind(effect_sizes, es_mortality)

#calculate marginal effect for mortality
cond_effects <- conditional_effects(fit_logHIndex_mortality , effect = "logHIndex", resp = "NMalignant")
cond_effects_plot <- data.frame(cond_effects[[1]])
cond_effects_plot$estimate__ <- cond_effects_plot$estimate__/cond_effects_plot$NNecropsies[1]
cond_effects_plot$lower__ <- cond_effects_plot$lower__/cond_effects_plot$NNecropsies[1]
cond_effects_plot$upper__ <- cond_effects_plot$upper__/cond_effects_plot$NNecropsies[1]

#plot marginal effect for mortality
p_logHIndex_mortality <- ggplot()+
  geom_ribbon(data = cond_effects_plot, mapping = aes(x = logHIndex, ymin = lower__, ymax = upper__), alpha = 0.25)+
  geom_line(data = cond_effects_plot, mapping = aes(x = logHIndex, y = estimate__), size = 1.5)+
  geom_point(data =data_mortality, mapping = aes(x = logHIndex, y = MortalityPrevalence), size = 3.5, alpha = 0.75, color = "limegreen")+
  theme_classic()+
  xlab("log10 H-index+1")+
  ylab("Neoplasia prevalence")+
  theme_classic()+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=20,face="bold"))+
  theme(legend.position = c(0.2, 0.5))+
  theme(plot.title = element_text(size = 20, face = "bold"))+
  theme(legend.text=element_text(size=12),
        legend.title = element_text(size = 12))+
  ggtitle("H-index")+
  xlab("log10 H-index+1")+
  ylab("Mortality prevalence")+
  ylim(0,0.70)+
  xlim(0,3)+
  annotate(geom="text", x=1.4, y=0.68, label=rsq_lbl, size = 8)
p_logHIndex_mortality


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Tumor prevalence ~ Number of institutions ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Neoplasia and malignancy ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#fit multivariate model predicting neoplasia and malignancy prevalence in function of logHIndex
fit_logNInstitutions_neo_mal <- brm(mvbind(NNeoplasia, NMalignant) | trials(NNecropsies)  ~ logNInstitutions + (1|gr(Species, cov = vcv_brms_OU)),
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

fit_logNInstitutions_neo_mal

#calculate odd ratios for neoplasia
fe <- format(round(exp(fixed.effects(fit_logNInstitutions_neo_mal)), digits = 1), nsmall = 1)
OR <- paste("OR = ", fe[3,1], " (95%CI:", fe[3,3], "-", fe[3,4], ")", sep = "")
rsq_lbl <- as.expression(bquote(~.(OR)))
es_neoplasia<-  data.frame(OR  = fe[3,1], Q2.5 = fe[3,3], Q97.5 = fe[3,4], Variable = "log10 Number of institutions", Tumour = "Neoplasia")
effect_sizes <- rbind(effect_sizes, es_neoplasia)

#calculate the marginal effect for neoplasia
cond_effects <- conditional_effects(fit_logNInstitutions_neo_mal , effect = "logNInstitutions", resp = "NNeoplasia")
cond_effects_plot <- data.frame(cond_effects[[1]])
cond_effects_plot$estimate__ <- cond_effects_plot$estimate__/cond_effects_plot$NNecropsies[1]
cond_effects_plot$lower__ <- cond_effects_plot$lower__/cond_effects_plot$NNecropsies[1]
cond_effects_plot$upper__ <- cond_effects_plot$upper__/cond_effects_plot$NNecropsies[1]

#plot marginal effect for neoplasia  (24 missing value)
p_logNInstitutions_neoplasia <- ggplot()+
  geom_ribbon(data = cond_effects_plot, mapping = aes(x = logNInstitutions, ymin = lower__, ymax = upper__), alpha = 0.25)+
  geom_line(data = cond_effects_plot, mapping = aes(x = logNInstitutions, y = estimate__), size = 1.5)+
  geom_point(data = data_neo_mal, mapping = aes(x = logNInstitutions, y = NeoplasiaPrevalence, color = Class), size = 3, alpha = 0.75)+
  theme_classic()+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=20,face="bold"))+
  theme(legend.position = c(0.2, 0.5))+
  theme(plot.title = element_text(size = 20, face = "bold"))+
  theme(legend.text=element_text(size=12),
        legend.title = element_text(size = 12))+
  ggtitle("Number of institutions")+
  xlab("log10 number of institutions")+
  ylab("Neoplasia prevalence")+
  scale_color_manual(values = my_colors)+
  ylim(0,0.70)+
  xlim(0,3)+
  annotate(geom="text", x=1.5, y=0.68, label=rsq_lbl, size = 8)
p_logNInstitutions_neoplasia


#calculate odd ratios for malignancy
fe <- format(round(exp(fixed.effects(fit_logNInstitutions_neo_mal)), digits = 1), nsmall = 1)
OR <- paste("OR = ", fe[4,1], " (95%CI:", fe[4,3], "-", fe[4,4], ")", sep = "")
rsq_lbl <- as.expression(bquote(~.(OR)))
es_malignant<-  data.frame(OR  = fe[4,1], Q2.5 = fe[4,3], Q97.5 = fe[4,4], Variable = "log10 Number of institutions", Tumour = "Malignancy")
effect_sizes <- rbind(effect_sizes, es_malignant)

#calculate the predicted values to plot them
cond_effects <- conditional_effects(fit_logNInstitutions_neo_mal , effect = "logNInstitutions", resp = "NMalignant")
cond_effects_plot <- data.frame(cond_effects[[1]])
cond_effects_plot$estimate__ <- cond_effects_plot$estimate__/cond_effects_plot$NNecropsies[1]
cond_effects_plot$lower__ <- cond_effects_plot$lower__/cond_effects_plot$NNecropsies[1]
cond_effects_plot$upper__ <- cond_effects_plot$upper__/cond_effects_plot$NNecropsies[1]

p_logNInstitutions_malignant <- ggplot()+
  geom_ribbon(data = cond_effects_plot, mapping = aes(x = logNInstitutions, ymin = lower__, ymax = upper__), alpha = 0.25)+
  geom_line(data = cond_effects_plot, mapping = aes(x = logNInstitutions, y = estimate__), size = 1.5)+
  geom_point(data = data_neo_mal, mapping = aes(x = logNInstitutions, y = MalignancyPrevalence, color = Class), size = 3.5, alpha = 0.75)+
  theme_classic()+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=20,face="bold"))+
  theme(legend.position = c(0.2, 0.5))+
  theme(plot.title = element_text(size = 20, face = "bold"))+
  theme(legend.text=element_text(size=12),
        legend.title = element_text(size = 12))+
  ggtitle("Number of institutions")+
  xlab("log10 number of institutions")+
  ylab("Malignancy prevalence")+
  scale_color_manual(values = my_colors)+
  ylim(0,0.70)+
  xlim(0,3)+
  annotate(geom="text", x=1.5, y=0.68, label=rsq_lbl, size = 8)
p_logNInstitutions_malignant

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Mortality ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#fit multivariate model predicting mortality prevalence in function of logHIndex
fit_logNInstitutions_mortality <- brm( NMortality | trials(NNecropsies)  ~ logNInstitutions + (1|gr(Species, cov = vcv_brms_OU)),
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

fit_logNInstitutions_mortality

#calculate odd ratios for mortality
fe <- format(round(exp(fixed.effects(fit_logNInstitutions_mortality)), digits = 1), nsmall = 1)
OR <- paste("OR = ", fe[2,1], " (95%CI:", fe[2,3], "-", fe[2,4], ")", sep = "")
rsq_lbl <- as.expression(bquote(~.(OR)))
es_mortality<-  data.frame(OR  = fe[2,1], Q2.5 = fe[2,3], Q97.5 = fe[2,4], Variable = "log10 Number of institutions", Tumour = "Mortality")
effect_sizes <- rbind(effect_sizes, es_mortality)

#calculate the marginal effect for mortality (3 missing values)
cond_effects <- conditional_effects(fit_logNInstitutions_mortality , effect = "logNInstitutions", resp = "NMalignant")
cond_effects_plot <- data.frame(cond_effects[[1]])

p_logNInstitutions_mortality <- ggplot()+
  geom_ribbon(data = cond_effects_plot, mapping = aes(x = logNInstitutions, ymin = lower__, ymax = upper__), alpha = 0.25)+
  geom_line(data = cond_effects_plot, mapping = aes(x = logNInstitutions, y = estimate__), size = 1.5)+
  geom_point(data = data_mortality, mapping = aes(x = logNInstitutions, y = MortalityPrevalence), color = "limegreen", size = 3.5, alpha = 0.75)+
  theme_classic()+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=20,face="bold"))+
  theme(legend.position = c(0.2, 0.5))+
  theme(plot.title = element_text(size = 20, face = "bold"))+
  theme(legend.text=element_text(size=12),
        legend.title = element_text(size = 12))+
  ggtitle("Number of institutions")+
  xlab("log10 number of institutions")+
  ylab("Mortality prevalence")+
  ylim(0,0.70)+
  xlim(0,3)+
  annotate(geom="text", x=1.5, y=0.68, label=rsq_lbl, size = 8)
p_logNInstitutions_mortality

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Tumor prevalence ~ Number of ExSitu individuals ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Neoplasia and malignancy ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#fit multivariate model predicting tumor prevalence in function of logNExSitu
fit_logNExSitu_neo_mal <- brm(mvbind(NNeoplasia, NMalignant) | trials(NNecropsies)  ~ logNExSitu + (1|gr(Species, cov = vcv_brms_OU)),
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

fit_logNExSitu_neo_mal

#calculate odd ratios for neoplasia
fe <- format(round(exp(fixed.effects(fit_logNExSitu_neo_mal)), digits = 1), nsmall = 1)
OR <- paste("OR = ", fe[3,1], " (95%CI:", fe[3,3], "-", fe[3,4], ")", sep = "")
rsq_lbl <- as.expression(bquote(~.(OR)))
es_neoplasia<-  data.frame(OR  = fe[3,1], Q2.5 = fe[3,3], Q97.5 = fe[3,4], Variable = "log10 Number of animals in captivity", Tumour = "Neoplasia")
effect_sizes <- rbind(effect_sizes, es_neoplasia)

#calculate the marginal effects for neoplasia
cond_effects <- conditional_effects(fit_logNExSitu_neo_mal , effect = "logNExSitu", resp = "NNeoplasia")
cond_effects_plot <- data.frame(cond_effects[[1]])
cond_effects_plot$estimate__ <- cond_effects_plot$estimate__/cond_effects_plot$NNecropsies[1]
cond_effects_plot$lower__ <- cond_effects_plot$lower__/cond_effects_plot$NNecropsies[1]
cond_effects_plot$upper__ <- cond_effects_plot$upper__/cond_effects_plot$NNecropsies[1]

#plot the marginal effects for neoplasia (24 missing values)
p_logNExSitu_neoplasia <- ggplot()+
  geom_ribbon(data = cond_effects_plot, mapping = aes(x = logNExSitu, ymin = lower__, ymax = upper__), alpha = 0.25)+
  geom_line(data = cond_effects_plot, mapping = aes(x = logNExSitu, y = estimate__), size = 1.5)+
  geom_point(data = data_neo_mal, mapping = aes(x = logNExSitu, y = NeoplasiaPrevalence, color = Class), size = 3.5, alpha = 0.75)+
  theme_classic()+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=20,face="bold"))+
  theme(legend.position = c(0.2, 0.5))+
  theme(plot.title = element_text(size = 20, face = "bold"))+
  theme(legend.text=element_text(size=12),
        legend.title = element_text(size = 12))+
  ggtitle("Number of animals in captivity")+
  xlab("log10 number of animals in captivity")+
  ylab("Neoplasia prevalence")+
  scale_color_manual(values = my_colors)+
  ylim(0,0.70)+
  xlim(0,4.5)+
  annotate(geom="text", x=2.1, y=0.68, label=rsq_lbl, size = 8)
p_logNExSitu_neoplasia


#calculate odd ratios for malignancy
fe <- format(round(exp(fixed.effects(fit_logNExSitu_neo_mal)), digits = 1), nsmall = 1)
OR <- paste("OR = ", fe[4,1], " (95%CI:", fe[4,3], "-", fe[4,4], ")", sep = "")
rsq_lbl <- as.expression(bquote(~.(OR)))
es_malignant<-  data.frame(OR  = fe[4,1], Q2.5 = fe[4,3], Q97.5 = fe[4,4], Variable = "log10 Number of animals in captivity", Tumour = "Malignancy")
effect_sizes <- rbind(effect_sizes, es_malignant)

#calculate the marginal effects for malignancy
cond_effects <- conditional_effects(fit_logNExSitu_neo_mal , effect = "logNExSitu", resp = "NMalignant")
cond_effects_plot <- data.frame(cond_effects[[1]])
cond_effects_plot$estimate__ <- cond_effects_plot$estimate__/cond_effects_plot$NNecropsies[1]
cond_effects_plot$lower__ <- cond_effects_plot$lower__/cond_effects_plot$NNecropsies[1]
cond_effects_plot$upper__ <- cond_effects_plot$upper__/cond_effects_plot$NNecropsies[1]

#plot the marginal effects for malignancy(24 missing values)
p_logNExSitu_malignant <- ggplot()+
  geom_ribbon(data = cond_effects_plot, mapping = aes(x = logNExSitu, ymin = lower__, ymax = upper__), alpha = 0.25)+
  geom_line(data = cond_effects_plot, mapping = aes(x = logNExSitu, y = estimate__), size = 1.5)+
  geom_point(data = data_neo_mal, mapping = aes(x = logNExSitu, y = MalignancyPrevalence, color = Class), size = 3.5, alpha = 0.75)+
  theme_classic()+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=20,face="bold"))+
  theme(legend.position = c(0.2, 0.5))+
  theme(plot.title = element_text(size = 20, face = "bold"))+
  theme(legend.text=element_text(size=12),
        legend.title = element_text(size = 12))+
  ggtitle("Number of animals in captivity")+
  xlab("log10 number of animals in captivity")+
  ylab("Maligancy prevalence")+
  scale_color_manual(values = my_colors)+
  ylim(0,0.70)+
  xlim(0,4.5)+
  annotate(geom="text", x=2.1, y=0.68, label=rsq_lbl, size = 8)
p_logNExSitu_malignant


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Mortality ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#fit multivariate model predicting mortality prevalence in function of logNExSitu
fit_logNExSitu_mortality <- brm( NMortality  | trials(NNecropsies)  ~ logNExSitu + (1|gr(Species, cov = vcv_brms_OU)),
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

fit_logNExSitu_mortality

#calculate odd ratios for mortality
fe <- format(round(exp(fixed.effects(fit_logNExSitu_mortality)), digits = 1), nsmall = 1)
OR <- paste("OR = ", fe[2,1], " (95%CI:", fe[2,3], "-", fe[2,4], ")", sep = "")
rsq_lbl <- as.expression(bquote(~.(OR)))
es_mortality<-  data.frame(OR  = fe[2,1], Q2.5 = fe[2,3], Q97.5 = fe[2,4], Variable = "log10 Number of animals in captivity", Tumour = "Mortality")
effect_sizes <- rbind(effect_sizes, es_mortality)

#calculate the marginal effects for mortality
cond_effects <- conditional_effects(fit_logNExSitu_mortality , effect = "logNExSitu", resp = "NMalignant")
cond_effects_plot <- data.frame(cond_effects[[1]])

#plot marginal effect for mortality (3 missing values)
p_logNExSitu_mortality <- ggplot()+
  geom_ribbon(data = cond_effects_plot, mapping = aes(x = logNExSitu, ymin = lower__, ymax = upper__), alpha = 0.25)+
  geom_line(data = cond_effects_plot, mapping = aes(x = logNExSitu, y = estimate__), size = 1.5)+
  geom_point(data = data_mortality, mapping = aes(x = logNExSitu, y = MortalityPrevalence), color = "limegreen", size = 3.5, alpha = 0.75)+
  theme_classic()+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=20,face="bold"))+
  theme(legend.position = c(0.2, 0.5))+
  theme(plot.title = element_text(size = 18, face = "bold"))+
  theme(legend.text=element_text(size=16),
        legend.title = element_text(size = 12))+
  ggtitle("Number of animals in captivity")+
  xlab("log10 number of animals in captivity")+
  ylab("Mortality prevalence")+
  ylim(0,0.70)+
  xlim(0,4.5)+
  annotate(geom="text", x=2.1, y=0.68, label=rsq_lbl, size = 8)
p_logNExSitu_mortality

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create a plot of effect sizes ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#change working directory
setwd("C:/Users/antoi/My Drive/Postdoc/Zoo data/Biases/Figures")

effect_sizes$Tumour <- factor(effect_sizes$Tumour, levels = c("Neoplasia", "Malignancy", "Mortality"))
effect_sizes$OR <- as.numeric(effect_sizes$OR)
effect_sizes$Q2.5 <- as.numeric(effect_sizes$Q2.5)
effect_sizes$Q97.5 <- as.numeric(effect_sizes$Q97.5)
effect_sizes$Variable[effect_sizes$Variable == "GDP per capita"] <- "Weighted average GDP per capita (10000's of USD)"
effect_sizes$Variable[effect_sizes$Variable == "log10 Number of views"] <- "log10 number of Wikipedia views"
effect_sizes$Variable[effect_sizes$Variable == "log10 H-Index + 1"] <- "log10 H-index + 1"
effect_sizes$Variable[effect_sizes$Variable == "log10 Residuals number of views"] <- "log10 residuals number of Wikipedia views"
effect_sizes$Variable[effect_sizes$Variable == "log10 Number of animals in captivity"] <- "log10 number of animals in captivity"
effect_sizes$Variable[effect_sizes$Variable == "log10 Number of institutions"] <- "log10 number of institutions"
  
fig <- ggplot(data = effect_sizes, aes(y = OR, x  = Variable))+
  geom_hline(yintercept=1, size = 1.5, color = "darkgrey")+
  geom_linerange(aes(ymin = Q2.5, ymax = Q97.5), color = "skyblue", size = 1.5)+
  geom_point(size = 4, color = "red")+
  facet_wrap(~Tumour)+
  theme_classic()+
  ylim(0,6)+
  theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust=1, face="bold"),
        panel.border = element_rect(colour = "black", fill = NA))+
  theme(axis.text=element_text(size=12),
         axis.title=element_text(size=16,face="bold"),
         strip.text = element_text(face="bold", size=13))+
  xlab("Confounding variable")+
  ylab("Odds ratio and 95% credible intervals")
fig

png("Effect_Size.png", width = 35, height = 25, unit = "cm",res = 300)
fig
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Combining the plots for supplementary material ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~
## Neoplasia ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~

pneoplasia <- ggarrange(
  
  p_mean_GDP_5years_neoplasia,
  p_logNViews_neoplasia,
  p_ResidualslogNViews_neoplasia,
  p_logNInstitutions_neoplasia,
  p_logNExSitu_neoplasia,
  p_logHIndex_neoplasia,
  common.legend = T,
  legend = "bottom",
  nrow = 3,
  ncol = 2,
  align = "h"
)
pneoplasia

png("Neoplasia.png", width = 40, height = 50, unit = "cm",res = 300)
pneoplasia
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~
## Malignancy ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~

pmalignancy<- ggarrange(
  
  p_mean_GDP_5years_malignant,
  p_logNViews_malignant,
  p_ResidualslogNViews_malignant,
  p_logNInstitutions_malignant,
  p_logNExSitu_malignant,
  p_logHIndex_malignant,
  common.legend = T,
  legend = "bottom",
  nrow = 3,
  ncol = 2,
  align = "h"
)
pmalignancy

png("Malignancy.png", width = 40, height = 50, unit = "cm",res = 300)
pmalignancy
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~
## Mortality ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~

pmortality<- ggarrange(
  
  p_mean_GDP_5years_mortality,
  p_logNViews_mortality,
  p_ResidualslogNViews_mortality,
  p_logNInstitutions_mortality,
  p_logNExSitu_mortality,
  p_logHIndex_mortality,
  common.legend = T,
  legend = "bottom",
  nrow = 3,
  ncol = 2,
  align = "h"
)
pmortality

png("Mortality.png", width = 40, height = 50, unit = "cm",res = 300)
pmortality
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create a plot for Figure 2 ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

p_logNViews_neoplasia_f <-  p_logNViews_neoplasia + ggtitle("(A) Wikipedia views - Neoplasia")
p_logHIndex_neoplasia_f <-  p_logHIndex_neoplasia + ggtitle("(B) H-index - Neoplasia")
p_logNViews_malignant_f <- p_logNViews_malignant +  ggtitle("(C) Wikipedia views - Malignancy")
p_logHIndex_malignant_f <- p_logHIndex_malignant + ggtitle("(D) H-index - Malignancy")
p_logNViews_mortality_f <- p_logNViews_mortality + ggtitle("(E) Wikipedia views - Mortality")
p_logHIndex_mortality_f <-   p_logHIndex_mortality + ggtitle("(F) H-index - Mortality")
pfigure <- ggarrange(
  
  p_logNViews_neoplasia_f,
  p_logHIndex_neoplasia_f,
  p_logNViews_malignant_f,
  p_logHIndex_malignant_f,
  p_logNViews_mortality_f,
  p_logHIndex_mortality_f,
  common.legend = T,
  legend = "bottom",
  nrow = 3,
  ncol = 2,
  align = "h"
)
pfigure


png("Figure 2.png", width = 40, height = 50, unit = "cm",res = 300)
pfigure
dev.off()