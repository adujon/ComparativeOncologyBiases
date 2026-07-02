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
library(ggtext)

#setup import and export folders
export_dir <- "C:/Users/antoi/My Drive/Postdoc/Zoo data/Peto/models"
import_dir <- "C:/Users/antoi/My Drive/Postdoc/Zoo data/Biases/Data"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Import phylogeny data ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#import the phylogenetic tree for mammals
setwd("C:/Users/antoi/My Drive/Postdoc/Zoo data/Biases/Phylogenies")
phylo_all <- read.tree("Phylogeny_all_species.tre") 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import Mortality data----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("C:/Users/antoi/My Drive/Postdoc/Zoo data/Biases/reanalysis pageview proxy")
data_mortality <- read_xlsx("Raw dataset confounding variables.xlsx", sheet = "Mortality_CMR", col_types = c("text", rep("numeric", 20), "text", "numeric", "numeric"))
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
obj_brms_mortality<- name.check(phylo_all, data_mortality)
obj_brms_mortality

#drop the species from the tree for which we have no cancer data
phylo_mortality <- drop.tip(phylo_all, obj_brms_mortality$tree_not_data)

#compute the var-covar matrix from the tree to use in brms
vcv_brms_OU_mortality <- vcv(corMartins(1, phylo_mortality , fixed=FALSE,form=~phylo_mortality$tip.label)) #OU model

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Tumor prevalence ~ mean_GDP_5years ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
OR <- paste("OR = ", fe[2,1], " (95% CI:", fe[2,3], ", ", fe[2,4], ")", sep = "")
rsq_lbl <- as.expression(bquote(~.(OR)))
es_mortality <-  data.frame(OR  = fe[2,1], Q2.5 = fe[2,3], Q97.5 = fe[2,4], Variable = "GDP per capita", Tumour = "Mortality")
effect_sizes <- rbind(effect_sizes, es_mortality)

#calculate marginal effect for mortality
cond_effects <- conditional_effects(fit_mean_GDP_5years_mortality , effect = "mean_GDP_5years", resp = "NMortality")
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
  ggtitle(expression("Weighted average GDP " * italic("per capita"))) +
  xlab(expression(atop("Weighted average GDP " * italic("per capita"),
                       "(2018-2023, 10 000's of USD)"))) +
  ylab("Mortality prevalence (CMR based)")+
  ylim(0,0.70)+
  xlim(1.5,7.5)+
  annotate(geom="text", x=4, y=0.68, label=rsq_lbl, size = 8)
p_mean_GDP_5years_mortality

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Tumor prevalence ~ logNViews ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Mortality - Wikimedia data ---- 
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
OR <- paste("OR = ", fe[2,1], " (95% CI:", fe[2,3], ", ", fe[2,4], ")", sep = "")
rsq_lbl <- as.expression(bquote(~.(OR)))
es_mortality<-  data.frame(OR  = fe[2,1], Q2.5 = fe[2,3], Q97.5 = fe[2,4], Variable = "log10 Number of views", Tumour = "Mortality")
effect_sizes <- rbind(effect_sizes, es_mortality)

#calculate marginal effect for mortality
cond_effects <- conditional_effects(fit_logNViews_mortality , effect = "logNViews", resp = "NMortality")
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
  ylab("Mortality prevalence (CMR based)")+
  ylim(0,0.70)+
  xlim(3,7)+
  annotate(geom="text", x=5, y=0.68, label=rsq_lbl, size = 8)
p_logNViews_mortality

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Tumor prevalence ~ ResidualslogNViews ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
OR <- paste("OR = ", fe[2,1], " (95% CI :", fe[2,3], ", ", fe[2,4], ")", sep = "")
rsq_lbl <- as.expression(bquote(~.(OR)))
es_mortality<-  data.frame(OR  = fe[2,1], Q2.5 = fe[2,3], Q97.5 = fe[2,4], Variable = "log10 Residuals number of views", Tumour = "Mortality")
effect_sizes <- rbind(effect_sizes, es_mortality)

#calculate marginal effect for mortality
cond_effects <- conditional_effects(fit_ResidualslogNViews_mortality , effect = "ResidualslogNViews", resp = "NMortality")
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
  ylab("Mortality prevalence (CMR based)")+
  xlim(-2.5,2)+
  annotate(geom="text", x= 0, y=0.68, label=rsq_lbl, size = 8)
p_ResidualslogNViews_mortality

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Tumor prevalence ~ logHIndex ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
OR <- paste("OR = ", fe[2,1], " (95% CI : ", fe[2,3], ", ", fe[2,4], ")", sep = "")
rsq_lbl <- as.expression(bquote(~.(OR)))
es_mortality<-  data.frame(OR  = fe[2,1], Q2.5 = fe[2,3], Q97.5 = fe[2,4], Variable = "log10 H-Index + 1", Tumour = "Mortality")
effect_sizes <- rbind(effect_sizes, es_mortality)

#calculate marginal effect for mortality
cond_effects <- conditional_effects(fit_logHIndex_mortality , effect = "logHIndex", resp = "NMortality")
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
  ylab("Mortality prevalence (CMR based)")+
  ylim(0,0.70)+
  xlim(0,3)+
  annotate(geom="text", x=1.4, y=0.68, label=rsq_lbl, size = 8)
p_logHIndex_mortality

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Tumor prevalence ~ Number of institutions ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
OR <- paste("OR = ", fe[2,1], " (95% CI:", fe[2,3], ", ", fe[2,4], ")", sep = "")
rsq_lbl <- as.expression(bquote(~.(OR)))
es_mortality<-  data.frame(OR  = fe[2,1], Q2.5 = fe[2,3], Q97.5 = fe[2,4], Variable = "log10 Number of institutions", Tumour = "Mortality")
effect_sizes <- rbind(effect_sizes, es_mortality)

#calculate marginal effect for mortality
cond_effects <- conditional_effects(fit_logNInstitutions_mortality , effect = "logNInstitutions", resp = "NMortality")
cond_effects_plot <- data.frame(cond_effects[[1]])
cond_effects_plot$estimate__ <- cond_effects_plot$estimate__/cond_effects_plot$NNecropsies[1]
cond_effects_plot$lower__ <- cond_effects_plot$lower__/cond_effects_plot$NNecropsies[1]
cond_effects_plot$upper__ <- cond_effects_plot$upper__/cond_effects_plot$NNecropsies[1]

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
  ylab("Mortality prevalence (CMR based)")+
  ylim(0,0.70)+
  xlim(0,3)+
  annotate(geom="text", x=1.5, y=0.68, label=rsq_lbl, size = 8)
p_logNInstitutions_mortality

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Tumor prevalence ~ Number of ExSitu individuals ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
OR <- paste("OR = ", fe[2,1], " (95% CI:", fe[2,3], ", ", fe[2,4], ")", sep = "")
rsq_lbl <- as.expression(bquote(~.(OR)))
es_mortality<-  data.frame(OR  = fe[2,1], Q2.5 = fe[2,3], Q97.5 = fe[2,4], Variable = "log10 Number of animals in captivity", Tumour = "Mortality")
effect_sizes <- rbind(effect_sizes, es_mortality)


#calculate marginal effect for mortality
cond_effects <- conditional_effects(fit_logNExSitu_mortality , effect = "logNExSitu", resp = "NMortality")
cond_effects_plot <- data.frame(cond_effects[[1]])
cond_effects_plot$estimate__ <- cond_effects_plot$estimate__/cond_effects_plot$NNecropsies[1]
cond_effects_plot$lower__ <- cond_effects_plot$lower__/cond_effects_plot$NNecropsies[1]
cond_effects_plot$upper__ <- cond_effects_plot$upper__/cond_effects_plot$NNecropsies[1]

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
  ylab("Mortality prevalence (CMR based)")+
  ylim(0,0.70)+
  xlim(0,4.5)+
  annotate(geom="text", x=2.1, y=0.68, label=rsq_lbl, size = 8)
p_logNExSitu_mortality

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create a plot of effect sizes ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#change working directory
#setwd("C:/Users/antoi/My Drive/Postdoc/Zoo data/Biases/Figures")

effect_sizes$Tumour <- factor(effect_sizes$Tumour, levels = c("Neoplasia", "Malignancy", "Mortality"))
effect_sizes$OR <- as.numeric(effect_sizes$OR)
effect_sizes$Q2.5 <- as.numeric(effect_sizes$Q2.5)
effect_sizes$Q97.5 <- as.numeric(effect_sizes$Q97.5)
effect_sizes$Variable[effect_sizes$Variable == "GDP per capita"] <- "Weighted average GDP <i>per capita</i> (10 000's of USD)"
effect_sizes$Variable[effect_sizes$Variable == "log10 Number of views"] <- "log10 number of Wikipedia views"
effect_sizes$Variable[effect_sizes$Variable == "log10 H-Index + 1"] <- "log10 H-index + 1"
effect_sizes$Variable[effect_sizes$Variable == "log10 Residuals number of views"] <- "log10 residuals number of Wikipedia views"
effect_sizes$Variable[effect_sizes$Variable == "log10 Number of animals in captivity"] <- "log10 number of animals in captivity"
effect_sizes$Variable[effect_sizes$Variable == "log10 Number of institutions"] <- "log10 number of institutions"
  
fig <- ggplot(data = effect_sizes, aes(y = OR, x  = Variable))+
  geom_hline(yintercept=1, size = 1.5, color = "darkgrey")+
  geom_linerange(aes(ymin = Q2.5, ymax = Q97.5), color = "skyblue", size = 1.5)+
  geom_point(size = 4, color = "red")+
  theme_classic()+
  ylim(0,3.5)+
  theme(axis.text.x = element_markdown(angle = 80, vjust = 1, hjust=1, face="bold"),
        panel.border = element_rect(colour = "black", fill = NA))+
  theme(axis.text=element_text(size=12),
         axis.title=element_text(size=16,face="bold"),
         strip.text = element_text(face="bold", size=13))+
  xlab("Confounding variable")+
  ylab("Odds ratio and 95% credible intervals")+
  ggtitle("Effects on CMR based mortality prevalence")
fig

png("Summary of effect sizes CMR.png", width = 35, height = 25, unit = "cm",res = 300)
fig
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Combining the plots for supplementary material ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

png("Univariate models mortality CMR.png", width = 40, height = 50, unit = "cm",res = 300)
pmortality
dev.off()

