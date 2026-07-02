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
library(scales)

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
data_mortality <- read_xlsx("Raw dataset confounding variables.xlsx", sheet = "Mortality_ICM", col_types = c("text", rep("numeric", 20), "text", "numeric", "numeric"))
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
# Agreement between two number of views metrics ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#function to help puting spaces between 000's
fmt_space <- function(x) format(x, big.mark = " ", scientific = FALSE, trim = TRUE)

#computing r-square
r_sq_val <- cor(data_mortality$Wikipedia_Views,
                data_mortality$Wikipedia_Views_Langviews,
                use = "complete.obs")^2

#plot the two Wikipedia views estimates
p_comp <- ggplot(data_mortality, aes(Wikipedia_Views, Wikipedia_Views_Langviews)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
  geom_point(size = 2) +
  annotate("text", x = -Inf, y = Inf,
           label = paste0("R² = ", round(r_sq_val, 3)),
           hjust = -0.1, vjust = 1.4, size = 8) +
  scale_x_continuous(labels = fmt_space) +
  scale_y_continuous(labels = fmt_space) +
  theme_classic() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 20, face = "bold"),
        plot.title = element_text(size = 20, face = "bold")) +
  labs(title = "Correlation between Wikipedia views from\ncustom R script and views from Langviews tool",
       x = "Wikipedia views (custom R script)",
       y = "Wikipedia views (Langviews tool)")

p_comp

png("Custom VS Langviews Wikipedia correlation ICM.png", width = 20, height = 20, unit = "cm",res = 300)
p_comp
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Tumor prevalence ~ logNViews ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Mortality - Wikimedia data ---- 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#fit multivariate model predicting mortality prevalence prevalence in function of logNViews (obtained from custom R script)
fit_logNViews_mortality_custom <- brm(NMortality | trials(NNecropsies)  ~ logNViews + (1|gr(Species, cov = vcv_brms_OU)),
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

fit_logNViews_mortality_custom

#calculate odd ratios for mortality
rsq <- rstantools::bayes_R2(fit_logNViews_mortality_custom, re.form = NA)
fe <- format(round(exp(fixed.effects(fit_logNViews_mortality_custom)), digits = 1), nsmall = 1)
OR <- paste("OR = ", fe[2,1], " (95% CI:", fe[2,3], ", ", fe[2,4], ")", sep = "")
rsq_lbl <- as.expression(bquote(~.(OR)))
es_mortality<-  data.frame(OR  = fe[2,1], Q2.5 = fe[2,3], Q97.5 = fe[2,4], Variable = "log10 Number of views", Tumour = "Mortality")
effect_sizes <- rbind(effect_sizes, es_mortality)

#calculate marginal effect for mortality
cond_effects <- conditional_effects(fit_logNViews_mortality_custom , effect = "logNViews", resp = "NMortality")
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
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=20,face="bold"))+
  theme(legend.position = c(0.2, 0.5))+
  theme(plot.title = element_text(size = 20, face = "bold"))+
  theme(legend.text=element_text(size=12),
        legend.title = element_text(size = 12))+
  ggtitle("Wikipedia views (based on custom R script)")+
  xlab("log10 number of Wikipedia views")+
  ylab("Mortality prevalence (ICM based)")+
  ylim(0,0.70)+
  xlim(3,8.5)+
  annotate(geom="text", x=5, y=0.68, label=rsq_lbl, size = 8)
p_logNViews_mortality

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Mortality - Pageviews data ---- 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#fit multivariate model predicting mortality prevalence prevalence in function of logNViews (obtained from langview tools)
p_logNViews_mortality_langviews  <- brm(NMortality | trials(NNecropsies)  ~ logNViews_Langviews + (1|gr(Species, cov = vcv_brms_OU)),
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

p_logNViews_mortality_langviews 

#calculate odd ratios for mortality
rsq <- rstantools::bayes_R2(p_logNViews_mortality_langviews , re.form = NA)
fe <- format(round(exp(fixed.effects(p_logNViews_mortality_langviews)), digits = 1), nsmall = 1)
OR <- paste("OR = ", fe[2,1], " (95% CI:", fe[2,3], ", ", fe[2,4], ")", sep = "")
rsq_lbl <- as.expression(bquote(~.(OR)))
es_mortality<-  data.frame(OR  = fe[2,1], Q2.5 = fe[2,3], Q97.5 = fe[2,4], Variable = "log10 Number of views", Tumour = "Mortality")
effect_sizes <- rbind(effect_sizes, es_mortality)

#calculate marginal effect for mortality
cond_effects <- conditional_effects(p_logNViews_mortality_langviews , effect = "logNViews_Langviews", resp = "NMortality")
cond_effects_plot <- data.frame(cond_effects[[1]])
cond_effects_plot$estimate__ <- cond_effects_plot$estimate__/cond_effects_plot$NNecropsies[1]
cond_effects_plot$lower__ <- cond_effects_plot$lower__/cond_effects_plot$NNecropsies[1]
cond_effects_plot$upper__ <- cond_effects_plot$upper__/cond_effects_plot$NNecropsies[1]

#plot marginal effect for mortality
p_logNViews_mortality_langviews <- ggplot()+
  geom_ribbon(data = cond_effects_plot, mapping = aes(x = logNViews_Langviews, ymin = lower__, ymax = upper__), alpha = 0.25)+
  geom_line(data = cond_effects_plot, mapping = aes(x = logNViews_Langviews, y = estimate__), size = 1.5)+
  geom_point(data =data_mortality, mapping = aes(x = logNViews_Langviews, y = MortalityPrevalence), size = 3.5, alpha = 0.75, color = "limegreen")+
  theme_classic()+
  theme_classic()+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=20,face="bold"))+
  theme(legend.position = c(0.2, 0.5))+
  theme(plot.title = element_text(size = 20, face = "bold"))+
  theme(legend.text=element_text(size=12),
        legend.title = element_text(size = 12))+
  ggtitle("Wikipedia views (based on Langviews tool)")+
  xlab("log10 number of Wikipedia views")+
  ylab("Mortality prevalence (ICM based)")+
  ylim(0,0.70)+
  xlim(3,8.5)+
  annotate(geom="text", x=5, y=0.68, label=rsq_lbl, size = 8)
p_logNViews_mortality_langviews 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create a plot to export ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pfigure <- ggarrange(
  
  p_logNViews_mortality, 
  p_logNViews_mortality_langviews,
  common.legend = T,
  legend = "bottom",
  nrow = 1,
  ncol = 2,
  align = "h"
)
pfigure


png("Custom VS Langviews Wikipedia views model comparison ICM.png", width = 40, height = 20, unit = "cm",res = 300)
pfigure
dev.off()