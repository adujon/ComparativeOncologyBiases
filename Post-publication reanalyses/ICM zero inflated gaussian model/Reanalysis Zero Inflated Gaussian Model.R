rm(list=ls()) #remove all objects into memmory
library(readxl)
library(ape)
library(nlme)
library(geiger)
library(caper)
library(phylolm)
library(phytools)
library(MuMIn)
library(ggplot2)
library(rr2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(grid)
library(parallel)
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
setwd("C:/Users/antoi/My Drive/Postdoc/Zoo data/Biases/reanalysis pageview proxy/Post-publication reanalyses")
dataICM <- read_xlsx("Raw dataset confounding variables.xlsx", sheet = "Mortality_ICM", col_types = c("text", rep("numeric", 20), "text", "numeric", "numeric", "numeric"))
dataICM <- data.frame(dataICM)
rownames(dataICM) <- dataICM$Species
dataICM$logitICM <- car::logit(dataICM$ICM) #logit transform the ICM
dataICM$W <- log(dataICM$NNecropsies) #calculate the weights of the regression
dataICM$occurance <- ifelse(dataICM$ICM == 0, 0 ,1)
dataICM$knownDeaths <-   dataICM$NNecropsies

#create an empty dataframe to store the effect size for each confounding variables
effect_sizes <- data.frame(OR  = NULL, Q2.5 = NULL, Q97.5 = NULL, Variable = NULL, Tumour = NULL)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Prepare var covar matrixes for brms models ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#check which species are both in the main dataframe and the phylogeny tree
obj_mortality<- name.check(phylo_all, dataICM)
obj_mortality

#drop the species from the tree for which we have no cancer data
phylo_mortality <- drop.tip(phylo_all, obj_mortality$tree_not_data)

#compute the var-covar matrix from the tree to use in brms
phy <- phylo_mortality

row.names(dataICM) <- dataICM$Species
phyICM <- treedata(phy, dataICM, warnings = FALSE)$phy
dataICM$occurrence <- ifelse(dataICM$ICM == 0, 0, 1)
dataICM <- dataICM[as.character(phyICM$tip.label), ] # organise data to follow same order as phylogeny


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Implement Zero-inflated phylogenetic model (from Vincze et al') ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Zero-inflated phylogenetic model
ZIlogis <- function(formbin,   # model formula for binomial model
                    formlogis, # model formula for logistic model
                    databin,   # dataset for binomial GLM repsonse variable
                    datalogis, # dataset with logistic response variable
                    phyloTree, # phylogeneitc tree for logistic regression
                    phylobin,  # phylogenetic tree for binomial regression
                    l0) {      # starting value of Pagel's lambda in PGLS regression
  out <- list()
  out$bin <- binaryPGLMM(formbin, data = databin, phy = phylobin)
  out$wlogis <- gls(formlogis,
                    correlation = corPagel(l0, phy = phyloTree, fixed = FALSE, form=~Species),
                    data = datalogis, weights = ~1/log(knownDeaths) )
  if(out$wlogis$modelStruct[1] < 0){ # refit model with lanbda=0 if lambda converged to negative
    out$wlogis <- gls(formlogis,
                      correlation = corPagel(0, phy = phyloTree, fixed = TRUE, form=~Species),
                      data = datalogis, weights = ~1/log(knownDeaths))}
  class(out) <- "ZILogis"
  return(out)}


# Summary function for the zero-inflated phylogenetic model
sumZILogis <- function(x,...){
  x1 <- cbind(as.data.frame(x$bin[[2]]), as.data.frame(x$bin[3]), as.data.frame(x$bin[5]), as.data.frame(x$bin[6]))
  x1[,1:3] <- round(x1[,1:3],2)
  x1[,4] <- round(x1[,4],4)
  x1[,1] <- paste(x1[,1]," (", x1[,2],")", sep='')
  x1 <- x1[,-2]
  x1b <- cbind(row.names(x1),x1[,1:3]);
  row.names(x1b) <- NULL
  x1b[nrow(x1)+1, 4] <- paste("n=", nrow(as.data.frame(x$bin[9])))
  x1b[nrow(x1)+1, 3] <- paste("s2=", round(as.numeric(x$bin[7]),2))
  x1b[nrow(x1)+1, 2] <- paste("P_s2=", round(as.numeric(x$bin[8]),4))
  x1b[,1] <- as.character(x1b[,1])
  x1b[nrow(x1)+1, 1] <- "ModelStats"
  x2w <- as.data.frame(summary(x$wlogis)$tTable); names(x2w) <- names(x1)
  x2w[,1:3] <- round(x2w[,1:3],2)
  x2w[,4] <- round(x2w[,4],4)
  x2w[,1] <- paste(x2w[,1]," (", x2w[,2],")", sep='')
  x2w <- x2w[,-2]
  x2wb <- cbind(row.names(x2w),x2w[,1:3]);
  row.names(x2wb) <- NULL
  x2wb[nrow(x2w)+1, 2] <- paste("AIC=", round(summary(x$wlogis)$AIC,2))
  x2wb[nrow(x2w)+1, 4] <- paste("n=", length(x$wlogis$residuals))
  x2wb[nrow(x2w)+1, 3] <- paste("L=", round(unlist(summary(x$wlogis))$modelStruct.corStruct,2))
  x2wb[,1] <- as.character(x2wb[,1])
  x2wb[nrow(x2w)+1, 1] <- "ModelStats"
  names(x2wb) <- names(x1b)
  names(x1b)[1] <- ""
  x3 <- x1b[1,]; x3[,1:4] <- ""; x3[1,1] <- "Probability of zeros"
  names(x2wb)[1] <- ""
  x5 <- x3; x5[1,1] <- "Weighted logistic reg."
  X <- rbind(x3,x1b,x5,x2wb)
  names(X)[1:4] <- c('_',"b (SE)","z/t-value", "p-value")
  return(X)        }

#function to get the effect sizes as a table
effectsZILogis <- function(x) {
  
  # ---- Binomial (probability of zero / occurrence) part ----
  B    <- as.data.frame(x$bin$B)
  B.se <- as.data.frame(x$bin$B.se)[,1]
  B.z  <- as.data.frame(x$bin$B.zscore)[,1]
  B.p  <- as.data.frame(x$bin$B.pvalue)[,1]
  B.n  <- length(x$bin$mu)   # number of species used in the binomial (PGLMM) model
  
  bin_df <- data.frame(
    Term      = rownames(B),
    Estimate  = B[,1],
    SE        = B.se,
    Statistic = B.z,
    P         = B.p,
    N         = B.n,
    Model     = "Probability of zeros",
    stringsAsFactors = FALSE
  )
  
  # response variable name, pulled straight from the binomial model formula
  bin_response <- all.vars(x$bin$formula)[1]
  
  # ---- Gaussian (weighted logistic/PGLS) part ----
  tT <- as.data.frame(summary(x$wlogis)$tTable)
  names(tT) <- c("Estimate", "SE", "Statistic", "P")
  gau_n <- length(x$wlogis$residuals)  # number of species used in the gls model
  
  gau_df <- data.frame(
    Term      = rownames(tT),
    Estimate  = tT$Estimate,
    SE        = tT$SE,
    Statistic = tT$Statistic,
    P         = tT$P,
    N         = gau_n,
    Model     = "Weighted logistic reg.",
    stringsAsFactors = FALSE
  )
  
  # response variable name, pulled straight from the gls model formula
  gau_response <- all.vars(formula(x$wlogis))[1]
  
  out <- rbind(bin_df, gau_df)
  out$Tumour <- ifelse(out$Model == "Probability of zeros", bin_response, gau_response)
  rownames(out) <- NULL
  return(out)
}

#empty data frame to store the effect sizes
effect_sizes <- data.frame(
  Term      = character(),
  Estimate  = numeric(),
  SE        = numeric(),
  Statistic = numeric(),
  P         = numeric(),
  Model     = character(),
  OR        = numeric(),
  Q2.5      = numeric(),
  Q97.5     = numeric(),
  Tumour    = character(),
  N = numeric(),
  stringsAsFactors = FALSE
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Tumor prevalence ~ mean_GDP_5years ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Zero-inflated COMPOSITE models
# --- NA-handling for the binomial part ---
complete_sp <- dataICM$Species[complete.cases(dataICM[, c("occurrence", "mean_GDP_5years")])]

phyICM_bin  <- drop.tip(phyICM, setdiff(phyICM$tip.label, complete_sp))
dataICM_bin <- dataICM[phyICM_bin$tip.label, ]

stopifnot(nrow(dataICM_bin) == length(phyICM_bin$tip.label))
stopifnot(all(dataICM_bin$Species == phyICM_bin$tip.label))
# --- end NA-handling (binomial) ---

# Data and phylogeny for linear phylogenetic model of non-zero cancer mortality risk
phy_nonzeroICM  <- drop.tip(phyICM, which(!dataICM$occurrence == 1))
data_nonzeroICM <- dataICM[dataICM$Species %in% phy_nonzeroICM$tip.label, ]

# --- NA / non-finite handling for the gls() part ---
ok_logis <- with(data_nonzeroICM,
                  is.finite(logitICM) &
                  is.finite(mean_GDP_5years) &
                  is.finite(knownDeaths) &
                  knownDeaths > 1)              # log(knownDeaths) must be > 0

complete_sp_logis <- data_nonzeroICM$Species[ok_logis]

phy_nonzeroICM  <- drop.tip(phy_nonzeroICM, setdiff(phy_nonzeroICM$tip.label, complete_sp_logis))
data_nonzeroICM <- data_nonzeroICM[phy_nonzeroICM$tip.label, ]

stopifnot(nrow(data_nonzeroICM) == length(phy_nonzeroICM$tip.label))
stopifnot(all(data_nonzeroICM$Species == phy_nonzeroICM$tip.label))
# --- end NA-handling (gls) ---

# quick diagnostics — worth checking once before fitting
cat("Binomial part:", nrow(dataICM_bin), "species (of", length(phyICM$tip.label), "original)\n")
cat("Gaussian part:", nrow(data_nonzeroICM), "species (of",
    length(drop.tip(phyICM, which(!dataICM$occurrence == 1))$tip.label), "non-zero species)\n")

#model
modICM_mean_GDP_5years <- ZIlogis(formbin   = occurrence ~ mean_GDP_5years+log(knownDeaths),
                   formlogis = logitICM ~ mean_GDP_5years,
                   databin   = dataICM_bin,
                   phylobin  = phyICM_bin,
                   datalogis = data_nonzeroICM,
                   phyloTree = phy_nonzeroICM,
                   l0 = 0)

#summary of the zero inflateed model
sumZILogis(modICM_mean_GDP_5years)
effectsZILogis(modICM_mean_GDP_5years)

#store effect size into a table
effect_sizes<- rbind(effect_sizes, effectsZILogis(modICM_mean_GDP_5years))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Tumor prevalence ~ logNViews ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Zero-inflated COMPOSITE models
# --- NA-handling for the binomial part ---
complete_sp <- dataICM$Species[complete.cases(dataICM[, c("occurrence", "logNViews")])]

phyICM_bin  <- drop.tip(phyICM, setdiff(phyICM$tip.label, complete_sp))
dataICM_bin <- dataICM[phyICM_bin$tip.label, ]

stopifnot(nrow(dataICM_bin) == length(phyICM_bin$tip.label))
stopifnot(all(dataICM_bin$Species == phyICM_bin$tip.label))
# --- end NA-handling (binomial) ---

# Data and phylogeny for linear phylogenetic model of non-zero cancer mortality risk
phy_nonzeroICM  <- drop.tip(phyICM, which(!dataICM$occurrence == 1))
data_nonzeroICM <- dataICM[dataICM$Species %in% phy_nonzeroICM$tip.label, ]

# --- NA / non-finite handling for the gls() part ---
ok_logis <- with(data_nonzeroICM,
                 is.finite(logitICM) &
                   is.finite(logNViews) &
                   is.finite(knownDeaths) &
                   knownDeaths > 1)              # log(knownDeaths) must be > 0

complete_sp_logis <- data_nonzeroICM$Species[ok_logis]

phy_nonzeroICM  <- drop.tip(phy_nonzeroICM, setdiff(phy_nonzeroICM$tip.label, complete_sp_logis))
data_nonzeroICM <- data_nonzeroICM[phy_nonzeroICM$tip.label, ]

stopifnot(nrow(data_nonzeroICM) == length(phy_nonzeroICM$tip.label))
stopifnot(all(data_nonzeroICM$Species == phy_nonzeroICM$tip.label))
# --- end NA-handling (gls) ---

# quick diagnostics — worth checking once before fitting
cat("Binomial part:", nrow(dataICM_bin), "species (of", length(phyICM$tip.label), "original)\n")
cat("Gaussian part:", nrow(data_nonzeroICM), "species (of",
    length(drop.tip(phyICM, which(!dataICM$occurrence == 1))$tip.label), "non-zero species)\n")

#model
modICM_logNViews <- ZIlogis(formbin   = occurrence ~ logNViews + log(knownDeaths),
                                  formlogis = logitICM ~ logNViews,
                                  databin   = dataICM_bin,
                                  phylobin  = phyICM_bin,
                                  datalogis = data_nonzeroICM,
                                  phyloTree = phy_nonzeroICM,
                                  l0 = 0)

#summary of the zero inflateed model
sumZILogis(modICM_logNViews)
effectsZILogis(modICM_logNViews)

#store effect size into a table
effect_sizes<- rbind(effect_sizes, effectsZILogis(modICM_logNViews))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Tumor prevalence ~ ResidualslogNViews ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Zero-inflated COMPOSITE models
# --- NA-handling for the binomial part ---
complete_sp <- dataICM$Species[complete.cases(dataICM[, c("occurrence", "ResidualslogNViews")])]

phyICM_bin  <- drop.tip(phyICM, setdiff(phyICM$tip.label, complete_sp))
dataICM_bin <- dataICM[phyICM_bin$tip.label, ]

stopifnot(nrow(dataICM_bin) == length(phyICM_bin$tip.label))
stopifnot(all(dataICM_bin$Species == phyICM_bin$tip.label))
# --- end NA-handling (binomial) ---

# Data and phylogeny for linear phylogenetic model of non-zero cancer mortality risk
phy_nonzeroICM  <- drop.tip(phyICM, which(!dataICM$occurrence == 1))
data_nonzeroICM <- dataICM[dataICM$Species %in% phy_nonzeroICM$tip.label, ]

# --- NA / non-finite handling for the gls() part ---
ok_logis <- with(data_nonzeroICM,
                 is.finite(logitICM) &
                   is.finite(ResidualslogNViews) &
                   is.finite(knownDeaths) &
                   knownDeaths > 1)              # log(knownDeaths) must be > 0

complete_sp_logis <- data_nonzeroICM$Species[ok_logis]

phy_nonzeroICM  <- drop.tip(phy_nonzeroICM, setdiff(phy_nonzeroICM$tip.label, complete_sp_logis))
data_nonzeroICM <- data_nonzeroICM[phy_nonzeroICM$tip.label, ]

stopifnot(nrow(data_nonzeroICM) == length(phy_nonzeroICM$tip.label))
stopifnot(all(data_nonzeroICM$Species == phy_nonzeroICM$tip.label))
# --- end NA-handling (gls) ---

# quick diagnostics — worth checking once before fitting
cat("Binomial part:", nrow(dataICM_bin), "species (of", length(phyICM$tip.label), "original)\n")
cat("Gaussian part:", nrow(data_nonzeroICM), "species (of",
    length(drop.tip(phyICM, which(!dataICM$occurrence == 1))$tip.label), "non-zero species)\n")

#model
modICM_ResidualslogNViews <- ZIlogis(formbin   = occurrence ~ ResidualslogNViews+log(knownDeaths),
                            formlogis = logitICM ~ ResidualslogNViews,
                            databin   = dataICM_bin,
                            phylobin  = phyICM_bin,
                            datalogis = data_nonzeroICM,
                            phyloTree = phy_nonzeroICM,
                            l0 = 0)

#summary of the zero inflateed model
sumZILogis(modICM_ResidualslogNViews)
effectsZILogis(modICM_ResidualslogNViews)

#store effect size into a table
effect_sizes<- rbind(effect_sizes, effectsZILogis(modICM_ResidualslogNViews))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Tumor prevalence ~ logHIndex ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Zero-inflated COMPOSITE models
# --- NA-handling for the binomial part ---
complete_sp <- dataICM$Species[complete.cases(dataICM[, c("occurrence", "logHIndex")])]

phyICM_bin  <- drop.tip(phyICM, setdiff(phyICM$tip.label, complete_sp))
dataICM_bin <- dataICM[phyICM_bin$tip.label, ]

stopifnot(nrow(dataICM_bin) == length(phyICM_bin$tip.label))
stopifnot(all(dataICM_bin$Species == phyICM_bin$tip.label))
# --- end NA-handling (binomial) ---

# Data and phylogeny for linear phylogenetic model of non-zero cancer mortality risk
phy_nonzeroICM  <- drop.tip(phyICM, which(!dataICM$occurrence == 1))
data_nonzeroICM <- dataICM[dataICM$Species %in% phy_nonzeroICM$tip.label, ]

# --- NA / non-finite handling for the gls() part ---
ok_logis <- with(data_nonzeroICM,
                 is.finite(logitICM) &
                   is.finite(logHIndex) &
                   is.finite(knownDeaths) &
                   knownDeaths > 1)              # log(knownDeaths) must be > 0

complete_sp_logis <- data_nonzeroICM$Species[ok_logis]

phy_nonzeroICM  <- drop.tip(phy_nonzeroICM, setdiff(phy_nonzeroICM$tip.label, complete_sp_logis))
data_nonzeroICM <- data_nonzeroICM[phy_nonzeroICM$tip.label, ]

stopifnot(nrow(data_nonzeroICM) == length(phy_nonzeroICM$tip.label))
stopifnot(all(data_nonzeroICM$Species == phy_nonzeroICM$tip.label))
# --- end NA-handling (gls) ---

# quick diagnostics — worth checking once before fitting
cat("Binomial part:", nrow(dataICM_bin), "species (of", length(phyICM$tip.label), "original)\n")
cat("Gaussian part:", nrow(data_nonzeroICM), "species (of",
    length(drop.tip(phyICM, which(!dataICM$occurrence == 1))$tip.label), "non-zero species)\n")

#model
modICM_logHIndex <- ZIlogis(formbin   = occurrence ~ logHIndex+log(knownDeaths),
                                     formlogis = logitICM ~ logHIndex,
                                     databin   = dataICM_bin,
                                     phylobin  = phyICM_bin,
                                     datalogis = data_nonzeroICM,
                                     phyloTree = phy_nonzeroICM,
                                     l0 = 0)

#summary of the zero inflateed model
sumZILogis(modICM_logHIndex)
effectsZILogis(modICM_logHIndex)

#store effect size into a table
effect_sizes<- rbind(effect_sizes, effectsZILogis(modICM_logHIndex))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Tumor prevalence ~ Number of institutions ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Zero-inflated COMPOSITE models
# --- NA-handling for the binomial part ---
complete_sp <- dataICM$Species[complete.cases(dataICM[, c("occurrence", "logNInstitutions")])]
phyICM_bin  <- drop.tip(phyICM, setdiff(phyICM$tip.label, complete_sp))
dataICM_bin <- dataICM[phyICM_bin$tip.label, ]

stopifnot(nrow(dataICM_bin) == length(phyICM_bin$tip.label))
stopifnot(all(dataICM_bin$Species == phyICM_bin$tip.label))
# --- end NA-handling (binomial) ---

# Data and phylogeny for linear phylogenetic model of non-zero cancer mortality risk
phy_nonzeroICM  <- drop.tip(phyICM, which(!dataICM$occurrence == 1))
data_nonzeroICM <- dataICM[dataICM$Species %in% phy_nonzeroICM$tip.label, ]

# --- NA / non-finite handling for the gls() part ---
ok_logis <- with(data_nonzeroICM,
                 is.finite(logitICM) &
                   is.finite(logNInstitutions) &
                   is.finite(knownDeaths) &
                   knownDeaths > 1)              # log(knownDeaths) must be > 0

complete_sp_logis <- data_nonzeroICM$Species[ok_logis]

phy_nonzeroICM  <- drop.tip(phy_nonzeroICM, setdiff(phy_nonzeroICM$tip.label, complete_sp_logis))
data_nonzeroICM <- data_nonzeroICM[phy_nonzeroICM$tip.label, ]

stopifnot(nrow(data_nonzeroICM) == length(phy_nonzeroICM$tip.label))
stopifnot(all(data_nonzeroICM$Species == phy_nonzeroICM$tip.label))
# --- end NA-handling (gls) ---

# quick diagnostics — worth checking once before fitting
cat("Binomial part:", nrow(dataICM_bin), "species (of", length(phyICM$tip.label), "original)\n")
cat("Gaussian part:", nrow(data_nonzeroICM), "species (of",
    length(drop.tip(phyICM, which(!dataICM$occurrence == 1))$tip.label), "non-zero species)\n")

#model
modICM_logNInstitutions <- ZIlogis(formbin   = occurrence ~ logNInstitutions+log(knownDeaths),
                            formlogis = logitICM ~ logNInstitutions,
                            databin   = dataICM_bin,
                            phylobin  = phyICM_bin,
                            datalogis = data_nonzeroICM,
                            phyloTree = phy_nonzeroICM,
                            l0 = 0)

#summary of the zero inflateed model
sumZILogis(modICM_logNInstitutions)
effectsZILogis(modICM_logNInstitutions)

#store effect size into a table
effect_sizes<- rbind(effect_sizes, effectsZILogis(modICM_logNInstitutions))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Tumor prevalence ~ Number of ExSitu individuals ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Zero-inflated COMPOSITE models
# --- NA-handling for the binomial part ---
complete_sp <- dataICM$Species[complete.cases(dataICM[, c("occurrence", "logNExSitu")])]
phyICM_bin  <- drop.tip(phyICM, setdiff(phyICM$tip.label, complete_sp))
dataICM_bin <- dataICM[phyICM_bin$tip.label, ]

stopifnot(nrow(dataICM_bin) == length(phyICM_bin$tip.label))
stopifnot(all(dataICM_bin$Species == phyICM_bin$tip.label))
# --- end NA-handling (binomial) ---

# Data and phylogeny for linear phylogenetic model of non-zero cancer mortality risk
phy_nonzeroICM  <- drop.tip(phyICM, which(!dataICM$occurrence == 1))
data_nonzeroICM <- dataICM[dataICM$Species %in% phy_nonzeroICM$tip.label, ]

# --- NA / non-finite handling for the gls() part ---
ok_logis <- with(data_nonzeroICM,
                 is.finite(logitICM) &
                   is.finite(logNExSitu) &
                   is.finite(knownDeaths) &
                   knownDeaths > 1)              # log(knownDeaths) must be > 0

complete_sp_logis <- data_nonzeroICM$Species[ok_logis]

phy_nonzeroICM  <- drop.tip(phy_nonzeroICM, setdiff(phy_nonzeroICM$tip.label, complete_sp_logis))
data_nonzeroICM <- data_nonzeroICM[phy_nonzeroICM$tip.label, ]

stopifnot(nrow(data_nonzeroICM) == length(phy_nonzeroICM$tip.label))
stopifnot(all(data_nonzeroICM$Species == phy_nonzeroICM$tip.label))
# --- end NA-handling (gls) ---

# quick diagnostics — worth checking once before fitting
cat("Binomial part:", nrow(dataICM_bin), "species (of", length(phyICM$tip.label), "original)\n")
cat("Gaussian part:", nrow(data_nonzeroICM), "species (of",
    length(drop.tip(phyICM, which(!dataICM$occurrence == 1))$tip.label), "non-zero species)\n")

#model
modICM_logNExSitu <- ZIlogis(formbin   = occurrence ~ logNExSitu+log(knownDeaths),
                                   formlogis = logitICM ~ logNExSitu,
                                   databin   = dataICM_bin,
                                   phylobin  = phyICM_bin,
                                   datalogis = data_nonzeroICM,
                                   phyloTree = phy_nonzeroICM,
                                   l0 = 0)

#summary of the zero inflateed model
sumZILogis(modICM_logNExSitu)
effectsZILogis(modICM_logNExSitu)
               
#store effect size into a table
effect_sizes<- rbind(effect_sizes, effectsZILogis(modICM_logNExSitu))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Clean the table ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#remove intercept and calculate odds ratio
effect_sizes <-subset(effect_sizes, Term != "(Intercept)")
effect_sizes <-subset(effect_sizes, Term != "log(knownDeaths)")
effect_sizes$OR <- exp(effect_sizes$Estimate)
effect_sizes$Q2.5 <- exp(effect_sizes$Estimate - 1.96 * effect_sizes$SE)
effect_sizes$Q97.5 <- exp(effect_sizes$Estimate + 1.96 * effect_sizes$SE)

#round values
effect_sizes[, 2:4] <- round(effect_sizes[, 2:4], digits = 2)
effect_sizes$P <- round(effect_sizes$P, digits = 5)
effect_sizes[, 9:11] <- round(effect_sizes[, 9:11], digits = 2)

#reorder dataframe
effect_sizes <- effect_sizes[, c(
  "Tumour",
  "Model",
  "Term",
  "Estimate",
  "SE",
  "Statistic",
  "P",
  "N",
  "OR",
  "Q2.5",
  "Q97.5"
)]

effect_sizes$Significant <- ifelse(effect_sizes$P < 0.05, "YES", "NO")
effect_sizes$Variable <- effect_sizes$Term
effect_sizes <- effect_sizes[order(effect_sizes$Tumour), ]
effect_sizes

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create a plot of effect sizes ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

effect_sizes <- effect_sizes %>%
  mutate(
    Tumour = factor(
      Tumour,
      levels = c("occurrence", "logitICM"),
      labels = c(
        "Probability of detecting cancer mortality\nin at least one individual",
        "Weighted phylogenetic regression\nfor species in which cancer mortality was detected"
      )
    ),
    across(c(OR, Q2.5, Q97.5), as.numeric),
    Variable = recode(
      Term,
      mean_GDP_5years    = "Weighted average GDP <i>per capita</i><br>(10,000s USD)",
      logNViews          = "log10 Wikipedia views",
      ResidualslogNViews = "log10 residual Wikipedia views",
      logHIndex          = "log10 H-index + 1",
      logNInstitutions   = "log10 institutions",
      logNExSitu         = "log10 animals in captivity"
    ),
    Variable = factor(
      Variable,
      levels = c(
        "log10 H-index + 1",
        "log10 animals in captivity",
        "log10 institutions",
        "log10 Wikipedia views",
        "log10 residual Wikipedia views",
        "Weighted average GDP <i>per capita</i><br>(10,000s USD)"
      )
    )
  )

fig <- ggplot(effect_sizes, aes(Variable, OR)) +
  geom_hline(yintercept = 1, linewidth = 1.2, colour = "darkgrey") +
  geom_linerange(
    aes(ymin = Q2.5, ymax = Q97.5),
    linewidth = 1.3,
    colour = "skyblue"
  ) +
  geom_point(aes(shape = Significant), size = 4, colour = "red") +
  facet_wrap(~Tumour, ncol = 1, scales = "free_y") +
  scale_shape_manual(values = c(YES = 16, NO = 1)) +
  labs(
    x = "Confounding variable",
    y = "Odds ratio and 95% confidence intervals",
    title = paste0(
      "Associations between popularity metrics and ICM\n",
      "(modelled as a continuous variable)"
    ),
    shape = "Significant"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_markdown(
      angle = 70, hjust = 1, vjust = 1,
      face = "bold", size = 11
    ),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 13, face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA),
    plot.title = element_text(size = 17, face = "bold", hjust = 0.5)
  )
fig


setwd("C:/Users/antoi/My Drive/Postdoc/Zoo data/Biases/reanalysis pageview proxy/Post-publication reanalyses/ICM zero inflated gaussian model")

ggsave(
  filename = "ICM_popularity_effect_sizes_zero_inflated_model.png",
  plot = fig,
  width = 8,
  height = 10,
  units = "in",
  dpi = 600
)

