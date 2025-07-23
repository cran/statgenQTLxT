## ----setup, include = FALSE---------------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.dim = c(7, 4)
)
library(statgenQTLxT)
op <- options(width = 110, 
              digits = 2)

## ----addPheno-----------------------------------------------------------------------------------------------
data(dropsPheno, package = "statgenGWAS")

## Convert phenotypic data to a list.
colnames(dropsPheno)[1] <- "genotype"
dropsPheno <- dropsPheno[c("Experiment", "genotype", "grain.yield", "grain.number",
                           "anthesis", "silking", "plant.height", "ear.height")]
## Split data by experiment.
dropsPhenoList <- split(x = dropsPheno, f = dropsPheno[["Experiment"]])
## Remove Experiment column.
## phenotypic data should consist only of genotype and traits.
dropsPhenoList <- lapply(X = dropsPhenoList, FUN = `[`, -1)

## ----createGdata--------------------------------------------------------------------------------------------
## Load marker data.
data("dropsMarkers", package = "statgenGWAS")
## Add genotypes as row names of dropsMarkers and drop Ind column.
rownames(dropsMarkers) <- dropsMarkers[["Ind"]]
dropsMarkers <- dropsMarkers[, -1]

## Load genetic map.
data("dropsMap", package = "statgenGWAS")
## Add genotypes as row names of dropsMap.
rownames(dropsMap) <- dropsMap[["SNP.names"]]
## Rename Chromosome and Position columns.
colnames(dropsMap)[2:3] <- c("chr", "pos")

## Create a gData object containing map, marker and phenotypic information.
gDataDrops <- statgenGWAS::createGData(geno = dropsMarkers,
                                       map = dropsMap, 
                                       pheno = dropsPhenoList)

## ----sumGData-----------------------------------------------------------------------------------------------
## Summarize gDataDrops.
summary(gDataDrops, trials = "Mur13W")

## ----removeDupMarkers---------------------------------------------------------------------------------------
## Set seed.
set.seed(1234)

## Remove duplicate SNPs from gDataDrops.
gDataDropsDedup <- statgenGWAS::codeMarkers(gDataDrops, 
                                            impute = FALSE, 
                                            verbose = TRUE) 

## ----mtg----------------------------------------------------------------------------------------------------
## Run multi-trait GWAS for 5 traits in trial Mur13W.
GWASDrops <- runMultiTraitGwas(gData = gDataDropsDedup, 
                               traits = c("grain.yield","grain.number",
                                          "anthesis", "silking" ,"plant.height"),
                               trials = "Mur13W", 
                               covModel = "fa")

## ----gwaRes-------------------------------------------------------------------------------------------------
head(GWASDrops$GWAResult$Mur13W)

head(GWASDrops$GWAResult$Mur13W[GWASDrops$GWAResult$Mur13W$trait == "grain.yield", ])

## ----signSnp------------------------------------------------------------------------------------------------
GWASDrops$signSnp$Mur13W

## ----sumMtg-------------------------------------------------------------------------------------------------
## Create summary of GWASDrops for the trait grain number.
summary(GWASDrops, traits = "grain.number")

## ----qqMtg--------------------------------------------------------------------------------------------------
## Plot a qq plot of GWAS Drops.
plot(GWASDrops, plotType = "qq")

## ----manhattanMtg-------------------------------------------------------------------------------------------
## Plot a manhattan plot of GWAS Drops.
plot(GWASDrops, plotType = "manhattan")

## ----qtlMtgNorm---------------------------------------------------------------------------------------------
## Plot a qtl plot of GWAS Drops for Mur13W.
## Set significance threshold to 5 and normalize effect estimates.
plot(GWASDrops, plotType = "qtl", yThr = 5, normalize = TRUE)

## ----mtgChrSpec---------------------------------------------------------------------------------------------
## Run multi-trait GWAS for trial 'Mur13W'.
## Use chromosome specific kinship matrices computed using method of van Raden.
GWASDropsChrSpec <- runMultiTraitGwas(gData = gDataDropsDedup, 
                                      trials = "Mur13W",
                                      GLSMethod = "multi",
                                      kinshipMethod = "vanRaden",
                                      covModel = "fa")

## ----addPhenoxE---------------------------------------------------------------------------------------------
## Reshape phenotypic data to data.frame in wide format containing only grain.yield.
phenoDat <- reshape(dropsPheno[, c("Experiment", "genotype", "grain.yield")], 
                    timevar = "Experiment", 
                    idvar = "genotype", 
                    direction = "wide", 
                    v.names = "grain.yield")
## Rename columns to trial name only.
colnames(phenoDat)[2:ncol(phenoDat)] <-
  gsub(pattern = "grain.yield.", replacement = "",
       x = colnames(phenoDat)[2:ncol(phenoDat)])

## ----createGdataxE------------------------------------------------------------------------------------------
## Create a gData object containing map, marker and phenotypic information.
gDataDropsxE <- statgenGWAS::createGData(geno = dropsMarkers,
                                         map = dropsMap, 
                                         pheno = phenoDat)
summary(gDataDropsxE)

## ----removeDupMarkersxE-------------------------------------------------------------------------------------
## Remove duplicate SNPs from gDataDrops.
gDataDropsDedupxE <- statgenGWAS::codeMarkers(gDataDropsxE, 
                                              impute = FALSE,
                                              verbose = TRUE) 

## ----mtgxE--------------------------------------------------------------------------------------------------
## Run multi-trial GWAS for one trait in all trials.
GWASDropsxE <- runMultiTraitGwas(gData = gDataDropsDedupxE, 
                                 covModel = "fa")

## ----signSnpxE----------------------------------------------------------------------------------------------
head(GWASDropsxE$signSnp$pheno, row.names = FALSE)

## ----sumMtgxE-----------------------------------------------------------------------------------------------
summary(GWASDropsxE, traits = c("Mur13W", "Kar12W"))

## ----qqMtgxE------------------------------------------------------------------------------------------------
plot(GWASDropsxE, plotType = "qq")

## ----manhattanMtgxE-----------------------------------------------------------------------------------------
plot(GWASDropsxE, plotType = "manhattan")

## ----qtlMtgNormxE-------------------------------------------------------------------------------------------
## Set significance threshold to 6 and do not normalize effect estimates.
plot(GWASDropsxE, plotType = "qtl", yThr = 6, normalize = FALSE)

## ----mtgSNPFixThr, eval=FALSE-------------------------------------------------------------------------------
# ## Run multi-trait GWAS for Mur13W.
# ## Use a fixed significance threshold of 4.
# GWASDropsFixThr <- runMultiTraitGwas(gData = gDataDropsDedup,
#                                      trials = "Mur13W",
#                                      covModel = "fa",
#                                      thrType = "fixed",
#                                      LODThr = 4)

## ----mtgSNPNR, eval=FALSE-----------------------------------------------------------------------------------
# ## Run multi-trait GWAS for for Mur13W.
# ## Use a factor analytic model for computing the variance components.
# GWASDropsFA <- runMultiTraitGwas(gData = gDataDropsDedup,
#                                  trials = "Mur13W",
#                                  covModel = "fa")
# 
# ## Rerun the analysis, using the variance components computed in the
# ## previous model as inputs.
# GWASDropsFA2 <- runMultiTraitGwas(gData = gDataDropsDedup,
#                                   trials = "Mur13W",
#                                   fitVarComp  = FALSE,
#                                   Vg = GWASDropsFA$GWASInfo$varComp$Vg,
#                                   Ve = GWASDropsFA$GWASInfo$varComp$Ve)

## ----mtgPar, eval = FALSE-----------------------------------------------------------------------------------
# ## Register parallel back-end with 2 cores.
# doParallel::registerDoParallel(cores = 2)
# 
# ## Run multi-trait GWAS for one trait in all trials.
# GWASDropsxEPar <- runMultiTraitGwas(gData = gDataDropsDedupxE,
#                                     covModel = "pw",
#                                     parallel = TRUE)

## ----mtgSNPCovar, eval=FALSE--------------------------------------------------------------------------------
# ## Run multi-trait GWAS for Mur13W.
# ## Use PZE-106021410, the most significant SNP, as SNP covariate.
# GWASDropsSnpCov <- runMultiTraitGwas(gData = gDataDropsDedup,
#                                      trials = "Mur13W",
#                                      snpCov = "PZE-106021410",
#                                      covModel = "fa")

## ----mtgMAF, eval=FALSE-------------------------------------------------------------------------------------
# ## Run multi-trait GWAS for Mur13W.
# ## Only include SNPs that have a MAF of 0.05 or higher.
# GWASDropsMAF <- runMultiTraitGwas(gData = gDataDropsDedup,
#                                   trials = "Mur13W",
#                                   covModel = "fa",
#                                   MAF = 0.05)

## ----mtgCommon----------------------------------------------------------------------------------------------
## Run multi-trait GWAS for Mur13W.
## Fit an additional common sNP effect model.
GWASDropsCommon <- runMultiTraitGwas(gData = gDataDropsDedup,
                                     trials = "Mur13W",
                                     covModel = "fa",
                                     estCom = TRUE)
head(GWASDropsCommon$GWAResult$Mur13W)

## ----winddown, include = FALSE------------------------------------------------
options(op)

