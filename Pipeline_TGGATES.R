##################################################
###
### Pipeline needs to be split into multiple parts:
### - PipelineA (function for TGGATES accession)
### downloads, normalizes, curates TGGATES data
### for specific id. Then performs fitting of linear
### models and pre-ranking of genes.
### - PipelineB (after PipelineA was run for all
### accessions) creates gene sets.
### - PipelineC (function again applied for accession)
### performs GSEA with gene sets from B.
### - PipelineD (function applied for accession)
### performs bi-clustering of normalized enrichment
### scores (NES) provided by GSEA, searches for
### overlap in modules across accessions
### 
### patrick@jimmy.harvard.edu, February 2014
###
### (Starting editing to include kidney data as
### response of revision on Nov 4, 2014.)
###
##################################################

#################################################
############ #### prerequisites ##### ###########
#################################################

rm(list=ls(all=TRUE))

dir.create("scripts", showWarnings=F)
dir.create(file.path("rdata", "intermediate"), recursive=T, showWarnings=F)
dir.create("plots", showWarnings=F)
dir.create(file.path("GSEA", "ranks"), recursive=T, showWarnings=F)
dir.create(file.path("GSEA", "output"), showWarnings=F)
dir.create(file.path("GSEA", "output_reactome"), showWarnings=F)
dir.create("data", showWarnings=F)
dir.create(file.path("plots", "heatmaps"), showWarnings=F)
dir.create("Supplementary_latex", showWarnings=F)

#-----------------------
### global variables ###
#-----------------------

## shortcuts
scriptDir <- file.path("scripts")
rdataDir <- file.path("rdata")
plotDir <- file.path("plots")
intermediateDir <- file.path(rdataDir, "intermediate")
gseaDir <- file.path("GSEA")
rankDir <- file.path(gseaDir,"ranks")
dataDir <- file.path("data")
latexDir <- file.path("Supplementary_latex")

# GSEA variables
gseaOutDir <- file.path(gseaDir,"output")
gseaOutDir_reactome <- file.path(gseaDir,"output_reactome")
gseaExecPath <- file.path("gsea2-2.0.14.jar")
gseaMinSize <- 15
gseaMaxSize <- 500
gseaNperm <- 1000

# aux function
sourcef <- function(file) source(file.path(scriptDir, file))

ncores <- 10

## booleans
installPharmacoGx <- FALSE
installMetaGx <- FALSE
installGenefu <- FALSE

cleanIntermediates <- FALSE
cleanGMTs <- FALSE
cleanGSEAout <- FALSE
cleanRdata <- FALSE
cleanRanks <- FALSE
cleanOnly <- FALSE # if TRUE, pipeline is aborted after cleaning step
cleanBiclusteringSteps <- FALSE

verbose <- TRUE

#--------------------
### libraries ###
#--------------------

## for loading TGGATES data, GSEA, ...
if (installPharmacoGx) {
  require(devtools)
  #   dev_mode(on=T)
  install_github("PharmacoGx", username="bhaibeka", ref="master")
  #   dev_mode(on=F)
}
if (installMetaGx) {
  require(devtools)
  #   dev_mode(on=T)
  install_github("MetaGx", username="bhaibeka", ref="master")
  #   dev_mode(on=F)
}
if (installGenefu) {
  require(devtools)
  #   dev_mode(on=T)
  install_github("genefu", username="bhaibeka", ref="master")
  #   dev_mode(on=F)
}

#-------------------
### requirements ###
#-------------------

## for loading TGGATES data, GSEA, ...
library("PharmacoGx")

#-----------------------------
### clean previous results ###
#-----------------------------

## clean intermediate results
if (cleanIntermediates) {
  unlink(file.path(intermediateDir,"*"))
  message("Cleaned all intermediate results!")
}
## clean gmt files
if (cleanGMTs) {
  unlink(file.path(gseaDir,"*.gmt"))
  message("Cleaned all gmt files!")
}
## clean gsea results
if (cleanGSEAout) {
  unlink(file.path(gseaOutDir,"*"), recursive=T)
  message("Cleaned all GSEA results recursively!")
}
## clean rdata
if (cleanRdata) {
  unlink(file.path(rdataDir,"*"))
  message("Cleaned all rdata files!")
}
## clean ranks
if (cleanRanks) {
  unlink(file.path(rankDir,"*.rnk"))
  message("Cleaned all ranks!")
}
## clean biclustering results
if (cleanBiclusteringSteps) {
  unlink(file.path(rdataDir,"overlaps*"))
  unlink(file.path(rdataDir,"modules*"))
  message("Cleaned all biclustering results!")
}

if (cleanOnly) stop("Only cleaning of files demanded!")

if (!file.exists(gseaExecPath)) stop("Please put gsea file into root dir!")

#----------------
### load data ###
#----------------

#----------------------------
### load common functions ###
#----------------------------

wd_tmp <- getwd()
# setwd("~/Dropbox/current_projects/717genesignatures/R/analysis/scripts")
# source("custom_functions_only.R")
setwd(wd_tmp)

#--------------------------------------
### additional function definitions ###
#--------------------------------------

#######################################################################################################################################

if (verbose) message(sprintf("Calling %s core(s)",ncores))

#######################################################################################################
################################## ###### START PIPELINE #### #########################################
#######################################################################################################

#------------------------------------------------------------------------------------
### Download TGGATES data, normalize, curate, and analyse predictive power of dose 
### response by the use of linear models. Calculate gene ranks according to results
### of linear models
#------------------------------------------------------------------------------------
if (verbose) message("Pipeline A started ...")
sourcef("perform.pipelineA_for_id.R")
## returns rank files (path) and gids
resA.797 <- perform.pipelineA_for_id("797")
resA.798 <- perform.pipelineA_for_id("798")
resA.799 <- perform.pipelineA_for_id("799", tissue="liver")
resA.800 <- perform.pipelineA_for_id("800", tissue="liver")

# adding kidney to compare in response to review EHP (Nov 4, 2014)
resA.799.kidney <- perform.pipelineA_for_id("799", tissue="kidney") 

#-----------------------------------------------------------------------------------------
### Generate gene set files for rat and human. Gene IDs retained from curation are used
#-----------------------------------------------------------------------------------------
if (verbose) message("Pipeline B started ...")
sourcef("perform.pipelineB.R")
## returns gmt file name
resB.commonGmts <- perform.pipelineB(gids.rat=resA.799$gids, gids.human=resA.798$gids, 
                                     genesetSize.min=gseaMinSize, genesetSize.max=gseaMaxSize)

#-----------------------------------------------------------------------------------------
### Perform GSEA on gene ranks
#-----------------------------------------------------------------------------------------
if (verbose) message("Pipeline C started ...")
sourcef("perform.pipelineC_for_id.R")

## extract information from gmt files to put in GSEA
# Reactome
ratReactomeGmt_fn <- resB.commonGmts$reactome$gmt.rat$fn.common.formatted
ratReactomeGsDescription <- resB.commonGmts$reactome$gmt.rat$gmt.common.formatted$description

humanReactomeGmt_fn <- resB.commonGmts$reactome$gmt.human$fn.common.formatted
humanReactomeGsDescription <- resB.commonGmts$reactome$gmt.human$gmt.common.formatted$description

## start GSEA analysis FOR REACTOME
res.gsea.797.reactome <-
  perform.pipelineC_for_id("797", resA.797$rankFiles, min.size=gseaMinSize, max.size=gseaMaxSize,
                           gmt.path=ratReactomeGmt_fn, genesetDescription=ratReactomeGsDescription,
                           out.dir=gseaOutDir_reactome, exec.path=gseaExecPath, nPerm=1000, cores=ncores,
                           lrdataDir=rdataDir, lintermediateDir=intermediateDir)

res.gsea.798.reactome <- 
  perform.pipelineC_for_id("798", resA.798$rankFiles, min.size=gseaMinSize, max.size=gseaMaxSize,
                           gmt.path=humanReactomeGmt_fn, genesetDescription=humanReactomeGsDescription,
                           out.dir=gseaOutDir_reactome, exec.path=gseaExecPath, nPerm=1000, cores=ncores,
                           lrdataDir=rdataDir, lintermediateDir=intermediateDir)

res.gsea.799.reactome <- 
  perform.pipelineC_for_id("799", resA.799$rankFiles, min.size=gseaMinSize, max.size=gseaMaxSize,
                           gmt.path=ratReactomeGmt_fn, genesetDescription=ratReactomeGsDescription,
                           out.dir=gseaOutDir_reactome, exec.path=gseaExecPath, nPerm=1000, cores=ncores,
                           lrdataDir=rdataDir, lintermediateDir=intermediateDir,
                           tissue="liver") 

res.gsea.800.reactome <-
  perform.pipelineC_for_id("800", resA.800$rankFiles, min.size=gseaMinSize, max.size=gseaMaxSize,
                           gmt.path=ratReactomeGmt_fn, genesetDescription=ratReactomeGsDescription,
                           out.dir=gseaOutDir_reactome, exec.path=gseaExecPath, nPerm=1000, cores=ncores,
                           lrdataDir=rdataDir, lintermediateDir=intermediateDir,
                           tissue="liver")

res.gsea.799.kidney.reactome <- 
  perform.pipelineC_for_id("799", resA.799.kidney$rankFiles, min.size=gseaMinSize, max.size=gseaMaxSize,
                           gmt.path=ratReactomeGmt_fn, genesetDescription=ratReactomeGsDescription,
                           out.dir=gseaOutDir_reactome, exec.path=gseaExecPath, nPerm=1000, cores=ncores,
                           lrdataDir=rdataDir, lintermediateDir=intermediateDir,
                           tissue="kidney") 

#-----------------------------------------------------------------------------------------
### Perform Bi-Clustering on enrichment by GSEA
#-----------------------------------------------------------------------------------------
if (verbose) message("Pipeline D started ...")
sourcef("perform.pipelineD_for_id.R")
## cutoff in hypergeometric test
pvalue <- 10^-3

# reactome thrrow <- seq(1.5,2,by=0.5)
thrrow <- list(thrRow797=seq(1.6, 2.6, by=0.2),
               thrRow798=seq(1.6, 2.6, by=0.2),
               thrRow799=seq(2, 3, by=0.2),
               thrRow800=seq(1.6, 2.6, by=0.2))
thrrow <- c(thrrow, list(thrRow799.kidney=thrrow$thrRow799))
gseaResults.reactome <- list(gsea797=res.gsea.797.reactome, gsea798=res.gsea.798.reactome, 
                             gsea799=res.gsea.799.reactome, gsea800=res.gsea.800.reactome,
                             gsea799.kidney=res.gsea.799.kidney.reactome)
overlappingModules.reactome <- perform.pipelineD_for_id(gseaResults.reactome, genesetSource="reactome",
                                                        seed=1, thrRow=thrrow, p.value=pvalue, cor.limit=0.5)

save(overlappingModules.reactome, file="PipelineResults_modulesANDoverlaps.reactome.rda")

# to display result objects of isa function better (does not print out everything when object called)
require(biclust)
detach("package:biclust", unload=T)
require(biclust) 

#-----------------------------------------------------------------------------------------
### CREATING SOME NICE FIGURES
#-----------------------------------------------------------------------------------------
if (verbose) message("Pipeline E started ...")
### Plot conservation of modules across datasets in barplots ###
sourcef("perform.pipelineE.R")

resE.reac <- perform.pipelineE(overlappingModules.reactome, title="Overlap in Reactome pathways", sig.thres=pvalue)

## output module numbers
require(xtable)
print(xtable(resE.reac$p.res.nonredundant$df,
             label="tab:fig2",
             caption="Number of modules for Figure 2 of the manuscript."),
      type="latex",
      file=file.path(latexDir, "resE.reac.overlaps-nonredundant.tex"))

### Plot heatmaps on nesTables ###
sourcef("perform.pipelineF.R")

### For plotting Euler diagrams on the number of genesets in our .GMT files ###

getGenesetNames <- function(genesetSource) {
  fn <- file.path(intermediateDir, sprintf("genesetNames_%s",genesetSource))
  if (!file.exists(fn)) {
    require(GSEABase)
    print(sprintf("Getting unformatted, full (not only common) genesets for %s",genesetSource))
    rat <- getGmt(resB.commonGmts[[genesetSource]]$gmt.rat$fn.gmt)
    human <- getGmt(resB.commonGmts[[genesetSource]]$gmt.human$fn.gmt)
    
    # list of genesets
    rat.gs <- lapply(rat, geneIds)
    names(rat.gs) <- names(rat)
    human.gs <- lapply(human, geneIds)
    names(human.gs) <- names(human)
    
    rat.filtered <- Filter(function(x) length(x) > gseaMinSize & length(x) < gseaMaxSize, rat.gs)
    human.filtered <- Filter(function(x) length(x) > gseaMinSize & length(x) < gseaMaxSize, human.gs)
    
    res <- lapply(list(rat=rat.filtered,
                       human=human.filtered), names)
    save(res, file=fn)
  } else load(fn)
  res
}

### Plot Euler diagrams and heatmaps for each module!!! ###

genesetNames.reactome <- getGenesetNames("reactome")
perform.pipelineF(overlappingModules.reactome, genesetSource="Reactome", genesetNames=genesetNames.reactome)

#-----------------------------------------------------------------------------------------
## Gather final results (figures, tables, gene lists, etc.)
#-----------------------------------------------------------------------------------------
dir.create(file.path("final", "supplementary"), recursive=T, showWarnings=F)  
dir.create("tables", recursive=T, showWarnings=F)  
tabs <- list.files(pattern="*xls")
file.copy(from=tabs, to=sapply(tabs, function(x) file.path("tables",x)), overwrite=T)

## manuscript figures
file.copy(from=file.path(plotDir,"Overlap_in_Reactome_pathwaysnonredundant.pdf"),
          to=file.path("final", "Figure2.pdf"), overwrite=T)
file.copy(from=file.path(plotDir,"heatmap_common_Reactome_non-redundant_RLV2vsPHH15vsPRH10_modulePRH.pdf"),
          to=file.path("final", "Figure3_1.pdf"), overwrite=T)
file.copy(from=file.path(plotDir,"heatmap_common_Reactome_non-redundant_RLV2vsPHH15vsPRH10_modulePHH.pdf"),
          to=file.path("final", "Figure3_2.pdf"), overwrite=T)
file.copy(from=file.path(plotDir,"heatmap_common_Reactome_non-redundant_RLV2vsPHH15vsPRH10_moduleRLV.pdf"),
          to=file.path("final", "Figure3_3.pdf"), overwrite=T)
# file.copy(from=file.path(plotDir,"heatmap_common_Reactome_non-redundant_RLV24vsPHH4_modulePHH.pdf"),
#           to=file.path("final", "Figure3_4.pdf"), overwrite=T)
file.copy(from=file.path(plotDir, "heatmaps", "heatmap_Reactome_modules798_non-redundant_module6.pdf"),
          to=file.path("final", "Figure4_1.pdf"), overwrite=T)
file.copy(from=file.path(plotDir, "heatmaps", "heatmap_Reactome_modules799_non-redundant_module5.pdf"),
          to=file.path("final", "Figure4_2.pdf"), overwrite=T)

# ## manuscript tables
# file.copy(from=file.path("tables", "NonredundantModuleOverlaps_reactome.xls"), to=file.path("final", "Table1.xls"), overwrite=T)

## supplementary venn diagram
file.copy(from=file.path("plots", "euler_Reactome.pdf"), to=file.path("final", "supplementary", "S1.pdf"), overwrite=T)

## supplementary heatmap figures
supplementaryHeatmaps <- list.files(file.path("plots","heatmaps"), 
                                    pattern="heatmap_Reactome_modules[78][90][7890](\\.kidney)*_non-redundant_module.+.pdf",
                                    full.names=T)
zip(file.path("final", "supplementary", "S3.zip"), supplementaryHeatmaps)
system(sprintf("gs -q -sPAPERSIZE=letter -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=%s/modules.pdf %s",
               latexDir, paste(supplementaryHeatmaps, collapse=" ")))
## supplementary special tables
supplementaryTables <- file.path("tables",c("NonredundantModuleOverlaps_reactome.xls",
                                            "NonredundantModule_special_PRH_PHH_reactome.xls",
                                            "NonredundantModule_special_RLV_RLV_reactome.xls",
                                            "NonredundantModule_special_RLVRepeated_RLVRepeated_reactome.xls",
                                            "NonredundantModule_special_RLV_RLV.kidney_reactome.xls"))
zip(file.path("final", "supplementary", "S4.zip"), supplementaryTables)

## supplementary leading edge tables
lE.PRH <- lapply(overlappingModules.reactome$modules.nonredundant$modules797$leading.genes, data.frame)
lE.PRH <- lapply(lE.PRH, function(x) {if (dim(x)[2]>0)colnames(x) <- "leading_genes"; x})
names(lE.PRH) <- paste0("module",1:length(lE.PRH))

lE.PHH <- lapply(overlappingModules.reactome$modules.nonredundant$modules798$leading.genes, data.frame)
lE.PHH <- lapply(lE.PHH, function(x) {if (dim(x)[2]>0)colnames(x) <- "leading_genes"; x})
names(lE.PHH) <- paste0("module",1:length(lE.PHH))

lE.RLV <- lapply(overlappingModules.reactome$modules.nonredundant$modules799$leading.genes, data.frame)
lE.RLV <- lapply(lE.RLV, function(x) {if (dim(x)[2]>0) colnames(x) <- "leading_genes"; x})
names(lE.RLV) <- paste0("module",1:length(lE.RLV))

lE.RLV.repeated <- lapply(overlappingModules.reactome$modules.nonredundant$modules800$leading.genes, data.frame)
lE.RLV.repeated <- lapply(lE.RLV.repeated, function(x) {if (dim(x)[2]>0) colnames(x) <- "leading_genes"; x})
names(lE.RLV.repeated) <- paste0("module",1:length(lE.RLV.repeated))

lE.RLV.kidney <- lapply(overlappingModules.reactome$modules.nonredundant$modules799.kidney$leading.genes, data.frame)
lE.RLV.kidney <- lapply(lE.RLV.kidney, function(x) {if (dim(x)[2]>0) colnames(x) <- "leading_genes"; x})
names(lE.RLV.kidney) <- paste0("module",1:length(lE.RLV.kidney))

WriteXLS("lE.PRH", file.path("tables", "leadingEdge_PRH.xls"))
WriteXLS("lE.PHH", file.path("tables", "leadingEdge_PHH.xls"))
WriteXLS("lE.RLV", file.path("tables", "leadingEdge_RLV.xls"))
WriteXLS("lE.RLV.repeated", file.path("tables", "leadingEdge_RLV.repeated.xls"))
WriteXLS("lE.RLV.kidney", file.path("tables", "leadingEdge_RLV.kidney.xls"))
supplementary_leadingEdges <- file.path("tables", c("leadingEdge_PRH.xls", "leadingEdge_PHH.xls", 
                                                    "leadingEdge_RLV.xls", "leadingEdge_RLV.repeated.xls",
                                                    "leadingEdge_RLV.kidney.xls"))
zip(file.path("final","supplementary","S5.zip"), supplementary_leadingEdges)

supplementary_histograms <- list.files("plots",pattern="wilcox.*", full.names=T)
zip(file.path("final","supplementary","S8.zip"), supplementary_histograms)
system(sprintf("gs -q -sPAPERSIZE=letter -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=%s/histograms.pdf %s",
               latexDir, paste(supplementary_histograms, collapse=" ")))

#-----------------------------------------------------------------------------------------
## save session info
#-----------------------------------------------------------------------------------------

write(toLatex(sessionInfo(), locale = FALSE), file=file.path(latexDir,"sessionInfoR.tex"), append=FALSE)

#-----------------------------------------------------------------------------------------
## copy supplementary latex to final directory
#-----------------------------------------------------------------------------------------

file.copy(from=latexDir, to=file.path("final","supplementary"), overwrite=T, recursive=T)