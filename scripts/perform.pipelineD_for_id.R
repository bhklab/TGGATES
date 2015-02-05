##############################################################
###
### Perform bi-clustering on the enrichment scores of GSEA
### to identify overlapping enrichments
###
### patrick@jimmy.harvard.edu
###
##############################################################

#' Perform bi-clustering on table of NES values
#' @param nesTable table of NES values (rows gene sets, columsn drugs)
#' @param return return bi-clustering results (modules)
biClustering <- function(nesTable, lplotDir, id,
                         seed=987654321,
                         thrRow=seq(1.5,2.5,by=0.5), # seq(1.5,2,by=0.5)
                         thrCol=1, plotHeatmap=F) {
  ## start bi-clustering here ##
  require(isa2)
  require(eisa)
  require(biclust)
  
  require(Biobase)
  
#   mm <- matrix(runif(100),ncol=5)
#   colnames(mm) <- letters[1:ncol(mm)]
#   rownames(mm) <- rev(letters)[1:nrow(mm)]
#   modules <- isa(mm)
#   eset <- ExpressionSet(assayData=mm)
#   bc <-isa.biclust(modules)
#   bc2 <- annotate(bc,eset)
#   modules2 <- as(bc2, "ISAModules")
  
  nesTable.orig <- nesTable
  
#   nesTable <- matrix(runif(length(nesTable.orig)),ncol=ncol(nesTable.orig))
#   rownames(nesTable) <- rownames(nesTable.orig)
#   colnames(nesTable) <- colnames(nesTable.orig)
  
  set.seed(seed)
  ## kind of tricky approach abusing data structures for gene expression
  
  # comment
  tmp <- tryCatch({
    modules.isa <- isa(nesTable, thr.row=thrRow, thr.col=thrCol)
    Bc <- isa.biclust(modules.isa)
    list(modules.isa=modules.isa, Bc=Bc)
  }, error=function(e) {
    modules.isa <- isa(nesTable)
    Bc <- isa.biclust(modules.isa)
    list(modules.isa=modules.isa, Bc=Bc)
  })

  modules.isa <- tmp$modules.isa
  Bc <- tmp$Bc
  
  # comment
  eset <- ExpressionSet(assayData=nesTable)
  Bc <- annotate(Bc, eset) # eisa package
  modules <- as(Bc, "ISAModules")
  
  if (plotHeatmap) {
    pdf(file.path(plotDir,sprintf("isa2heatmap_%s.pdf",id)))
    sapply(1:length(modules), function(x)      
      ISA2heatmap(modules, 1, eset, margins = c(10, 7), cexCol = 0.8, cexRow= 0.7)
    )
    dev.off()
  }
  
  list(modules=modules, modules.isa=modules.isa)
}

#' Get NES table from list of gsea results
#' @param id tggates accession
#' @param gsea.results list of gsea results. gsea.results object
#' outputted from perform.pipelineC_for_id
#' @param genesetNames vector of gene set names. will be used to subset gsea.results 
#' (use common drugs)
getNEStable <- function(id, gsea.results, drugNames, lintermediateDir) {  
  fn.nes <- file.path(lintermediateDir,sprintf("%s_nesTable.rda",id))
  if (!file.exists(fn.nes)) {
    gsea.results <- gsea.results[drugNames]
    #   drugnames <- names(gsea.results)
    
    ## create NES table (rows genesets cols drugs)
    nesTable <- sapply(gsea.results, function(x) {
      orderedRWNs <- rownames(x)[order(rownames(x))] # to bring all in same order
      nes <- x[orderedRWNs,"NES"]
      names(nes) <- orderedRWNs # smarter way to do in R?
      nes
    })
  } else {
    print(sprintf("%s skipped because exists",fn.nes))
    load(fn.mns)
  }  
  nesTable
}

#' put together
getNESandBiclustering <- function(gseaResults, commonDrugNames, id, lintermediateDir, 
                                  lplotDir, seed, thrRow, thrCol, shouldSubset, subsetBy="cancer") {
  nesTable <- getNEStable(id, gseaResults, commonDrugNames, lintermediateDir)
  if (shouldSubset)
    nesTable <- nesTable[grep(subsetBy,rownames(nesTable),ignore.case=T),]
  res.biclustering <- 
    biClustering(nesTable, lplotDir, seed=seed, id=id, 
                 thrRow=thrRow, thrCol=thrCol)
  list(nesTable=nesTable,
       res.biclustering=res.biclustering)
}

#' Compute overlapping subsets of two modules
#' @param modules1 result of biclustering from biClustering() or getNESandBiclustering() 
#' @param modules2 result of..
#' @param name1 name used for rows of dimnames of result matrix for module 1
#' @param name2 ... used for column names ...
#' @return matrix 
computeOverlap <- function(modules1, modules2,
                           name1, name2, pValue.cutoff=10^-6) {
#   geneTable1 <- module1@genes
#   geneTable2 <- module2@genes
  
  modules1.length <- length(modules1)
  modules2.length <- length(modules2)
  
  # note to go on with working: names(modules) is NULL
  overlaps <- sapply(1:modules2.length, function(j) sapply(1:modules1.length, function(i) {
    themodule1 <- modules1[[i]]
    themodule2 <- modules2[[j]]
    
#     tryCatch({
      total <- length(union(rownames(themodule1@genes), rownames(themodule2@genes)))
#     }, error=function(e) {
#      print(e) 
#      browser()
#     })
#     
    
    # get feature names
    tmpf <- function(mods) unique(unlist(getFeatureNames(mods)))
    
    features1 <- tmpf(themodule1)
    features2 <- tmpf(themodule2)
    
    intersLength <- length(intersect(features1,features2))
    l1 <- length(features1)
    l2 <- length(features2)
    
    #minxx <- min(1 - cumsum(dhyper(0:intxx-1, l1, totalx-l1, l2)))
    ifelse(intersLength > 0, 1-phyper(intersLength-1, l1, total-l1, l2), 1)
  }))

########################################################################################################
################ Adding this May 12, 2014. Purpose: Last analysis is analysis of overlap ###############
################ to hepatotoxicants and cancer terms ###################################################

  #' compute also overlap of hepatotoxicants
  computeOverlap.chems <- function(themodules, commonName, chems = c("ethionine", "thioacetamide", "methapyrilene",
                                                                     "coumarin", "ethionamide", "monocrotaline", 
                                                                     "haloperidol", "hexachlorobenzene", "gemfibrozil", 
                                                                     "fenofibrate", "clofibrate", "WY14643", "phenobarbital",
                                                                     "phenytoin", "rifampicin", "ethanol", "carbamazepine", 
                                                                     "sulfasalazine", "griseofulvin", "tamoxifen", "tannicacid",
                                                                     "ethinylestradiol", "diazepam", "acetamidofluorene", 
                                                                     "nitrosodiethylamine")) {
    
    res <- sapply(1:length(themodules), function(i, b) {
      themodule <- themodules[[i]]
      total <- length(rownames(themodule@conditions))
      a <- unique(unlist(eisa::getSampleNames(themodule)))
      intersLength <- length(intersect(a, b))
      l1 <- length(a)
      l2 <- length(b)
      ifelse(intersLength > 0, 1-phyper(intersLength-1, l1, total-l1, l2), 1)
    }, b=intersect(rownames(themodules@conditions), chems))
    tt <- data.frame(res)
    colnames(tt) <- commonName
    tt
  }

  computeOverlap.cancer <- function(themodules) {
    
    cancerterms <- grep("CANCER", rownames(themodules@genes), val=T)
    res <- sapply(1:length(themodules), function(i, b) {
      themodule <- themodules[[i]]
      total <- length(rownames(themodule@genes))
      a <- unique(unlist(eisa::getFeatureNames(themodule)))
      intersLength <- length(intersect(a, b))
      l1 <- length(a)
      l2 <- length(b)
      ifelse(intersLength > 0, 1-phyper(intersLength-1, l1, total-l1, l2), 1)
    }, b=intersect(rownames(themodules@genes), cancerterms))
    data.frame(CancerPathways=res)
  }

  addDimNames <- function(overlapTable, theprefix) {
    if (!is.null(dim(overlaps))) {
#       colnames(overlapTable) <- paste(name2, 1:modules2.length, sep=".") # ask Nehme if he assumes ncol(modules) == length(modules), my example shows ncol() < length()
      rownames(overlapTable) <- paste(theprefix, 1:nrow(overlapTable), sep=".")
    }
    overlapTable
  }
  
  ## hepatoxic
  hepatoxicOverlap1 <- addDimNames(computeOverlap.chems(modules1, commonName="Hepatoxicants"), theprefix=name1)
  hepatoxicOverlap2 <- addDimNames(computeOverlap.chems(modules2, commonName="Hepatoxicants"), theprefix=name2)

  ## non hepatoxic
  nonHepatoxicants <-  c("allylalcohol", "aspirin", "erythromycinethylsuccinate", "disulfiram", "ibuprofen", "nicotinicacid", "cimetidine","lornoxicam",
                         "tolbutamide", "promethazine", "chlorpheniramine", "chlorpropamide") 
  nonHepatoxicOverlap1 <- addDimNames(computeOverlap.chems(modules1, commonName="nonHepatoxicants", chems=nonHepatoxicants), theprefix=name1)
  nonHepatoxicOverlap2 <- addDimNames(computeOverlap.chems(modules2, commonName="nonHepatoxicants", chems=nonHepatoxicants), theprefix=name2)

  ## environmental chemicals
  ppara <- c("benziodarone","benzbromarone","fenofibrate", "ibuprofen","clofibrate","WY14643","gemfibrozil")
  environmentalOverlap1 <- addDimNames(computeOverlap.chems(modules1, commonName="PPARA", chems=ppara), theprefix=name1)
  environmentalOverlap2 <- addDimNames(computeOverlap.chems(modules2, commonName="PPARA", chems=ppara), theprefix=name2)

  ## cancer related pathways
  cancerOverlap1 <- addDimNames(computeOverlap.cancer(modules1), theprefix=name1)
  cancerOverlap2 <- addDimNames(computeOverlap.cancer(modules2), theprefix=name2)

  ## gather results
  tmp.list <- list(hepatoxicOverlap1, hepatoxicOverlap2, nonHepatoxicOverlap1, nonHepatoxicOverlap2, cancerOverlap1, cancerOverlap2,
                   environmentalOverlap1, environmentalOverlap2)
  names(tmp.list) <- c(sprintf("Hepatotoxicants_%s",name1), sprintf("Hepatotoxicants_%s",name2), 
                       sprintf("NonHepatotoxicants_%s",name1), sprintf("NonHepatotoxicants_%s",name2), 
                       sprintf("CancerPathways_%s",name1), sprintf("CancerPathways_%s",name2),
                       sprintf("PPARA_%s",name1), sprintf("PPARA_%s",name2))


#   overlaps.special <- Reduce(cbind, tmp.list)
#   colnames(overlaps.special) <- c(sprintf("Hepatotoxicants_%s",name1), sprintf("Hepatotoxicants_%s",name2), 
#                                   sprintf("CancerPathways_%s",name1), sprintf("CancerPathways_%s",name2))

########################################################################################################

  if (!is.null(dim(overlaps))) {
    colnames(overlaps) <- paste(name2, 1:modules2.length, sep=".") # ask Nehme if he assumes ncol(modules) == length(modules), my example shows ncol() < length()
    rownames(overlaps) <- paste(name1, 1:modules1.length, sep=".")
    significant.rows <- apply(overlaps, 2, function(x) (which (x < pValue.cutoff)))
  } else {
    significant.rows <- sapply(overlaps, function(x) (which (x < pValue.cutoff)))
  }
  
  list(overlaps=overlaps,
       special=tmp.list,
       significant.rows=significant.rows)

}

#' Get for each module it's highest overlapping pendantt
#' @param modules result of bi-clustering
removeRedundantModule <- function(res.isa, nesTable, cor.limit) {
  require(isa2)
  
  isa.norm <- isa.normalize(nesTable, prenormalize=FALSE)
  isa.norm.unique  <- isa.unique(isa.norm, res.isa, method=c("cor"), 
                                 ignore.div=TRUE, cor.limit=cor.limit, neg.cor=TRUE, 
                                 drop.zero=TRUE)
  
  Bc <- isa.biclust(isa.norm.unique)
  fakeEset <- ExpressionSet(assayData=nesTable)
  Bc <- annotate(Bc, fakeEset)
  res.modules <- as(Bc, "ISAModules")
  
  res.modules
}

#' 
getModuleLeadingEdge <- function(mod, leading, species) {
  if (species=="rat") {
    require(org.Rn.eg.db)
    db <- 'org.Rn.eg'
  } else if (species=="human") {
    db <- 'org.Hs.eg.db'
    require(org.Hs.eg.db)
  } else stop ("species not implemented")
  require(annotate)
  require(eisa)
  
  drugs <- eisa::getSampleNames(mod)[[1]]
  pathways <- eisa::getFeatureNames(mod)[[1]]
  
  extractLeading <- function(leading.=leading, drugs.=drugs, pathways.=pathways) {
    Reduce(union, lapply(leading.[drugs.], function(x) Reduce(union, x[!sapply(x,is.null)][pathways.])))
  }
  
  rr <- extractLeading()
  
  if (!is.null(rr)) {
    symbx <- lookUp(as.character(rr), db, 'SYMBOL')
    symbx <- unname(symbx, force = FALSE)
    return(toupper(symbx))
  } else return(NULL)
}

addModuleLeadingEdge <- function(modObj, leading, isNonRedundant, species) {
  if (!isNonRedundant)
    mods <- modObj$res.biclustering$modules
  else
    mods <- modObj
  res <- lapply(1:length(mods), function(i) {
    genes <- getModuleLeadingEdge(mod=mods[[i]], leading=leading, species=species)
    genes
  })
  if (!isNonRedundant)
    return(c(modObj, leading.genes=list(res)))
  else
    return(list(modules=modObj, leading.genes=res))
}

#' @param gsea.results list of list of gsea results. gsea.results object
#' first level of the list are the accessions, the second level corresponds
#' to drugs, objects outputted from perform.pipelineC_for_id
#' @param genesetSource string either go or reactome, etc. used for naming
#' intermediate results (special cases are reactome_cancer and gprofiler_cancer!)
#' @return returns bi-clustering results
#' @param thrRow list of isa row thresholds (names refer to TGGATES accessions)
perform.pipelineD_for_id <- function(gsea.results, 
                                     genesetSource,
                                     p.value, cor.limit,
                                     thrRow, seed=987654321,
                                     thrCol=1, lrdataDir=rdataDir,
                                     lintermediateDir=intermediateDir,
                                     lplotDir=plotDir,
                                     ldataDir=dataDir) {
  
  #############################################################
  ### use gprofiler to get list of reactome cancer pathways ###
  
  if (genesetSource=="gprofiler_cancer") {
    fn.cancerGenes <- file.path(ldataDir,"cancer_gene_census.tsv")
    if (!file.exists(fn.cancerGenes)) {
      status <- download.file("ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/data_export/cancer_gene_census.tsv",
                              destfile=fn.cancerGenes)
      if (status != 0) stop("Downloading of ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/data_export/cancer_gene_census.tsv failed, re-run!")
    } 
    cancerGenes.table <- read.csv(fn.cancerGenes, sep="\t")

    fn.gprofilerOut <- file.path(lrdataDir,"gprofiler_cancer.rda")
    if (!file.exists(fn.gprofilerOut)) {
      require(gProfileR)
      genes <- as.vector(unique(cancerGenes.table$Symbol))
      gprofiler.human <- gprofiler(genes, organism="hsapiens")
      gprofiler.rat <- gprofiler(genes, organism="rnorvegicus")
      save(list=c("gprofiler.human","gprofiler.rat"), file=fn.gprofilerOut)
    } else load(fn.gprofilerOut)
    
    a <- unique(gprofiler.human[gprofiler.human$p.value<0.05,])
    a <- a[grep("REAC",a$term.id),"term.name"]
    b <- unique(gprofiler.rat[gprofiler.rat$p.value<0.05,])
    b <- b[grep("REAC",b$term.id),"term.name"]
    gprofiler.intersectGenesets <- intersect(a,b)
    write(gprofiler.intersectGenesets, file="gprofiler.intersectGenesets_reactome.txt")
  }
  
  ############################################################
  ### Get bi-clustering by NES scores of overlapping drugs ###
  
  res.gsea.results.797 <- gsea.results$gsea797$res.gseas
  res.gsea.results.798 <- gsea.results$gsea798$res.gseas
  res.gsea.results.799 <- gsea.results$gsea799$res.gseas
  res.gsea.results.800 <- gsea.results$gsea800$res.gseas
  
  ## leading edge genes
  res.gsea.results.797.leading <- gsea.results$gsea797$res.leadingEdges
  res.gsea.results.798.leading <- gsea.results$gsea798$res.leadingEdges
  res.gsea.results.799.leading <- gsea.results$gsea799$res.leadingEdges
  res.gsea.results.800.leading <- gsea.results$gsea800$res.leadingEdges
  
  ## add 799 kidney
  res.gsea.results.799.kidney <- gsea.results$gsea799.kidney$res.gseas
  res.gsea.results.799.kidney.leading <- gsea.results$gsea799.kidney$res.leadingEdges
  
  ## parameters
  thrRow.797 <- thrrow$thrRow797
  thrRow.798 <- thrrow$thrRow798
  thrRow.799 <- thrrow$thrRow799
  thrRow.800 <- thrrow$thrRow800
  
  ## add  799 kidney
  thrRow.799.kidney <- thrrow$thrRow799.kidney
  
  ## take overlapping drugs in 797, 798, 799, and 800. exclude 799 kidney
  iix <- !grepl("799.kidney", names(gsea.results))
  commonDrugs <- Reduce(intersect, sapply(sapply(gsea.results[iix], "[[", "res.gseas"), names))
  
  fn.modules <- file.path(lrdataDir,sprintf("modules_%s.rda",genesetSource))
#   if (!file.exists(fn.modules)) {
    print("Performing bi-clustering")
    shouldSubset <- ifelse(genesetSource=="reactome_cancer",TRUE,FALSE)
    subsetBy <- ifelse(genesetSource=="reactome_cancer","cancer","") # TODO CHECK IF "" IS CORRECT BEHAVIOR
    # contain nesTable and results of biclustering (modules)
    modules797 <- getNESandBiclustering(res.gsea.results.797, commonDrugs, id="797", 
                                        lintermediateDir, lplotDir, seed=seed,
                                        shouldSubset=shouldSubset,
                                        subsetBy=subsetBy,
                                        thrRow=thrRow.797,
                                        thrCol=thrCol) # rat in vitro
    modules798 <- getNESandBiclustering(res.gsea.results.798, commonDrugs, id="798", 
                                        lintermediateDir, lplotDir, seed=seed,
                                        shouldSubset=shouldSubset,
                                        subsetBy=subsetBy,
                                        thrRow=thrRow.798,
                                        thrCol=thrCol) # human in vitro
    modules799 <- getNESandBiclustering(res.gsea.results.799, commonDrugs, id="799",
                                        lintermediateDir, lplotDir, seed=seed,
                                        shouldSubset=shouldSubset,
                                        subsetBy=subsetBy,
                                        thrRow=thrRow.799,
                                        thrCol=thrCol) # rat in vivo
    modules800 <- getNESandBiclustering(res.gsea.results.800, commonDrugs, id="800",
                                        lintermediateDir, lplotDir, seed=seed,
                                        shouldSubset=shouldSubset,
                                        subsetBy=subsetBy,
                                        thrRow=thrRow.800,
                                        thrCol=thrCol) # rat in vivo repeated dose
    modules799.kidney <- getNESandBiclustering(res.gsea.results.799.kidney, 
                                               intersect(commonDrugs,names(res.gsea.results.799.kidney)), id="799.kidney",
                                               lintermediateDir, lplotDir, seed=seed,
                                               shouldSubset=shouldSubset,
                                               subsetBy=subsetBy,
                                               thrRow=thrRow.799.kidney,
                                               thrCol=thrCol) # rat in vivo kidney
    modules797 <- addModuleLeadingEdge(modules797,res.gsea.results.797.leading, isNonRedundant=F, species="rat")
    modules798 <- addModuleLeadingEdge(modules798,res.gsea.results.798.leading, isNonRedundant=F, species="human")
    modules799 <- addModuleLeadingEdge(modules799,res.gsea.results.799.leading, isNonRedundant=F, species="rat")
    modules800 <- addModuleLeadingEdge(modules800,res.gsea.results.800.leading, isNonRedundant=F, species="rat")
    modules799.kidney <- addModuleLeadingEdge(modules799.kidney,res.gsea.results.799.kidney.leading, isNonRedundant=F, species="rat")

    save(list=c("modules797","modules798","modules799","modules800", "modules799.kidney"), file=fn.modules)
#   } else {
#     message("Computing bi-clustering skipped because exists")
#     load(fn.modules)
#   }
  
  ###############################
  ### Get overlapping modules ###
  
  ## get overlapping modules for rat in-vivo
  
  fn.overlaps <- file.path(lrdataDir,sprintf("overlaps_%s.rda",genesetSource))
#   if (!file.exists(fn.overlaps)) {
    print("Computing signifance of module overlaps")
    overlaps.RLV.RLV <-
      computeOverlap(modules799$res.biclustering$modules, modules799$res.biclustering$modules, "RLV1", "RLV2", pValue.cutoff=p.value)
    overlaps.PRH.PHH <- 
      computeOverlap(modules797$res.biclustering$modules, modules798$res.biclustering$modules, "PRH", "PHH", pValue.cutoff=p.value)
    overlaps.RLV.PHH <- 
      computeOverlap(modules799$res.biclustering$modules, modules798$res.biclustering$modules, "RLV", "PHH", pValue.cutoff=p.value)
    overlaps.RLV.PRH <- 
      computeOverlap(modules799$res.biclustering$modules, modules797$res.biclustering$modules, "RLV", "PRH", pValue.cutoff=p.value)  
    ## repeated dose
    overlaps.RLVRepeated.RLVRepeated <- 
      computeOverlap(modules800$res.biclustering$modules, modules800$res.biclustering$modules, "RLV.repeated1", "RLV.repeated2", pValue.cutoff=p.value)    
    overlaps.RLVRepeated.PHH <- 
      computeOverlap(modules800$res.biclustering$modules, modules798$res.biclustering$modules, "RLV.repeated", "PHH", pValue.cutoff=p.value)
    overlaps.RLVRepeated.PRH <- 
      computeOverlap(modules800$res.biclustering$modules, modules797$res.biclustering$modules, "RLV.repeated", "PRH", pValue.cutoff=p.value)  
    ## kidney 
    overlaps.RLV.RLV.kidney <- 
      computeOverlap(modules799$res.biclustering$modules, modules799.kidney$res.biclustering$modules, "RLV", "RLV.kidney", pValue.cutoff=p.value)    
    
    save(list=c("overlaps.RLV.RLV","overlaps.PRH.PHH",
                "overlaps.RLV.PRH","overlaps.RLV.PHH",
                "overlaps.RLVRepeated.RLVRepeated",
                "overlaps.RLVRepeated.PHH",
                "overlaps.RLVRepeated.PRH",
                "overlaps.RLV.RLV.kidney"),
         file=fn.overlaps)
#   } else {
#     message("Overlaps skipped because existing")
#     load(fn.overlaps)
#   }

  ################################
  ### Remove redundant modules ###

  fn.modules.nonredundant <- file.path(lrdataDir,sprintf("modules.nonredundant_%s",genesetSource))
  print("Remove redundant modules")
  modules797.nonredundant <- removeRedundantModule(modules797$res.biclustering$modules.isa, modules797$nesTable, cor.limit=cor.limit)
  modules798.nonredundant <- removeRedundantModule(modules798$res.biclustering$modules.isa, modules798$nesTable, cor.limit=cor.limit)
  modules799.nonredundant <- removeRedundantModule(modules799$res.biclustering$modules.isa, modules799$nesTable, cor.limit=cor.limit)
  modules800.nonredundant <- removeRedundantModule(modules800$res.biclustering$modules.isa, modules800$nesTable, cor.limit=cor.limit)
  modules799.kidney.nonredundant <- removeRedundantModule(modules799.kidney$res.biclustering$modules.isa, modules799.kidney$nesTable, cor.limit=cor.limit)

  modules797.nonredundant <- addModuleLeadingEdge(modules797.nonredundant,res.gsea.results.797.leading, isNonRedundant=T, species="rat")
  modules798.nonredundant <- addModuleLeadingEdge(modules798.nonredundant,res.gsea.results.798.leading, isNonRedundant=T, species="human")
  modules799.nonredundant <- addModuleLeadingEdge(modules799.nonredundant,res.gsea.results.799.leading, isNonRedundant=T, species="rat")
  modules800.nonredundant <- addModuleLeadingEdge(modules800.nonredundant,res.gsea.results.800.leading, isNonRedundant=T, species="rat")
  modules799.kidney.nonredundant <- addModuleLeadingEdge(modules799.kidney.nonredundant,res.gsea.results.799.leading, isNonRedundant=T, species="rat")

  ############################################
  ### Compute overlap in redundant modules ###

  fn.overlaps <- file.path(lrdataDir,sprintf("overlaps.nonredundant_%s.rda",genesetSource))
  #   if (!file.exists(fn.overlaps)) {
  print("Computing signifance of nonredundant module overlaps")
  overlaps.RLV.RLV.nonredundant <-
    computeOverlap(modules799.nonredundant$modules, modules799.nonredundant$modules, "RLV1", "RLV2", pValue.cutoff=p.value)
  overlaps.PRH.PHH.nonredundant <- 
    computeOverlap(modules797.nonredundant$modules, modules798.nonredundant$modules, "PRH", "PHH", pValue.cutoff=p.value)
  overlaps.RLV.PHH.nonredundant <- 
    computeOverlap(modules799.nonredundant$modules, modules798.nonredundant$modules, "RLV", "PHH", pValue.cutoff=p.value)
  overlaps.RLV.PRH.nonredundant <- 
    computeOverlap(modules799.nonredundant$modules,modules797.nonredundant$modules, "RLV", "PRH", pValue.cutoff=p.value)  
  ## repeated dose
  overlaps.RLVRepeated.RLVRepeated.nonredundant <- 
    computeOverlap(modules800.nonredundant$modules,modules800.nonredundant$modules, "RLV.repeated1", "RLV.repeated2", pValue.cutoff=p.value) 
  overlaps.RLVRepeated.PHH.nonredundant <- 
    computeOverlap(modules800.nonredundant$modules, modules798.nonredundant$modules, "RLV.repeated", "PHH", pValue.cutoff=p.value)
  overlaps.RLVRepeated.PRH.nonredundant <- 
    computeOverlap(modules800.nonredundant$modules,modules797.nonredundant$modules, "RLV.repeated", "PRH", pValue.cutoff=p.value)  
  ## kidney
  overlaps.RLV.RLV.kidney.nonredundant <- 
    computeOverlap(modules799.kidney.nonredundant$modules,modules799.kidney.nonredundant$modules, "RLV", "RLV.kidney", pValue.cutoff=p.value)  

  ## print XLS sheets ##

  ## overlap of modules
  overlaps.RLVvsRLV <- as.data.frame(overlaps.RLV.RLV.nonredundant$overlaps)
#   overlaps.RLVvsRLV[is.na(overlaps.RLVvsRLV)] <- 1
  overlaps.PRHvsPHH <- as.data.frame(overlaps.PRH.PHH.nonredundant$overlaps)
#   overlaps.PRHvsPHH[is.na(overlaps.PRHvsPHH)] <- 1
  overlaps.RLVvsPRH <- as.data.frame(overlaps.RLV.PRH.nonredundant$overlaps)
#   overlaps.RLVvsPRH[is.na(overlaps.RLVvsPRH)] <- 1
  overlaps.RLVvsPHH <- as.data.frame(overlaps.RLV.PHH.nonredundant$overlaps)
#   overlaps.RLVvsPHH[is.na(overlaps.RLVvsPHH)] <- 1
  ## repeated dose
  overlaps.RLVRepeatevsRLVRepeate <- as.data.frame(overlaps.RLVRepeated.RLVRepeated.nonredundant$overlaps)
  overlaps.RLVRepeatedvsPHH <- as.data.frame(overlaps.RLVRepeated.PHH.nonredundant$overlaps)
  overlaps.RLVRepeatedvsPRH <- as.data.frame(overlaps.RLVRepeated.PRH.nonredundant$overlaps)
  ## kidney
  overlaps.RLVvsRLV.kidney <- as.data.frame(overlaps.RLV.RLV.kidney.nonredundant$overlaps)

  require(WriteXLS)
  WriteXLS(c("overlaps.RLVvsRLV", "overlaps.PRHvsPHH", "overlaps.RLVvsPRH", "overlaps.RLVvsPHH", 
             "overlaps.RLVRepeatevsRLVRepeate", "overlaps.RLVRepeatedvsPHH",
             "overlaps.RLVRepeatedvsPRH", "overlaps.RLVvsRLV.kidney"),
           sprintf("NonredundantModuleOverlaps_%s.xls",genesetSource),
           row.names=TRUE)

  printXLS4Special <- function(listOfDfs, fn) {
    env <- environment()
    ## if names are longer than 31 characters (won't work with the xls .pl script)
    names(listOfDfs) <- sapply(names(listOfDfs), function(x) {
      splittedName <- strsplit(x, "")[[1]]
      if (length(splittedName) > 31) {
        nn <- length(splittedName) 
        thename <- paste0(splittedName[(nn-30):nn], collapse="")
      } else thename <- x
      thename
    })
    ## assign names to local workspace and print xls
    for (dfname in names(listOfDfs)) {
      assign(dfname, listOfDfs[[dfname]], envir=env)
    }
    WriteXLS(names(listOfDfs), fn, row.names=TRUE)
  }

  printXLS4Special(overlaps.RLV.RLV.nonredundant$special,
                   fn=sprintf("NonredundantModule_special_RLV_RLV_%s.xls",genesetSource))
  printXLS4Special(overlaps.PRH.PHH.nonredundant$special,
                   fn=sprintf("NonredundantModule_special_PRH_PHH_%s.xls",genesetSource))
  printXLS4Special(overlaps.RLV.PRH.nonredundant$special,
                   fn=sprintf("NonredundantModule_special_RLV_PRH_%s.xls",genesetSource))
  printXLS4Special(overlaps.RLV.PHH.nonredundant$special,
                   fn=sprintf("NonredundantModule_special_RLV_PHH_%s.xls",genesetSource))
  ## repeated dose
  printXLS4Special(overlaps.RLVRepeated.RLVRepeated.nonredundant$special,
                  fn=sprintf("NonredundantModule_special_RLVRepeated_RLVRepeated_%s.xls",genesetSource))
  printXLS4Special(overlaps.RLVRepeated.PRH.nonredundant$special,
                   fn=sprintf("NonredundantModule_special_RLVRepeated_PRH_%s.xls",genesetSource))
  printXLS4Special(overlaps.RLVRepeated.PHH.nonredundant$special,
                  fn=sprintf("NonredundantModule_special_RLVRepeated_PHH_%s.xls",genesetSource))
  ## kidney
  printXLS4Special(overlaps.RLV.RLV.kidney.nonredundant$special,
                   fn=sprintf("NonredundantModule_special_RLV_RLV.kidney_%s.xls",genesetSource))
    
## below code could be used if number of modules were the same (Reduce((rbind, tmp.list in computeOverlap)))
#   overlaps.RLVvsRLV.special <- as.data.frame(overlaps.RLV.RLV.nonredundant$special)
#   overlaps.PRHvsPHH.special <- as.data.frame(overlaps.PRH.PHH.nonredundant$special)
#   overlaps.RLVvsPRH.special <- as.data.frame(overlaps.RLV.PRH.nonredundant$special)
#   overlaps.RLVvsPHH.special <- as.data.frame(overlaps.RLV.PHH.nonredundant$special)
# 
#   WriteXLS(c("overlaps.RLVvsRLV.special", "overlaps.PRHvsPHH.special", "overlaps.RLVvsPRH.special", "overlaps.RLVvsPHH.special"),
#            sprintf("NonredundantModuleOverlaps_special_%s.xls",genesetSource))

  ## save results
  save(list=c("overlaps.RLV.RLV.nonredundant","overlaps.PRH.PHH.nonredundant",
              "overlaps.RLV.PRH.nonredundant","overlaps.RLV.PHH.nonredundant",
              "overlaps.RLVRepeated.RLVRepeated.nonredundant",
              "overlaps.RLVRepeated.PRH.nonredundant", "overlaps.RLVRepeated.PHH.nonredundant",
              "overlaps.RLVvsRLV.kidney"),
       file=fn.overlaps)

  ###############################
  ### Gather results together ###

  modules <- list(modules797=modules797,
                  modules798=modules798,
                  modules799=modules799,
                  modules800=modules800,
                  modules799.kidney=modules799.kidney)

  modules.nonredundant <- list(modules797=modules797.nonredundant,
                               modules798=modules798.nonredundant,
                               modules799=modules799.nonredundant,
                               modules800=modules800.nonredundant,
                               modules799.kidney=modules799.kidney.nonredundant)
  
  overlaps <- list(overlaps.RLV.RLV=overlaps.RLV.RLV,
                   overlaps.PRH.PHH=overlaps.PRH.PHH,
                   overlaps.RLV.PRH=overlaps.RLV.PRH,
                   overlaps.RLV.PHH=overlaps.RLV.PHH,
                   overlaps.RLVRepeated.RLVRepeated=overlaps.RLVRepeated.RLVRepeated,
                   overlaps.RLVRepeated.PRH=overlaps.RLVRepeated.PRH,
                   overlaps.RLVRepeated.PHH=overlaps.RLVRepeated.PHH,
                   overlaps.RLV.RLV.kidney=overlaps.RLV.RLV.kidney)

  overlaps.nonredundant <- list(overlaps.RLV.RLV=overlaps.RLV.RLV.nonredundant,
                                overlaps.PRH.PHH=overlaps.PRH.PHH.nonredundant,
                                overlaps.RLV.PRH=overlaps.RLV.PRH.nonredundant,
                                overlaps.RLV.PHH=overlaps.RLV.PHH.nonredundant,
                                overlaps.RLVRepeated.RLVRepeated=overlaps.RLVRepeated.RLVRepeated.nonredundant,
                                overlaps.RLVRepeated.PRH=overlaps.RLVRepeated.PRH.nonredundant,
                                overlaps.RLVRepeated.PHH=overlaps.RLVRepeated.PHH.nonredundant,
                                overlaps.RLV.RLV.kidney=overlaps.RLV.RLV.kidney.nonredundant)
  
  list(modules=modules,
       overlaps=overlaps,
       modules.nonredundant=modules.nonredundant,
       overlaps.nonredundant=overlaps.nonredundant)
  
}
