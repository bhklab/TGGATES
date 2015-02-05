#################################################################################################################
###
### These are the main custom functions used in the pipeline. The pipeline is generic for TGGATES accession
### number. The pipeline is written as a function that expects this number. This function uses global variables
### given in the call script (Pipeline_TGGATES.R). Intermediate results are stored to rdata/intermediate.
###
### patrick@jimmy.harvard.edu, Mar 31, 2014
###
#################################################################################################################

#' Keep only drugs with n experiments
#' @param eset expression set returned from normalize.TGGATES
#' @param id TGGATES accession 
#' @return subsetted expression set with only drugs that have n experiments and list of these
#' objects in case data set has to be split in multiple sets (different number of experiments
#' for same id as in 798 human, also see runGLM)
keepDrugs <- function(eset, id, tissue) {
  # number of experiments require for OUR analysis. USER INPUT!!
#   n <- ifelse(id=="799", 48, 24)  
  
  nvec <- switch(id,
                 "797"=24,
                 "798"=c(12,24),
                 "799"=48,
                 "800"=48,
                 "799.kidney"=48) # 48 is in this implementation only a dummy value (could be used though, change down)
  
  res <- lapply(nvec, function(n) {
    pd <- pData(eset)
    # only liver in case of 799 (rat in vivo also measured in kidney)
#     if (id=="799") pd <- pd[pd$Factor.Value.OrganismPart.=="liver",] 
#     if (id=="799.kidney") pd <- pd[pd$Factor.Value.OrganismPart.=="kidney",] 
    # vector named with drugs holding number of entries for each drug
    inn <- tapply(pd$Parameter.Value.Compound., pd$Parameter.Value.Compound., length) # equiv to table()
    if (tissue=="kidney") 
      inn <- names(which(inn >= 5)) # all that are < 5 have 1 instance only
    else if (id!="800")
      # get names with only drugs with n entries. critical!
      inn <- names(which(inn==n))
    else
      inn <- names(inn) # we consider 800 as validation set, and take ALL instances!
    pd2 <- pd[pd$Parameter.Value.Compound.%in%inn,]
    exp <- exprs(eset)
    exp2 <- exp[,rownames(pd2)]
    varmetadata <- data.frame(labelDescription=colnames(pd2), row.names=colnames(pd2))
    return(new("ExpressionSet", exprs=exp2, 
               phenoData=new("AnnotatedDataFrame", data=pd2, varMetadata=varmetadata), 
               annotation=annotation(eset)))  
  })    
  names(res) <- nvec
  res
}

# splitData <- function(eset) {
#   geneex <- exprs(eset)
#   phenod <- pData(eset)
#   mtable <- cbind(t(geneex)
# }

#' Run linear model on list of data.frames
#' @param dataframes list of data.frames, where each data.frame corresponds to one drug
#' @return list of summaries of glm
runGLM_old <- function(dataframes, family="gaussian", formula=X ~ D + T + T * D) {
  drugnames <- unique(pdata(EMTAB.filtered)$Parameter.Value.Compound.)
  genenames <- rownames(exprs(EMTAB.filtered))
  
  fun.glm <- function(df, family=family) {
    express.m1 <- d[,6:ncol(d)]
    #     group <- df$groupf
    dose <- df$Parameter.Value.Dose.
    time <- df$Parameter.Value.TimeOfCollection.
    df2 <- data.frame("D"=dose, "T"=time, "X"=express.m1)
    apply(df2[,3:ncol(df2)], function(x) summary(glm(formula, data=df2, family=family)))
  }
}

#' Run linear model on eSet for each drug and gene combination to test regression of gene expression 
#' to dose response
#' @param leset list of expression set
#' @return named list of list (drugs times genes) of glm coefficients vector 
#' (merged if input was >1-dimensional,
#' in case of 798 see keepDrugs)
runGLM <- function(leset, id, family="gaussian", formula=G ~ D + T + T * D, cores, intermediateDir, tissue) {
  require(parallel)
  
  ## do for each eSet obtained from keepDrugs (atm only relevant for 798 (human))
  esetRes <- lapply(leset, function(eset) {
    time.startGlm <- proc.time()[["elapsed"]]  
    
    pheno <- pData(eset)
    geneex <- exprs(eset)
    drugnames <- unique(pheno$Parameter.Value.Compound.)
    genenames <- rownames(geneex)
    
    drugnames <- drugnames
    genenames <- genenames
    
    # time parameter predicting gene expression
    timeParam <- ifelse(id=="799" | id=="800", "Parameter.Value.TimeOfSacrifice.", "Parameter.Value.TimeOfCollection.")
    # result list
    drugs.glm <- vector(mode = "list", length = length(drugnames))
    names(drugs.glm) <- drugnames
    ## run glm for each drug and each gene (in parallel)
    ## store intermediate results in intermediate dir
    for (d in drugnames) {
      # store intermediate results
      if (tissue=="")
        fn.d <- file.path(intermediateDir,sprintf("%s_%s_glm_res.d.rda",id,d))
      else
        fn.d <- file.path(intermediateDir,sprintf("%s.%s_%s_glm_res.d.rda",id,tissue,d))
      if (!file.exists(fn.d)) {
        res.d <-
          mclapply(genenames, function(g) { # parallel for each drug
            # rows with drug d
            # as.character required because rownames(799 dataset) are integers..
            insamples <- as.character(rownames(pheno)[which(pheno$Parameter.Value.Compound.==d)])
            # create minimal data.frame; _!TODO!_ could be more efficient to outsource in 
            # huge df first and subset from here on 
            df <- data.frame("G"=geneex[g,insamples],
                             "D"=pheno[insamples,"Parameter.Value.Dose."],
                             "T"=pheno[insamples,timeParam])
            # fit linear model
            tmp <- summary(glm(formula, data=df, family=family))
            # get results
            tmp$coefficients["D", c("Estimate", "Pr(>|t|)")]
          }, mc.cores=ncores)
        names(res.d) <- genenames
        save(res.d, file=fn.d)
      } else load(fn.d)
      # assign glm results
      drugs.glm[[d]] <- res.d
    } ## end drug loop
    
    time.endGlm <- proc.time()[["elapsed"]]
    print(sprintf("Took %s seconds to compute n_drugXn_genes GLMs",time.endGlm-time.startGlm))
    
    return(drugs.glm)
  })
  
  ## prepare results for processCoeff function
  
  ## concatenate first elements of esetRes
  res <- unlist(esetRes, recursive=F)
  tmp_names <- unlist(sapply(esetRes, function(drugs) names(drugs), simplify=F)) # names are drugs  
  stopifnot(length(tmp_names) != unique(tmp_names)) # if there are duplicates
  names(res) <- tmp_names
  
  res  
}

#' Process coefficients of glm summary results (merge results stored in list to big df and compute
#' weights of genes for rank)
#' @param res.glm list of coefficients of glm results (for each drug and gene) obtained from runGLM()
#' @return list of only coefficients of glm results
processCoeff <- function(res.glm, cores) {
#   res <- lapply(res.glm[1:2], function(d) {   
#     mclapply(d, function(g) g$coefficients["D", c("Estimate", "Pr(>|t|)")], mc.cores=cores)
#   })

  res <- res.glm
  
  # merge results
  res.merge <- do.call(rbind, mclapply(res, data.frame, stringsAsFactors=FALSE, mc.cores=cores))
  res.merge <- t(res.merge)
#   rownames(res.merge) <- gsub("X.", "", rownames(res.merge)) # Nehme?
  
  ## relevant columns
  justestim <- grep("Estimate", colnames(res.merge))
  justpr <- grep("Pr(>|t|)", colnames(res.merge))
  
  #### build 2 dataframes: 1 for the estimate and 1 for p values and assign drugnames for columns
  df.estimate <- res.merge[, justestim ]
  df.pvalue <- res.merge[,justpr]
  colnames(df.estimate) <- names(res)
  colnames(df.pvalue) <- names(res)
  
  ### applying the weighting formula on a data frame (sign of estimate * p value)
  signes <- sign(df.estimate)
  prob <- -log10(df.pvalue)
  df.weight <- signes * prob
  
  ## return ordered data.frame
  df.weight[,order(colnames(df.weight))]  
}

#' Compute gene ranks (for each drug) according to processed results of linear models
#' @param res.glm.coeff data.frame as from getCoeff()
#' @param rankDir directory where to store ranks
#' @param cores number of cores
#' @return named vector of rank file names (path), names are drugnames
makeRanks <- function(res.glm.coeff, id, rankDir, cores, tissue) {
  require(parallel)
  
  colnames(res.glm.coeff) <- gsub("[ -]", "", colnames(res.glm.coeff)) # replace problematic DRUG names for GSEA
  res <- mclapply(colnames(res.glm.coeff), function(d) {   
    if (tissue=="")
      fff <- file.path(rankDir, sprintf("%s_%s.rnk", id, d))
    else
      fff <- file.path(rankDir, sprintf("%s.%s_%s.rnk", id, tissue, d))
    if (!file.exists(fff)) {
      ss <- res.glm.coeff[,d]
      ss <- sort(ss, decreasing=TRUE, na.last=NA)
      ss[ss == Inf] <- .Machine$double.xmax
      ss[ss == -Inf] <- -.Machine$double.xmax
      # delete geneid. prepend for correct rnk format for GSEA
      names(ss) <-gsub("geneid\\.(.+)", "\\1", names(ss))
      rankg<-cbind(names(ss), ss)      
      write.table(rankg, file=fff, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
    } else message(sprintf("Skipping %s because exists",fff))
    fff
  }, mc.cores=cores)
  names(res) <- colnames(res.glm.coeff)
  unlist(res) # return vector of filenames
}

compare.GLMgenes <- function(res.glm, fn,
                              hepatotoxicants=c("ethionine", "thioacetamide", "methapyrilene",
                                                "coumarin", "ethionamide", "monocrotaline",
                                                "haloperidol", "hexachlorobenzene", "gemfibrozil",
                                                "fenofibrate", "clofibrate", "WY-14643", "phenobarbital",
                                                "phenytoin", "rifampicin", "ethanol", "carbamazepine",
                                                "sulfasalazine", "griseofulvin", "tamoxifen", "tannic acid",
                                                "ethinylestradiol", "diazepam", "acetamidofluorene",
                                                "nitrosodiethylamine"),
                              fdr=0.1, title="") {
  
  forOneChem <- function(chem) {
    pp <- sapply(chem, function(x) x["Pr(>|t|)"]) # p-values
    qq <- p.adjust(pp, method="fdr")
    sum(pp <= fdr)
  }
  
  numberSigGenes <- sapply(res.glm, forOneChem)
  numberSigGenes.hepato <- numberSigGenes[hepatotoxicants%in%names(numberSigGenes)]
  numberSigGenes.nonHepato <- numberSigGenes[!names(numberSigGenes)%in%hepatotoxicants]
  res.wilcox <- wilcox.test(numberSigGenes.hepato, numberSigGenes.nonHepato)
  
  pdf(fn)
  hist(numberSigGenes.hepato, breaks=20, type="l", xlim=range(c(numberSigGenes.hepato,numberSigGenes.nonHepato)),
       col="red", xlab="Number of significant genes", probability=TRUE, freq=FALSE, main=title)
  hist(numberSigGenes.nonHepato, breaks=20, col=rgb(0, 0, 1, 0.5), add=T, freq=FALSE)
  legend("top", legend=c("Hepatotoxic", "Not hepatotoxic"), fil=c("red","blue"))
  rp <- vector('expression',1)
  rp[1] <- substitute(expression(italic(p) == MYOTHERVALUE), 
                      list(MYOTHERVALUE = format(res.wilcox$p.value, digits = 2)))[2]
  legend("topright", legend=rp, bty='n')
  dev.off()
  list(numberSigGenes.hepato=numberSigGenes.hepato,
       numberSigGenes.nonHepato=numberSigGenes.nonHepato)
}

#' Performs full analysis pipeline for one data set. Includes downloading,
#' normalizing, curating, subsetting, running glm, pre-ranking genes.
#' @param id array express ID either 797, 798, or 799
#' @param rdataDir directory to store rdata
#' @param intermediateDir directory to store intermediate results
#' @return list of rank file names (full path), and gids after curation
perform.pipelineA_for_id <- function(id, lrdataDir=rdataDir, lintermediateDir=intermediateDir, 
                                    lrankDir=rankDir, lgseaDir=gseaDir, cores=ncores, tissue="") {
  
  #################################################################################
  ######################### ##### prepare analysis  ##### #########################
  #################################################################################
  
  #-----------------------------------------------
  ### load TGGATES data, normalize, and curate ###
  #-----------------------------------------------

  ## directories to store TGGATES data
  fd <- file.path(lrdataDir,id)
  
  ## name of eset object
  if (tissue=="")
    objname <- sprintf("EMTAB%s.eSet",id)
  else
    objname <- sprintf("EMTAB%s.eSet_%s",id,tissue)
  
  ## full file names, these are the names from normalize.TGGATES function
  fdn <- file.path(fd, sprintf("%s.rda", objname))
  
  ## download and process TGGATES data, object names are set in normalize.TGGATES function
  if (!file.exists(fdn)) {
    tissue.orig <- tissue
    if (tissue.orig=="")
      tissue <- "liver"
    accession <- paste0("E-MTAB-",id)
    srcdir <- NULL
    if (file.exists(file.path(fd,paste0(accession,"_eSet_orig.rda")))) srcdir <- fd
#     srcdir <- ifelse(file.exists(file.path(fd,paste0(accession,"_eSet_orig.rda"))), fd, NULL)
    EMTAB <- normalize.TGGATES(accession, outdir=fd, sourcedir=srcdir, verbose=T, tissue=tissue)
    tissue <- tissue.orig
    assign(objname, EMTAB)
  } else {
    print(sprintf("Loading %s ...",fdn))
    load(fdn)
    EMTAB <- get(objname)
  }
  
  # test if eSet names were assigned correctly
  if (tissue=="")
    stopifnot(sapply(c(sprintf("EMTAB%s.eSet",id)), exists, envir=environment()))
  else 
    stopifnot(sapply(c(sprintf("EMTAB%s.eSet_%s",id,tissue)), exists, envir=environment()))

  print("Expression set ready")
  
  #-------------------------------------------------
  ### only use drugs with full range of experiments
  #-------------------------------------------------
  
  EMTAB.filtered <- keepDrugs(EMTAB, id=id, tissue=tissue)
  if (tissue=="")
    save(EMTAB.filtered, file=sprintf("%s_EMTAB.filtered",id))
  else 
    save(EMTAB.filtered, file=sprintf("%s_%s_EMTAB.filtered",id,tissue))
  print("TGGATES drugs filtered ready")
  
  ###############################################################################
  ######################### ##### start analysis  ##### #########################
  ###############################################################################
  
  #----------------
  ###  run GLM  ###
  #----------------
  
  if (tissue=="")
    res.glm.fn <- file.path(rdataDir, sprintf("%s_res.glm.rda",id))
  else
    res.glm.fn <- file.path(rdataDir, sprintf("%s.%s_res.glm.rda",id,tissue))
  if (!file.exists(res.glm.fn)) {
    print("Fit linear models")
    res.glm <- runGLM(EMTAB.filtered, id=id, cores=cores, intermediateDir=lintermediateDir, tissue=tissue)
    save(res.glm, file=res.glm.fn)
  } else load(res.glm.fn)
  print("Linear models ready")

  #----------------------------------------------------------------------------------------------------------
  ### Compare distributions of number of significant genes for hepatocarcinogens vs non-hepatocarcinogens ###
  ### This analysis is added on Jan 13, 2015 in response to reviewers of EHP                              ###
  #----------------------------------------------------------------------------------------------------------
  
  if (tissue=="")
    fn <- sprintf("wilcox.%s.pdf",id)
  else
    fn <- sprintf("wilcox.%s.%s.pdf",id,tissue)

  title <- switch(id,
                  "797"=paste("PRH", tissue),
                  "798"=paste("PHH", tissue),
                  "799"=paste("RLV", tissue),
                  "800"=paste("RLV repeated", tissue))
  xx <- compare.GLMgenes(res.glm, fn=file.path("plots",fn), title=title)
  print("ID:")
  print(id)
  print("Number of genes hepatotoxicant:")
  print(length(xx$numberSigGenes.hepato))
  print("Number of genes non-hepatotoxicant:")
  print(length(xx$numberSigGenes.nonHepato))
  print("Number of genes intersecting:")
  print(length(intersect(xx$numberSigGenes.hepato,xx$numberSigGenes.nonHepato)))
  print("Overlap %:")
  print(length(intersect(xx$numberSigGenes.hepato,xx$numberSigGenes.nonHepato))*100/length(union(xx$numberSigGenes.hepato,xx$numberSigGenes.nonHepato)))

  #-----------------------------
  ###  process GLM  results ###
  #-----------------------------
  
  if (tissue=="")
    res.glm.coeff.fn <- file.path(rdataDir, sprintf("%s_res.glm.coeff.rda",id)) 
  else
    res.glm.coeff.fn <- file.path(rdataDir, sprintf("%s.%s_res.glm.coeff.rda",id,tissue)) 
  if (!file.exists(res.glm.coeff.fn)) {
    print("Process linear models")
    res.glm.coeff <- processCoeff(res.glm, cores=cores)
    save(res.glm.coeff, file=res.glm.coeff.fn)
  } else load(res.glm.coeff.fn)
  print("Processing of linear models ready")
  
  #----------------------------
  ###  compute gene ranks ###
  #----------------------------
  
  ## Compute gene ranks ##
  rankFiles <- makeRanks(res.glm.coeff, id=id, rankDir=lrankDir, cores=cores, tissue=tissue)
  print("Ranks ready")
  
  ###############
  ### results ###
  ###############
  
  return(list(rankFiles=rankFiles,
#                gids=gsub(".+\\.(.+)", "\\1", rownames(res.glm.coeff))))
              gids=gsub(".+\\.(.+)", "\\1", rownames(exprs(EMTAB)))))
  
  ###############
  
  ## Generate gmt file ##
  
#   gmtFileName <- makeGMT(id, gseaDir=lgseaDir, gids.rat, gids.human,
  
  #--------------------------------------------
  ###  perform Gene Set Enrichment Analysis ###
  #--------------------------------------------
  
  ###########################################################
  ###  ####
  
  ##################
  ### blablabla ####
  
  
  ### blabla ###
  
  
  ### bla
  
  ##
}