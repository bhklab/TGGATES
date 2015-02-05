###################################################################################
###
### Perform Gene Set Enrichment Analysis for set of rank files in parallel
### for each accession number
###
### patrick@jimmy.harvard.edu, Jan 23, 2014
###
###################################################################################

#' Perform GSEA on set of rank files
#' @param ranks vector of rank file names (full path)
#' @param id TGGATES accession number (used to name GSEA output folder)
#' @param minSize maxSize minimal/maximal gene set size to consider
#' @param gmtPath full path to gmt file
#' @param genesetDescription description of gene set obtained from formatGMT function
#' @param nPerm number of permutations to be performed for permutation test in GSEA
#' @param cores number of cores for parallelization
#' @return list of gathered gsea results
perform.pipelineC_for_id <- function(id, ranks, min.size, max.size, gmt.path, genesetDescription,
                              out.dir=gseaOutDir, exec.path=gseaExecPath, nPerm=1000, cores=ncores,
                              lrdataDir=rdataDir, lintermediateDir=intermediateDir, tissue="") {
  require(MetaGx) # use Ben's preRankedGSEA function in MetaGx package
  
  id <- ifelse(tissue=="",id,paste(id,tissue,sep="."))
  
  # store results in separate folders
  outDir.id <- file.path(out.dir,id)
  dir.create(outDir.id, showWarnings=F)
  
#   # make full path for files
#   ranks <- sapply(ranks, function(x) file.path(getwd(), x))
#   gmt.path <- file.path(getwd(), gmt.path)
  
  # perform GSEA in parallel
  fn.res.gsea <- file.path(lrdataDir,sprintf("%s_%s_res.gsea.rda",id,basename(gmt.path)))
  if (!file.exists(fn.res.gsea)) {
    print("Running GSEA")
#     message("Please note, GSEA on rank files is not in parallel (JAVA VM memory issue)")
    res.gsea <- mclapply(ranks, function(rankFile) {
      fn.tmp.rankFile <- file.path(lintermediateDir,sprintf("%s_%s_gsea_tt.rda",basename(rankFile),basename(gmt.path)))
      if (!file.exists(fn.tmp.rankFile)) {
        ## perform GSEA via gsea.prerank of MetaGx
        res.prerankGSEA <- MetaGx::prerankGSEA(exe.path=exec.path, gmt.path=gmt.path, rank.path=rankFile, gsea.collapse=FALSE, 
                                               nperm=nPerm, scoring.scheme="weighted", make.sets=TRUE, include.only.symbols=FALSE, 
                                               plot.top.x=20, set.max=max.size, set.min=min.size, zip.report=FALSE,
                                               gsea.out=outDir.id, replace.res=FALSE, gsea.seed=987654321)
        ## separate leading edge from table
        tt.leading <- res.prerankGSEA$geneset.core
        tt <- res.prerankGSEA$geneset.results
        ## gather results
        # columsn 1 2 3 are NAME, GS.br..follow.link.to.MSigDB, GS.DETAILS, respectively. 12 doesn't even exist?
        tt <- data.frame("geneset"=rownames(tt), tt[ , -c(1, 2, 3, 12)]) #, "description"=genesetDescription[rownames(tt)])
        ## adjust P-Values !!! note please in methods section or supplements
        tt[!is.na(tt[ , "NOM.p.val"]) & tt[ , "NOM.p.val"] == 0, "NOM.p.val"] <- 1 / (nPerm + 1)
        res.prerankGSEA.orig <- res.prerankGSEA
        GSEA <- list(res.gsea=tt, leadingEdge=tt.leading)
        save(GSEA,  file=fn.tmp.rankFile)
      } else {
        message(sprintf("Skipping %s because saved previously", fn.tmp.rankFile))
        load(fn.tmp.rankFile)
      }
      return(GSEA)
    }, mc.cores=cores)
    names(res.gsea) <- names(ranks) # drug names
    res.gseas <- lapply(res.gsea, "[[", "res.gsea")
    res.leadingEdges <- lapply(res.gsea, "[[", "leadingEdge")
    gsea.finalRes <- list(res.gseas.full=res.gsea, res.gseas=res.gseas, res.leadingEdges=res.leadingEdges)
    save(gsea.finalRes, file=fn.res.gsea)
    print(sprintf("Saved res.gsea to %s",fn.res.gsea))
  } else {
    print(sprintf("Skipping %s because exists. Loading from file instead..", fn.res.gsea))
    load(fn.res.gsea)
  }  
  gsea.finalRes
}