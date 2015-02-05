################################################################
###
### Creates gmt files from gene ids resulting from Pipeline A
###
### patrick@jimmy.harvard.edu, Jan 22, 2014
###
################################################################

### credits to Ben ###
formatGMT <- function(infile, outfile, replace=FALSE, verbose=TRUE) {
  ## this function read a gmt file and replace all weird characters in gene set names by "_"
  
  if(file.exists(outfile)) {
    if(!replace) { stop("Output file already exists!") }
    file.remove(outfile)
  }
  
  suppressPackageStartupMessages(require(GSA)) || stop("Library GSA is not available!")
  
  badchars <- "[,][:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
  
  if(verbose) { message(sprintf("reading %s", basename(infile))) }
  tempff <- file.path(dirname(outfile), basename(tempfile(pattern="formatGMT_", tmpdir="", fileext=".tmp")))
  sink(tempff)
  rr <- GSA::GSA.read.gmt(filename=infile)
  sink()
  file.remove(tempff)
  dupln <- sum(duplicated(rr$geneset.names))
  rr$geneset.names <- gsub(pattern=badchars, replacement="_", x=toupper(rr$geneset.names))
  dupln <- c(dupln, sum(duplicated(rr$geneset.names)))
  if(dupln[1] > 0) { warning("Duplicated geneset names in the original gmt file") }
  if(any(dupln[-1] > dupln[1])) { warning("duplicated geneset names due to formatting") }
  names(rr$genesets) <- names(rr$geneset.descriptions) <- rr$geneset.names
  golist <- mapply(c, rr$geneset.names, rr$geneset.descriptions, rr$genesets)  
  
  ## write gmt file
  if(verbose) { message(sprintf("writing %s to %s", basename(outfile), dirname(outfile))) }
  rr2 <- lapply(golist, function(x, file) {
    write(paste(c(x[1], x[2], unique(x[3:length(x)])), collapse="\t"), file=file, append=TRUE)
  }, file=path.expand(outfile))
  invisible(list("geneset"=rr$genesets, "description"=rr$geneset.descriptions))
}
#####

#' Create gmt file for human, rat, and for common using biomaRt. Return common gmt file name
#' CAUTION: uses sourcef function defined in Pipeline_TGGATES. Sources
#' ufoo.R in scripts/ folder. TODO: insource createGMT function.
#' @param id TGGATES accession number, used to return correct gmt file (rat or human)
#' @param gids.rat gene ids for rat
#' @param gids.human gene ids for human
#' @param minn minimal gene set size
#' @param maxx max gene set size
#' @return gmt file name (full path)
makeGMT.go <- function(id, gseaDir, gids.rat, gids.human, minn=15, maxx=500) {
  require(biomaRt)
  sourcef("ufoo.R") ## credits to Ben
  
  #' @param ds dataset for biomaRt depending on species
  #' @param gids gids depending on species
  #' @param outName name of gmt out file for species
  getOneGmt <- function(ds, gids, outName) {
    ## create gmt file
    ensemblrat <- useMart('ENSEMBL_MART_ENSEMBL', dataset=ds, host="www.ensembl.org")
#     go <- getBM(attributes=c("entrezgene", "go_id", "name_1006", "definition_1006"),
#                 filters="entrezgene",values=gids, mart=ensemblrat)
    gobp <- createGMT(gid=gids, value="entrezgene", gokeep="biological_process",
                      mart.db=ensemblrat, outfile=outName, verbose=TRUE)
    
    ## filter according to gene set size ()
    gobp[sapply(gobp, function(x) length(x) > minn+2 & length(x) < maxx+2)]
  }
    
  ## get rat and human gmt files ##
  fn.ratCommon <- file.path(gseaDir,"rat_common_bp.gmt")
  fn.humanCommon <- file.path(gseaDir,"human_common_bp.gmt")
  ## outsourced from big if clause below (see make.reactome)
  fn.rat <- file.path(gseaDir,"rat_bp.gmt")
  fn.human <- file.path(gseaDir,"human_bp.gmt")
  
  if (!file.exists(fn.ratCommon) || !file.exists(fn.humanCommon)) {
    # returns list of gene sets (gmt)
    
    b1 <- !file.exists(fn.rat)
    if (b1) {
      print("Compute rat gmt file")
      time.start <- proc.time()[["elapsed"]]  
      gmt.rat <- getOneGmt("rnorvegicus_gene_ensembl", gids.rat, fn.rat)
      print(sprintf("Took %s seconds to create rat gmt file", proc.time()[["elapsed"]]-time.start))
    } else message("Skipping rat gmt, exists already")
    
    b2 <- !file.exists(fn.human)
    if (b2) {
      print("Compute human gmt file")
      time.start <- proc.time()[["elapsed"]] 
      gmt.human <- getOneGmt("hsapiens_gene_ensembl", gids.human, fn.human)
      print(sprintf("Took %s seconds to create human gmt file", proc.time()[["elapsed"]]-time.start))
    } else message("Skipping human gmt, exists already")   
    # get common
    if (b1 && b2) { # if did not exist yet, but should exist after previous if clauses
      commonNames <- intersect(names(gmt.rat), (names(gmt.human)))
      gmt.rat.common <- gmt.rat[commonNames]
      gmt.human.common <- gmt.human[commonNames]
      lapply(gmt.rat.common, write, file=fn.ratCommon, append=TRUE, ncolumns=1000, sep="\t")  
      lapply(gmt.human.common, write, file=fn.humanCommon, append=TRUE, ncolumns=1000, sep="\t")  
    }    
  }  else message("Skipping common go gmt creation, as rat and human gmt were skipped")
  
  ## format gmt files ##
  fn.ratCommon.formatted <- file.path(gseaDir, "rat_common_bp.formatted.gmt")
  fn.humanCommon.formatted <- file.path(gseaDir, "human_common_bp.formatted.gmt")
#   fn.rat.formatted <- file.path(gseaDir, "rat_bp.formatted.gmt")
#   fn.human.formatted <- file.path(gseaDir, "human_bp.formatted.gmt")
#   
#   if (!file.exists(fn.ratCommon.formatted) || !file.exists(fn.humanCommon.formatted)) {
    gmt.rat.common.formatted <-
      formatGMT(infile=fn.ratCommon, outfile=fn.ratCommon.formatted, replace=TRUE)
    gmt.human.common.formatted <-
      formatGMT(infile=fn.humanCommon, outfile=fn.humanCommon.formatted, replace=TRUE)
#   }

#   ## we need this to visualize number of overlapping genesets
#   gmt.rat.formatted <-
#     formatGMT(infile=fn.rat, outfile=fn.rat.formatted, replace=TRUE)
#   gmt.human.formatted <-
#     formatGMT(infile=fn.human, outfile=fn.human.formatted, replace=TRUE)
  
  ## return correct gmt name
  switch(id,
         "797"=list(fn.common=fn.ratCommon, fn.common.formatted=fn.ratCommon.formatted,
                    gmt.common.formatted=gmt.rat.common.formatted, fn.gmt=fn.rat),
         "798"=list(fn.common=fn.humanCommon, fn.common.formatted=fn.humanCommon.formatted,
                    gmt.common.formatted=gmt.human.common.formatted, fn.gmt=fn.human),
        "799"=list(fn.common=fn.ratCommon, fn.common.formatted=fn.ratCommon.formatted,
                    gmt.common.formatted=gmt.rat.common.formatted, fn.gmt=fn.rat),
        "800"=list(fn.common=fn.ratCommon, fn.common.formatted=fn.ratCommon.formatted,
                    gmt.common.formatted=gmt.rat.common.formatted, fn.gmt=fn.rat))
}

## same parameters as makeGMT.go. Most parts look similar, but makeGTM.go code was
## difficult to reuse (don't wonder when re-reading your own code, Patrick ;)
makeGMT.reactome <- function(id, gseaDir, gids.rat, gids.human, minn=15, maxx=500) {
  require(biomaRt)
  sourcef("ufoo.R") ## credits to Ben
  
  #' @param lspecies either TGGATES accession ID
  #' @param outName name (full path) of gmt out file for species
  getOneGmt <- function(lspecies, outName) {
       ## query biomart
    mart <- useMart('REACTOME',dataset="pathway")
    mart.ds <- useDataset("pathway",mart=mart)
    # give correct species
    species <- switch(lspecies,
                      "rat"="Rattus norvegicus",
                      "human"="Homo sapiens")
    
    attrs <- getBM(attributes=c("referencedatabase_ncbi_gene", "_displayname"), 
                   filters="species_selection", values=species, mart=mart.ds)
    # help> attrs[,1] are gene ids, attrs[,2] pathway names
    
    ## restrict genes to those inputted (the ones from TGGATES experiments)
    gids <- switch(lspecies,
                   "rat"=gids.rat,
                   "human"=gids.human)
    attrs <- attrs[attrs[,1]%in%gids,] 
    
    ## create gmt of desired size
    geneSets <- tapply(attrs[,1], attrs[,2], identity)    
    res.gmt <- geneSets[sapply(geneSets, function(x) length(x) > minn+1 & length(x) < maxx+1)]
    names(res.gmt) <- gsub(" ", "_", names(res.gmt)) # replace white space
    
    ## print gmt file
    res.gmt2 <- lapply(names(res.gmt), function(gsname) c(gsname, gsname, res.gmt[[gsname]]))
    names(res.gmt2) <- names(res.gmt)
    lapply(res.gmt2, write, file=outName, append=TRUE, ncolumns=1000, sep="\t")      
    res.gmt2
  }
    
  ## get gmt for rat and human ##
  
  fn.ratCommon <- file.path(gseaDir, "rat_common_reactome.gmt")  
  fn.humanCommon <- file.path(gseaDir, "human_common_reactome.gmt")  
  
  # following filenames are outsourced from big clause below to obtain
  # the full geneset names of each dataset (not only common) later (euler diagrams)
  fn.rat <- file.path(gseaDir,"rat_reactome.gmt")
  fn.human <- file.path(gseaDir,"human_reactome.gmt")
  
  if (!file.exists(fn.ratCommon) || !file.exists(fn.humanCommon)) {
    # rat  # below fn was defined (see above)
    
    b1 <- !file.exists(fn.rat)
    if (b1) {
      print("Compute rat reactome gmt file")
      time.start <- proc.time()[["elapsed"]]  
      gmt.rat <- getOneGmt("rat", fn.rat)
      print(sprintf("Took %s seconds to create rat reactome gmt file", proc.time()[["elapsed"]]-time.start))
    } else message("Skipping rat reactome gmt, exists already")
    # human # below fn was defined (see above)
    
    b2 <- !file.exists(fn.human)
    if (b2) {
      print("Compute human reactome gmt file")
      time.start <- proc.time()[["elapsed"]]  
      gmt.human <- getOneGmt("human", fn.human)
      print(sprintf("Took %s seconds to create human reactome gmt file", proc.time()[["elapsed"]]-time.start))
    } else message("Skipping human reactome gmt, exists already")
    if (b1 && b2) { # if did not exist yet, but should exist after previous if clauses
      commonNames <- intersect(names(gmt.rat), (names(gmt.human)))
      gmt.rat.common <- gmt.rat[commonNames]
      gmt.human.common <- gmt.human[commonNames]
      lapply(gmt.rat.common, write, file=fn.ratCommon, append=TRUE, ncolumns=1000, sep="\t")  
      lapply(gmt.human.common, write, file=fn.humanCommon, append=TRUE, ncolumns=1000, sep="\t")  
    }   
  } else message("Skipping common reactome gmt creation, as rat and human gmt were skipped")
  
  ## format gmt files ##
  fn.ratCommon.formatted <- file.path(gseaDir, "rat_common_reactome.formatted.gmt")
  fn.humanCommon.formatted <- file.path(gseaDir, "human_common_reactome.formatted.gmt")
#   fn.rat.formatted <- file.path(gseaDir, "rat_reactome.formatted.gmt")
#   fn.human.formatted <- file.path(gseaDir, "human_reactome.formatted.gmt")
  
  #   if (!file.exists(fn.ratCommon.formatted) || !file.exists(fn.humanCommon.formatted)) {
  gmt.rat.common.formatted <-
    formatGMT(infile=fn.ratCommon, outfile=fn.ratCommon.formatted, replace=TRUE)
  gmt.human.common.formatted <-
    formatGMT(infile=fn.humanCommon, outfile=fn.humanCommon.formatted, replace=TRUE)
  #   }
  
  # we need this to visualize number of overlapping genesets
#   gmt.rat.formatted <-
#     formatGMT(infile=fn.rat, outfile=fn.rat.formatted, replace=TRUE)
#   gmt.human.formatted <-
#     formatGMT(infile=fn.human, outfile=fn.human.formatted, replace=TRUE)
  
  ## return correct file names depending on species
  switch(id,
         "797"=list(fn.common=fn.ratCommon, fn.common.formatted=fn.ratCommon.formatted,
                    gmt.common.formatted=gmt.rat.common.formatted, fn.gmt=fn.rat),
         "798"=list(fn.common=fn.humanCommon, fn.common.formatted=fn.humanCommon.formatted,
                    gmt.common.formatted=gmt.human.common.formatted, fn.gmt=fn.human),
         "799"=list(fn.common=fn.ratCommon, fn.common.formatted=fn.ratCommon.formatted,
                    gmt.common.formatted=gmt.rat.common.formatted, fn.gmt=fn.rat),
         "800"=list(fn.common=fn.ratCommon, fn.common.formatted=fn.ratCommon.formatted,
                    gmt.common.formatted=gmt.rat.common.formatted, fn.gmt=fn.rat))
}

##################
### PIPELINE B ###
perform.pipelineB <- function(gids.rat, gids.human, lgseaDir=gseaDir, genesetSize.min=gseaMinSize, genesetSize.max=gseaMaxSize) {
#   gmt.rat <- makeGMT.go("797", gseaDir=lgseaDir, gids.rat=gids.rat, gids.human=gids.human)
  print("Compute GO biological process gmt files")
  gmt.go.human <- makeGMT.go("798", gseaDir=lgseaDir, gids.rat=gids.rat, gids.human=gids.human, minn=genesetSize.min, maxx=genesetSize.max)
  gmt.go.rat <- makeGMT.go("799", gseaDir=lgseaDir, gids.rat=gids.rat, gids.human=gids.human, minn=genesetSize.min, maxx=genesetSize.max)
  print("Compute reactome gmt files")
  gmt.reactome.human <- makeGMT.reactome("798", gseaDir=lgseaDir, gids.rat=gids.rat, gids.human=gids.human, minn=genesetSize.min, maxx=genesetSize.max)
  gmt.reactome.rat <- makeGMT.reactome("799", gseaDir=lgseaDir, gids.rat=gids.rat, gids.human=gids.human, minn=genesetSize.min, maxx=genesetSize.max)
  print("Gmt files ready")
  list(go=list(gmt.rat=gmt.go.rat, gmt.human=gmt.go.human),
       reactome=list(gmt.rat=gmt.reactome.rat, gmt.human=gmt.reactome.human))
}