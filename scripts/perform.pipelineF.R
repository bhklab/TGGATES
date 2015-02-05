######################################################
###
### Cluster NES scores
###
### patrick@jimmy.harvard.edu, March 4, 2014
###
######################################################

#' @param modules annotated result of isa as returned from biClustering()
#' function in pipeline D
#' @param title used to name heatmap pdf
#' @param min.val assumed minimal value in heatmap (in NES scores)
plot.NESheatmaps <- function(modules, nesTable, title, plotDir, 
                             lowCol="green", highCol="blue",
                             min.val=-2.5, max.val=2.5) {
  require(gplots)
  require(eisa)
  
  genesets <- getFeatureNames(modules) # because rows genesets
  drugs <- getSampleNames(modules)
  
  bk <- seq(min.val, max.val, by=0.25)
  cols <- colorpanel(length(bk)-1, lowCol, highCol)
  
  # that's why functional code is cool actually, no need for taking care of index pitfalls
  if (length(modules) > 1)
  for (i in 1:length(modules)) {
    dir.create(plotDir, "heatmaps", showWarnings=F)
    pdf(file.path(plotDir, "heatmaps", sprintf("heatmap_%s_module%s.pdf",title,i)), 
        height=18, width=63)
    tt <- nesTable[genesets[[i]], drugs[[i]]]
    dataset <-  gsub("\\.", " ", gsub("modules(\\d\\d\\d\\.*.*)","\\1",title))
    dataset <- gsub("797", "PRH", dataset)
    dataset <- gsub("798", "PHH", dataset)
    dataset <- gsub("799", "RLV", dataset)
    dataset <- gsub("800", "RLV repeated", dataset)
    dataset <- gsub("_", "", dataset)
    dataset <- gsub("Reactome", "", dataset)
    dataset <- gsub("non-redundant", "", dataset)
    main <- paste(dataset, "module ", i)
    if (nrow(tt) >= 2 && ncol(tt) >= 2) {
      heatmap.2(tt, breaks=bk, col=cols, trace="none", srtCol=45, main=main, cex.main=2.7,
                dendrogram="none", margins = c(19, 87), cexCol = 1.8, cexRow= 1.75, key=F,
                labRow=gsub("_"," ",rownames(tt)), labCol=gsub("_"," ",colnames(tt)))
    } else {
      plot(1,1,type="n",main="MODULE LESS THAN TWO COLUMNS OR ROWS")
    }
    dev.off()
  }
  else warning("Length of modules is zero, nothing to plot.")
}

#' Takes commong pathways, applies ordering according to clustering of first element given
#' 
#' @modules modules list(!) of modules, names should be module names
#' @param nesTables list(!) of corresponding nesTables, element at pos i matches modules at pos i
plot.NESheatmaps.common <- function(modules, nesTables, title, lplotDir, lowCol="green", highCol="blue",
                                    min.val=-2.5, max.val=2.5) {
  require(gplots)
  require(eisa)
  
  genesets <- lapply(lapply(modules, eisa::getFeatureNames), "[[", 1)
  
  genesets.inters <- Reduce(intersect, genesets)
  
  bk <- seq(min.val, max.val, by=0.25)
  cols <- colorpanel(length(bk)-1, lowCol, highCol)
  
  if (length(modules) > 1)
    for (i in 1:length(modules)) {
      thename <- ifelse(!is.null(names(modules)), names(modules)[[i]], i)
      pdf(file.path(lplotDir, sprintf("heatmap_common_%s_module%s.pdf",title,thename)), height=12, width=48)
      drugs <- eisa::getSampleNames(modules[[i]])[[1]]
      tt <- nesTables[[i]][genesets.inters, drugs]
      if (nrow(tt) >= 2 && ncol(tt) >= 2) {
        heatmap.2(tt, breaks=bk, col=cols, trace="none", Rowv=F, 
                  dendrogram="none", margins = c(25, 72), cexCol = 1.9, cexRow= 1.7, key=F, srtCol=45,
                  labRow=gsub("_"," ",rownames(tt)), labCol=gsub("_"," ",colnames(tt)))
      } else {
        plot(1,1,type="n",main="MODULE LESS THAN TWO COLUMNS OR ROWS")
      }
      dev.off()
    }
  else warning("Length of modules is zero, nothing to plot.")
}

#' Takes commong drugs and pathways, applies ordering according to clustering of first element given
#' 
#' @modules modules list(!) of modules, names should be module names
#' @param nesTables list(!) of corresponding nesTables, element at pos i matches modules at pos i
plot.NESheatmaps.common_2 <- function(modules, nesTables, title, lplotDir, lowCol="green", highCol="blue",
                                    min.val=-2.5, max.val=2.5) {
  require(gplots)
  require(eisa)
  
  genesets <- lapply(lapply(modules, eisa::getFeatureNames), "[[", 1)
  drugs <- lapply(lapply(modules, eisa::getSampleNames), "[[", 1)
  
  genesets.inters <- Reduce(intersect, genesets)
  drugs.inters <- Reduce(intersect, drugs)
  
  bk <- seq(min.val, max.val, by=0.25)
  cols <- colorpanel(length(bk)-1, lowCol, highCol)
  
  ## ordering according to clustering of first element given
  tmpClust <- heatmap.2(nesTables[[1]][genesets.inters, drugs.inters], breaks=bk, col=cols, trace="none", Rowv=F, dendrogram="none") # don't cluster pathways
  
  genesets.ordered <- genesets.inters
  drugs.ordered <- drugs.inters[rev(tmpClust$colInd)]
  
  if (length(modules) > 1)
    for (i in 1:length(modules)) {
      thename <- ifelse(!is.null(names(modules)), names(modules)[[i]], i)
      pdf(file.path(lplotDir, sprintf("heatmap_common_%s_module%s.pdf",title,thename)), height=12, width=45)
      heatmap.2(nesTables[[i]][genesets.ordered, drugs.ordered], breaks=bk, col=cols, trace="none",
                dendrogram="none", margins = c(25, 72), cexCol = 2.4, cexRow= 1.7, key=F)
      dev.off()
    }
  else warning("Length of modules is zero, nothing to plot.")
}

#### euler diagrams ####

plot.vennEuler <- function(genesetNames.rat, genesetNames.human, genesetSource, plotDir) {
  require(VennDiagram)
  
  ## see example(vennDiagram, with filename = "Venn_4set_pretty.tiff",) ##
  
  thelist <- list(Rat=genesetNames.rat, Human=genesetNames.human)
# thelist=list(A=A,B=B)
  
  pdf(file.path(plotDir,sprintf("euler_%s.pdf",genesetSource)))
  plot.new()
  venn.plot <- venn.diagram(
    x = thelist,
    filename=NULL,
    col = "transparent",
    fill = c("cornflowerblue", "green"),
    alpha = 0.630,
    label.col = c("black", "white", "black"),
#                   , "white", 
#                   "white", "white", "white", "white", "darkblue", "white", 
#                   "white", "white", "white", "darkgreen", "white"),
    cex = 1.7,
    fontfamily = "serif",
    fontface = "bold",
    cat.col = c("darkblue", "darkgreen"),
    cat.cex = 1.5,
    cat.pos = 0,
    cat.dist = 0.09,
    cat.fontfamily = "serif",
#     rotation.degree = 270,
    margin = 0.2,
    main=genesetSource
#     euler.d=TRUE,
#     scaled=T
  );
  grid.draw(venn.plot)
  dev.off()
}

#' Plot heatmaps for all modules
#' Plot number of overlapping genesets (venn/euler diagram)
#' @param modulesForGenesetSource list of biclustering results including also
#' NES tables
#' @param genesetNames list A_x, where elements of A_x are human and rat
perform.pipelineF <- function(modulesForGenesetSource, genesetSource, genesetNames, lplotDir=plotDir) {
  
  ## HEATMAPS ##
  
  ## redundant modules
  themodules <- modulesForGenesetSource$modules
                 
  # dataset will equal module797, etc.
  for (dataset in names(themodules)) {
    modules <- themodules[[dataset]]$res.biclustering$modules
    nesTable <- themodules[[dataset]]$nesTable
    plot.NESheatmaps(modules, nesTable, title=sprintf("%s_%s",genesetSource,dataset), plotDir=lplotDir)
  }
  
  ## non-redundant modules
  themodules.nonredundant <- modulesForGenesetSource$modules.nonredundant
  
  # dataset will equal modules797, etc.
  for (dataset in names(themodules.nonredundant)) {
    modules <- themodules.nonredundant[[dataset]]$modules
    nesTable <- themodules[[dataset]]$nesTable
    plot.NESheatmaps(modules, nesTable, title=sprintf("%s_%s_non-redundant",genesetSource,dataset), plotDir=lplotDir)
  }
  
  ## plot module 2 in the datasets with common pathways
  if (genesetSource=="Reactome") {
    ## first panel figure 3
    mymodules <- list(RLV=themodules.nonredundant$modules799$modules[[2]], 
                      PHH=themodules.nonredundant$modules798$modules[[15]],
                      PRH=themodules.nonredundant$modules797$modules[[10]])
    mynestables <- list(RLV=themodules$modules799$nesTable,
                        PHH=themodules$modules798$nesTable,
                        PRH=themodules$modules797$nesTable)
    plot.NESheatmaps.common(modules=mymodules, nesTables=mynestables,
                            title=sprintf("%s_non-redundant_RLV2vsPHH15vsPRH10",genesetSource), 
                            lplotDir=lplotDir)
#     ## second panel figure 3
#     mymodules <- list(RLV=themodules.nonredundant$modules799$modules[[24]], PHH=themodules.nonredundant$modules798$modules[[4]])
#     mynestables <- list(RLV=themodules$modules799$nesTable, PHH=themodules$modules798$nesTable)
#     plot.NESheatmaps.common(modules=mymodules, nesTables=mynestables, title=sprintf("%s_non-redundant_RLV24vsPHH4",genesetSource), lplotDir=lplotDir)
  }
  
  ## VENN/EULER ##
  
#   for (genesetSource in names(genesetNames)) {
    plot.vennEuler(genesetNames.rat= genesetNames$rat,
                   genesetNames.human= genesetNames$human,
                   genesetSource= genesetSource, plotDir=lplotDir)
#   }
  
}

