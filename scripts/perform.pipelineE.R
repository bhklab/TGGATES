#########################################################################################################
###
### Create some nice figures for paper. Start with barplot representing mutually overlapping modules
###
### patrick@jimmy.harvard.edu, Feb 14, 2014
###
#########################################################################################################


#' For two matrices after testing for significance in overlapping gene set names of modules
#' (result hypergeometric tests), find overlapping rows (assuming rows of both tables
#' correspond to same dataset) for both matrices
#' @return data.frame with rownames and colnames according to rownames(m1)==rownames(m2) and 
#' colnames according to m1.name and m2.name organisms containing TRUE and FALSE for overlap in both
#' directions (organisms)
getOverlappingRows <- function(m1, m2, m1.name, m2.name, sig.thres) {
  # if row names not correspond to same organism
  # note that rownames are assigned in function where hypergeometric tests are performed
  stopifnot(all(rownames(m1) == rownames(m2))) 
  
  df <- data.frame(apply(m1, 1, function(x) any(x < sig.thres, na.rm=T)), 
                   apply(m2, 1, function(x) any(x < sig.thres, na.rm=T)),
                   row.names=rownames(m1))
  colnames(df) <- c(m1.name,m2.name)
  df
}

#' Count how many Overlaps of one dataset (rows) in others (columns)
#' @param data frame (two colums, n rows where n is number of modules) 
#' with TRUE and FALSE. Result of getOverlappingRows()
#' @return vector of counts where number of elements is length(criteria)
getCriteriaCounts <- function(binary.df) {
  stopifnot(ncol(binary.df)==2) # these functions are implemented for comparison of one ds with two others (2 columns)
  
  criteria <- c("Unique", "Conserved", "Dataset1_only", "Dataset2_only")
  
  res <- sapply(criteria, function(x) {
    fToApply <- switch(x,
                       "Unique"=function(y) all(!y),
                       "Conserved"=function(y) all(y),
                       "Dataset1_only"=function(y) y[1] & !y[2],
                       "Dataset2_only"=function(y) !y[1] & y[2])
    apply(binary.df, 1, fToApply)                       
  })
  # res is a data frame with nrow(binary.df) rows and length(criteria) colums. make it
  # a count vector of length four counting all TRUE in columns
  if (!is.null(dim(res))) # if only one module exists at all
    return(apply(res, 2, sum))
  else {
    names(res) <- criteria
    return(res)
  }
    
}

#' This function is only generic for the number of accessions, i.e. number of elements
#' in 'overlaps' parameter (one element for one accession). Criteria on overlaps (
#' i.e. what to count) is within the function (criteria variable)
#' @param overlaps list of overlaps to plot (e.g. 797, 798, 799). should be named!
#' @param adjust.method method used to correct for multiple testings. defaults in no correction
#' @param sig.thres threshold for p/q value (alpha for pValue)
#' @return result data.frame to be plotted!
getNumberOfModules.table <- function(overlaps, adjust.method="none", sig.thres) {
    
  # rat vivo VS rat vitro and human vitro
  df.rvRvt_rvHvt <- getOverlappingRows(m1=overlaps$RLV.PRH, m2=overlaps$RLV.PHH,
                                       m1.name="RLV.PRH", m2.name="RLV.PHH", sig.thres=sig.thres)
  # human vitro VS rat vivo and rat vitro
  df.HvtRv_HvtRvt <- getOverlappingRows(m1=t(overlaps$RLV.PHH), m2=t(overlaps$PRH.PHH),
                                       m1.name="RLV.PHH", m2.name="PRH.PHH", sig.thres=sig.thres)
  # rat vitro VS rat vivo and human vitro
  df.RvtRv_RvtHvt <- getOverlappingRows(m1=t(overlaps$RLV.PRH), m2=overlaps$PRH.PHH,
                                        m1.name="RLV.PRH", m2.name="PRH.PHH", sig.thres=sig.thres)
  
  # named vectors
  df.rvRvt_rvHvt.Counts <- getCriteriaCounts(df.rvRvt_rvHvt) 
  df.HvtRv_HvtRvt.Counts <- getCriteriaCounts(df.HvtRv_HvtRvt)
  df.RvtRv_RvtHvt.Counts <- getCriteriaCounts(df.RvtRv_RvtHvt)
 
#   res.df.rownames <- c("Unique", "Conserved", "RLV_vs._PRH", "PRH_vs._PHH", "PHH_vs._RLV")
#   res.df.colnames <- c("PRH", "PHH", "RLV") # 797, 798, 799
  
  res.RLV <- c(Unique=df.rvRvt_rvHvt.Counts[["Unique"]],
                    "Conserved"=df.rvRvt_rvHvt.Counts[["Conserved"]],
                    "RLV_vs._PRH"=df.rvRvt_rvHvt.Counts[["Dataset1_only"]],
                    "PRH_vs._PHH"=0,
                    "PHH_vs._RLV"=df.rvRvt_rvHvt.Counts[["Dataset2_only"]])
  
  res.PHH <- c(Unique=df.HvtRv_HvtRvt.Counts[["Unique"]],
                       "Conserved"=df.HvtRv_HvtRvt.Counts[["Conserved"]],
                       "RLV_vs._PRH"=0,
                       "PRH_vs._PHH"=df.HvtRv_HvtRvt.Counts[["Dataset2_only"]],
                       "PHH_vs._RLV"=df.HvtRv_HvtRvt.Counts[["Dataset1_only"]])
  
  res.PRH <- c(Unique=df.RvtRv_RvtHvt.Counts[["Unique"]],
                     "Conserved"=df.RvtRv_RvtHvt.Counts[["Conserved"]],
                     "RLV_vs._PRH"=df.RvtRv_RvtHvt.Counts[["Dataset1_only"]],
                     "PRH_vs._PHH"=df.RvtRv_RvtHvt.Counts[["Dataset2_only"]],
                     "PHH_vs._RLV"=0)
  
  stopifnot(all(names(res.RLV)==names(res.PHH)) &
              all(names(res.PHH)==names(res.PRH))) # maybe a typo bug somewhere, but really bug then, not user input
  
  res.df <- data.frame(PRH=res.PRH,
                       PHH=res.PHH,
                       RLV=res.RLV)
  res.df
  
}

#' function to plot barplot by df
barplot1 <- function(df, title, plotPercentage=T, ...) {
  require(ggplot2) # for plotting
  require(reshape2) # for melt
  require(scales) # for percent labels
  colnames(df) <- gsub("_"," ",colnames(df))
  rownames(df) <- gsub("_"," ",rownames(df))
#   if (plotPercentage) # outcommented so I can use scale_y_continuous(labels = percent_format())
#     df <- prop.table(as.matrix(df), 2)
  mat <- melt(as.matrix(df))
  tmp <- mat
  names(tmp) <- c("Modules", "Dataset", "Count")
  
  tmp$Dataset <- factor(tmp$Dataset, levels=c("PHH", "PRH", "RLV"))
  tmp$Modules <- factor(tmp$Modules, levels=c("Conserved", "PHH vs. RLV", "PRH vs. PHH", "RLV vs. PRH", "Unique"))

  if (plotPercentage) {
    require(plyr)
#     tmp <- ddply(tmp, .(Dataset), transform, pos = cumsum(Count) - 0.5*Count)
#     p <- ggplot(tmp, aes(x = Dataset, y = Count)) +
#       geom_bar(aes(fill = Modules)) +
#       geom_text(aes(label = Count, y = pos), size = 3) + scale_y_continuous(labels = percent_format())
    p <- ggplot(tmp, aes(x=Dataset, y=Count, fill=Modules, order=as.numeric(Modules)))
    p <- p + geom_bar(position="fill", stat="identity") +
      scale_y_continuous(labels = percent_format()) #+ stat_bin(geom = "text", aes(label = Count)) 
  }
  else {
    p <- ggplot(tmp, aes(x=Dataset, y=Count, fill=Modules))
    p <- p + geom_bar(stat="identity")
  }
  # http://stackoverflow.com/questions/9563368/create-stacked-percent-barplot-in-r
#     geom_bar(stat="identity") + # outcommented so I can use scale_y_continuous(labels = percent_format()),  better would be prop.table m aybe?

  p <-  p + ylab("Proportion of modules") +
  theme(axis.text.x = element_text(colour="grey20",size=14,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=13,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.title.y = element_text(size=13),
        axis.title.x = element_blank()) +
  theme(legend.text=element_text(size=12)) +
  scale_fill_brewer(type="qual", palette=1) +
  ggtitle(title)

  list(p=p, df=df)
#   ggplot(tmp, aes(x=Type, y=Rate, fill=factor(Year))) +
#     geom_bar(stat="identity", position="dodge", colour="black") +
#     scale_fill_brewer(type="qual", palette=1)
}

plot.theBarplot <- function(listOfOverlaps, title, barplot.fun=barplot1, ...) {
  df <- getNumberOfModules.table(listOfOverlaps, ...)
  print(df)
  barplot.fun(df, title=title, ...)
}


#' @param themodules list of module results for each organism, i.e. output of perform.pipelineD_for_id()
perform.pipelineE <- function(themodules, title, lplotDir=plotDir, ...) {
  
  ################################
  ### PRINT CONSERVED BARPLOTS ###
  
  ### full modules ###
  
  # these are tables of hypergeometric tests
  overlaps.RLV.PRH <- themodules$overlaps$overlaps.RLV.PRH$overlaps
  overlaps.RLV.PHH <- themodules$overlaps$overlaps.RLV.PHH$overlaps
  overlaps.PRH.PHH <- themodules$overlaps$overlaps.PRH.PHH$overlaps
  
  listOfOverlaps <- list(RLV.PRH=overlaps.RLV.PRH,
                         RLV.PHH=overlaps.RLV.PHH,
                         PRH.PHH=overlaps.PRH.PHH)
  
  p.res <- plot.theBarplot(listOfOverlaps, title=title, ...)  
  p <- p.res$p
  ggsave(p, filename=file.path(lplotDir,sprintf("%s.pdf",gsub(" ","_",title))))
  
  ### nonredundant modules ###
  
  # these are tables of hypergeometric tests
  overlaps.RLV.PRH.nonredundant <- themodules$overlaps.nonredundant$overlaps.RLV.PRH$overlaps
  overlaps.RLV.PHH.nonredundant <- themodules$overlaps.nonredundant$overlaps.RLV.PHH$overlaps
  overlaps.PRH.PHH.nonredundant <- themodules$overlaps.nonredundant$overlaps.PRH.PHH$overlaps
  
  listOfOverlaps.nonredundant <- list(RLV.PRH.nonredundant=overlaps.RLV.PRH.nonredundant,
                                      RLV.PHH.nonredundant=overlaps.RLV.PHH.nonredundant,
                                      PRH.PHH.nonredundant=overlaps.PRH.PHH.nonredundant)
  
  title <- paste0(title,".nonredundant")
  title <- gsub("\\.","",title)
  p.res.nonredundant <- plot.theBarplot(listOfOverlaps.nonredundant, title=title, ...)  
  p.nonredundant <- p.res.nonredundant$p
  ggsave(p.nonredundant, filename=file.path(lplotDir,sprintf("%s.pdf",gsub(" ","_",title))))
  
  ### repeated dose nonredundant modules ###
  
  # these are tables of hypergeometric tests
  overlaps.RLVRepeated.PRH.nonredundant <- themodules$overlaps.nonredundant$overlaps.RLVRepeated.PRH$overlaps
  overlaps.RLVRepeated.PHH.nonredundant <- themodules$overlaps.nonredundant$overlaps.RLVRepeated.PHH$overlaps
  overlaps.PRH.PHH.nonredundant <- themodules$overlaps.nonredundant$overlaps.PRH.PHH$overlaps
  
  listOfOverlaps.nonredundant <- list(RLV.PRH.nonredundant=overlaps.RLVRepeated.PRH.nonredundant,
                                      RLV.PHH.nonredundant=overlaps.RLVRepeated.PHH.nonredundant,
                                      PRH.PHH.nonredundant=overlaps.PRH.PHH.nonredundant)
  
  title <- paste0(title,".nonredundant_repeatedDose")
  title <- gsub("\\.","",title)
  p.res.nonredundant_repeatedDose <- plot.theBarplot(listOfOverlaps.nonredundant, title=title, ...)  
  p.nonredundant_repeatedDose <- p.res.nonredundant_repeatedDose$p
  ggsave(p.nonredundant_repeatedDose, filename=file.path(lplotDir,sprintf("%s.pdf",gsub(" ","_",title))))
  
  list(p.res=p.res,
       p.res.nonredundant=p.res.nonredundant,
       p.res.nonredundant_repeatedDose=p.nonredundant_repeatedDose)
}