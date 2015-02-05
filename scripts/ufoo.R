#################################################
## functions
#################################################

## gene set enrichment analysis
# gsea.prerank <- function(exe.path, gmt.path, rank.path, chip.path, gsea.collapse=FALSE, nperm=1000, scoring.scheme=c("weighted", "weighted_p2", "weighted_p1.5", "classic"), make.sets=TRUE, include.only.symbols=FALSE, plot.top.x=20, set.max=500, set.min=15, zip.report=FALSE, gsea.report, gsea.out, replace.res=FALSE, gsea.seed=987654321) {
# 	exe.path <- path.expand(exe.path)
# 	gmt.path <- path.expand(gmt.path)
# 	rank.path <- path.expand(rank.path)
# 	gsea.seed <- as.integer(gsea.seed)
# 	nperm <- as.integer(nperm)
# 	plot.top.x <- as.integer(plot.top.x)
# 	set.max <- as.integer(set.max)
# 	set.min <- as.integer(set.min)
# 	if(missing(gsea.out)) { gsea.out <- "." }
# 	if(missing(chip.path)) { chip.path <- "" } else { chip.path <- path.expand(chip.path) }
# 	if(!gsea.collapse) { gsea.collapse <- "false" } else { gsea.collapse <- "true" }
# 	if(!make.sets) { make.sets <- "false" } else { make.sets <- "true" }
# 	if(!include.only.symbols) { include.only.symbols <- "false" } else { include.only.symbols <- "true" }
# 	if(!zip.report) { zip.report <- "false" } else { zip.report <- "true" }
# 	if(missing(gsea.report)) { gsea.report <- paste("gsea_report", gsub("[.]", "_", gsub("[.]rnk", "", basename(rank.path))), sep="_") }
# 	scoring.scheme <- match.arg(scoring.scheme)
# 	rest <- dir(gsea.out)
# 	rest <- rest[grep(pattern=sprintf("^%s", gsea.report), x=rest)[1]]
# 	if(!replace.res && (length(rest) > 0 && !is.na(rest))) { warning("output directory already exists!") } else {
# 		gsea.cmd <- sprintf("java -Xmx10240m -cp %s xtools.gsea.GseaPreranked -gmx %s -chip %s -collapse %s -nperm %i -rnk %s -scoring_scheme %s -rpt_label %s -include_only_symbols %s -make_sets %s -plot_top_x %i -rnd_seed %i -set_max %i -set_min %i -zip_report %s -out %s -gui false", exe.path, gmt.path, chip.path, gsea.collapse, nperm, rank.path, scoring.scheme, gsea.report, include.only.symbols, make.sets, plot.top.x, gsea.seed, set.max, set.min, zip.report, gsea.out)
# 		system(gsea.cmd)
# 		## read results
# 		rest <- dir(gsea.out)
# 		rest <- rest[grep(pattern=gsea.report, x=rest)[1]]
# 		restn <- sapply(strsplit(rest, "[.]"), function(x) { return(x[length(x)]) })
# 		tt <- rbind(read.csv(file.path(gsea.out, rest, sprintf("gsea_report_for_na_pos_%s.xls",restn)), stringsAsFactors=FALSE, sep="\t", header=TRUE), read.csv(file.path(gsea.out, rest, sprintf("gsea_report_for_na_neg_%s.xls",restn)), stringsAsFactors=FALSE, sep="\t", header=TRUE))
# 		rownames(tt) <- as.character(tt[ ,"NAME"])
# 		## rename results directory
# 		file.rename(from=file.path(gsea.out, rest), to=file.path(gsea.out, gsub(sprintf("[.]GseaPreranked[.]%s", restn), "", rest)))
# 	}
# 	return(tt)
# }

createGMT <- function(gid, value, gokeep=c("biological_process", "molecular_function", "cellular_component"), mart.db, outdir, outfile, verbose=TRUE, ...) {
	#Define the GO terms and create a GMT file for GSEA

	#####################
	## Parameters
	######################
	## gid is a vector of gene ids
	## value is a string specifying what are the gene ids (for example "ensembl_gene_id" or "ensembl_transcript_id")
	## gokeep is the catergory to keep.  Only one GO category
	## mart.db is mart database such as mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
	## outdir is the path to write the output file
	## verbose=TRUE will display some messages
	## ... are additional parameters to be passed to createGMT
	#############################

	## create GO terms for ensembl ids
	require(biomaRt)

	gokeep <- match.arg(gokeep)
	if(missing(gid)) { gid <- getBM(attributes=value, filters="", values="", mart.db, ...)[ ,1] }
	if(missing(outdir)) { outdir <- getwd() }
	dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
	if(missing(outfile)) { outfile <- "GO_TERM.gmt" }
	
	#for all genes, extract GO terms
	gene.an <- getBM(attributes=c(value, "go_id", "name_1006", "definition_1006", "go_linkage_type", "namespace_1003"), filters=value, values=gid, mart=mart.db, ...)
	gene.an[gene.an == "" | gene.an == " "] <- NA
	gene.an <- gene.an[!is.na(gene.an[ ,"namespace_1003"]) & is.element(gene.an[ ,"namespace_1003"], gokeep), ,drop=FALSE]
	gene.an <- data.frame(gene.an, "GONAME"=gsub(pattern="[ ]|[\\]|[/]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]", replacement="_", x=toupper(gene.an[ , "name_1006"])))

	goo <- sort(unique(gene.an[ ,"go_id"]))
	names(goo) <- as.character(gene.an[match(goo, gene.an[ ,"go_id"]), "GONAME"])
	goo2 <- cbind(names(goo), goo)
	rownames(goo2) <- names(goo)
	golist <- apply(goo2, 1, function(x, z) {
	res <- c(x[1], x[2], unique(z[is.element(z[ ,"go_id"], x[2]), value]))
	names(res) <- NULL
	return(res)
	}, z=gene.an)
	names(golist) <- rownames(goo2)

	## write gmt file
	if(verbose) { message(sprintf("writing %s to %s", outfile, outdir)) }
	rr <- lapply(golist, function(x, file) { write(sprintf("%s\thttp://www.ebi.ac.uk/QuickGO/GTerm?id=%s\t%s", x[1], x[2], paste(unique(x[3:length(x)]), collapse="\t")), file=file, append=TRUE) }, file=file.path(outdir, outfile))
	invisible(golist)
}

setcolclass.df <- function(df, colclass, factor.levels) {
	ww <- options()$warn
	options(warn=-1)
	toCls <- function(x, cls) { do.call(paste("as", cls, sep = "."), list(x)) }
	df <- replace(df, , Map(toCls, x=df, cls=colclass))
	options(warn=ww)
	iix <- FALSE
	if(!missing(factor.levels)) { iix <- colclass == "factor" & !is.null(factor.levels) }
	if(any(iix)) {
		for(i in which(iix)) { levels(df[[i]]) <- factor.levels[[i]] }
	}
	return(df)
}
## if ss is the original data.frame and we try to copy it in ss2
## ss2 <- setcolclass.df(df=ss, colclass=sapply(ss, class), factor.levels=sapply(ss, levels)))

## courtesy of Matthew McCall
celfileDateHour <- function(filename) {
	require(affyio)
	h <- affyio::read.celfile.header(filename, info="full")
	#ddate <- grep("/", strsplit(h$DatHeader, " ")[[1]], value=TRUE)
	#ddate <- strsplit(ddate, split="/")[[1]]
	#CC <- ifelse(substr(ddate[3],1,1)=="9", "19", "20")
	if(length(h$ScanDate) > 0) {
	    h$ScanDate <- gsub(pattern="T", replacement=" ", x=h$ScanDate)
	    ddate <- strsplit(h$ScanDate, " ")[[1]]
    } else { ddate <- rep(NA, 2)}
    names(ddate) <- c("day", "hour")
	return(ddate)
}

celfileChip <- function(filename) {
	require(affyio)
	h <- affyio::read.celfile.header(filename, info="full")
	return(as.character(h$cdfName))
}


spearmanCI <- function(x, n, alpha=0.05) {
    require(survcomp)
    zz <- sqrt((n-3)/1.06) * survcomp::fisherz(x)
    zz.se <- 1/sqrt(n - 3)
    ll <- zz - qnorm(p=alpha/2, lower.tail=FALSE) * zz.se
    ll <- survcomp::fisherz(ll / sqrt((n-3)/1.06), inv=TRUE)
    uu <- zz + qnorm(p=alpha/2, lower.tail=FALSE) * zz.se
    uu <- survcomp::fisherz(uu/ sqrt((n-3)/1.06), inv=TRUE)
    pp <- pnorm(q=zz, lower.tail=x<0)
    res <- c("lower"=ll, "upper"=uu, "p.value"=pp)
    return(res)
}

## intersection of more than 2 sets
fold <- function(f, x, y, ...){
    if (missing(...)) { f(x, y) } else { f(x, fold(f, y, ...)) }
}