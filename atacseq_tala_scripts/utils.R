# Author: Florian Halbritter
# https://github.com/cancerbits/lib
# source(paste0(Sys.getenv("CODEBASE"),"/lib/lib/utils.R"))

# apply the same function to each of a given list of named arguments and
# bundle the results in a data.table with the argument names given in an ID column:
rblapply <- function(args, fun, id="id", ..., cores=1) {
  require(data.table)
  if(cores>1) {
    require(parallel)
    res <- parallel::mclapply(X=args, FUN=fun, ..., mc.cores=cores)
    names(res) <- names(args)
    res <- rbindlist(res, idcol=id, fill=T)
    res[,paste(id):=args[get(id)]] # args have been converted to numbers --> convert them back			
  } else {
    res <- rbindlist(sapply(X=args, FUN=fun, simplify=FALSE, ...), idcol=id, fill=T)
  }
  return(res)
}
# same as rblapply but with two sets of arguments:
rblapply2 <- function(args1, args2, fun, id1="id1", id2="id2", cores=1) {
  require(data.table)
  rblapply(args1, function(x) {
    rblapply(args2, function(y) {
      fun(x,y)
    }, id2)
  }, id1, cores=cores)
}

# dt <- rblapply(1:5, function(x) {	data.table(y=x*2) }, "x")
# print(dt)
# # x  y
# # 1: 1  2
# # 2: 2  4
# # 3: 3  6
# # 4: 4  8
# # 5: 5 10
# dt <- rblapply(1:5, function(x) {	data.table(y=x*2) }, "x", cores=4)
# print(dt)
# # x  y
# # 1: 1  2
# # 2: 2  4
# # 3: 3  6
# # 4: 4  8
# # 5: 5 10
# dt <- rblapply(1:5, function(x, another_argument) {	data.table(y=x*another_argument) }, "x", another_argument=5)
# print(dt)
# # x  y
# # 1: 1  5
# # 2: 2 10
# # 3: 3 15
# # 4: 4 20
# # 5: 5 25
# dt <- rblapply2(c("gene1","gene2","gene3"), 1:3, function(g,x) { data.table(y=x*2) }, "gene", "x", cores=4)
# print(dt)
# # gene x y
# # 1: gene1 1 2
# # 2: gene1 2 4
# # 3: gene1 3 6
# # 4: gene2 1 2
# # 5: gene2 2 4
# # 6: gene2 3 6
# # 7: gene3 1 2
# # 8: gene3 2 4
# # 9: gene3 3 6

# # read in a bunch of CSVs in one go:
# f <- list.files("/data_synology_rg3/cancerbits/resources/regions/chipseq_selected/mm10/takahashi_2019/regions", full.names=T, pattern=".bed")
# dt <- rblapply(f, fread, "input_file")


rblapplyDT <- function(dt, fun, idcol) {
  res <- apply(dt, 1, function(x) {
    res <- data.table(fun(as.list(x)))
    res[, paste(idcol):=x[idcol]]
    res
  })
  if(is.list(res)) res <- rbindlist(res)
  return(as.data.table(res))
}

lib <- list(
  projectName = "analysis",
  currentAnalysis = "",
  
  setCurrentAnalysis = function(curAnalysis) {
    lib$currentAnalysis <<- curAnalysis
    if(!dir.exists(lib$plotDir())) dir.create(lib$plotDir(), recursive=TRUE)
  },
  
  getCurrentAnalysis = function() {
    return(lib$currentAnalysis)
  },
  
  privateDir = function(...) {
    if(lib$projectName!="analysis") {
      baseDir <- paste0(Sys.getenv("CODEBASE"), "/", lib$projectName)
      return(paste0(baseDir, "/", ...))		
    }
    else paste0(paste0(getwd(), "/"),...)
  },
  
  commonDir = function(...) {
    if(lib$projectName!="analysis") {
      baseDir <- paste0(Sys.getenv("OUT"), "/", lib$projectName)
      return(paste0(baseDir, "/", ...))		
    }
    else paste0(paste0(getwd(), "/"),...)
  },
  
  analysisDir = function(...) {
    lib$privateDir("src/",...)
  },
  
  dataDir = function(...) {
    lib$commonDir("data/",...)
  },
  
  metaDir = function(...) {
    lib$privateDir("metadata/",...)
  },
  
  resultsDir = function(...) {
    lib$baseResultsDir(lib$currentAnalysis, "/", ...)
  },
  
  baseResultsDir = function(...) {
    lib$commonDir("results/", ...)
  },	
  
  pipelineDir = function(...) {
    lib$resultsDir(...)
  },
  
  plotDir = function(...) {
    lib$resultsDir(...)
  },
  
  resourceDir = function(...) {
    paste0(Sys.getenv("RESOURCES"), "/",...)
  },
  
  toolsDir = function(...) {
    paste0(Sys.getenv("TOOLS"), "/",...)
  },
  
  scriptPath = function(analysisName="main", moduleName=NULL) {
    if(!is.null(moduleName)) {	
      f <- lib$analysisDir(analysisName, "/", moduleName, ".R")
    } else {	
      f <- lib$analysisDir(analysisName, ".R")
    }	
    if(!file.exists(f)) f <- lib$analysisDir("modules/", analysisName, ifelse(is.null(moduleName), "", paste0("_",moduleName)), ".R")			
    return(f)
  },
  
  run = function(analysisName="main", moduleNames=NULL) {
    if(length(moduleNames)>1) {
      lib$forEach(moduleNames, function(n) lib$run(analysisName, n))
    } else {	
      mn <- NULL 
      if(!is.null(moduleNames)) mn <- moduleNames[1]
      f <- lib$scriptPath(analysisName, mn)
      if(!file.exists(f)) stop(sprintf("Cannot find analysis script for %s/%s.", analysisName, moduleNames))
      lib$statusMessageFormatted("+++++ starting %s %s +++++", analysisName, ifelse(is.null(moduleNames), "", moduleNames[1]))
      source(f)
      lib$statusMessageFormatted("----- completed %s %s -----", analysisName, ifelse(is.null(moduleNames), "", moduleNames[1]))
      analysisX <- paste0(analysisName, ":::", mn)
      lib$.completedAnalyses <<- c(lib$.completedAnalyses, analysisX)
    }
  },
  
  .completedAnalyses = c(),
  
  dependOn = function(analysisName="main", moduleNames=NULL, verbose=TRUE) {
    if(length(moduleNames)>1) {
      lib$forEach(moduleNames, function(n) lib$dependOn(analysisName, n))
    } else {	
      mn <- NULL 
      if(!is.null(moduleNames)) mn <- moduleNames[1]
      if(length(intersect(lib$.completedAnalyses, paste0(analysisName, ":::", mn)))==0) {
        lib$run(analysisName, mn)
      } else if(verbose) {
        lib$statusMessageFormatted("++ dependency already satisfied: %s %s", analysisName, mn)
      }		
    }
  },
  
  # runBox = function(analysisName="main", fun) {
  # e <- environment()
  # sp <- lib$scriptPath(analysisName, c("init"))
  # #cmd <- sprintf("source(\"%s\")", paste0(Sys.getenv("CODEBASE"),"/lib/lib/utils.R"))
  # #msg(cmd)
  # #eval(parse(text=cmd), envir=e)
  # cmd <- sprintf("{ source(\"%s\"); %s }", sp, fun)
  # res <- eval(parse(text=cmd), envir=e)
  # return(res)
  # },
  # #runBox("single_cell", "readCombos()")
  
  statusMessage = function(...) {
    message(format(Sys.time(), "%Y-%m-%d %H:%M"), "\t", ...)
  },
  
  statusMessageFormatted = function(...) {
    message(format(Sys.time(), "%Y-%m-%d %H:%M"), "\t", sprintf(...))
  },
  
  # extracted from RnBeads:
  rnb.beta2mval= function(betas, epsilon = 0.00001) {
    betas[betas < epsilon] <- epsilon
    betas[betas > (1 - epsilon)] <- 1 - epsilon
    return(log2(betas / (1 - betas)))
  },
  
  rblapply = rblapply,
  
  splitApply = function(dt, f, fun) {
    sapply(split(dt, by=f), fun)
  },
  
  rowApply = function(mat, fun, ...) {
    x <- apply(mat, 1, fun, ...)
    if(is.null(dim(x))) {
      names(x) <- rownames(mat)
    } else {
      x <- t(x)
      dimnames(x) <- dimnames(mat)
    }		
    x
  },
  
  trimColNames = function(dt,trimPattern,repWith="") {
    colnames(dt) <- gsub(trimPattern,repWith,colnames(dt))
    dt
  },
  
  bitmapPlot = function(name,width,height,pointsize=12,res=150,type="png16m",suffix="png",sub=NULL) {
    bitmap(lib$plotDir(sub,"/",lib$escapeStringInFileName(name),".",suffix), type=type, height=height, width=width, res=res, units="in", pointsize=pointsize)
  },
  
  tiffPlot = function(name,width,height,pointsize=12,res=150,sub=NULL) {
    lib$bitmapPlot(lib$plotDir(lib$escapeStringInFileName(name)),width=width,height=height,pointsize=pointsize,res=res,type="tifflzw",suffix=".tiff")
  },
  
  svgPlot = function(name,width,height,pointsize=12,res=NULL,sub=NULL) {
    svg(lib$plotDir(lib$escapeStringInFileName(name),".svg"), height=height, width=width, pointsize=pointsize)
  },
  
  pdfPlot = function(name,width,height,pointsize=12,res=NULL,useDingbats=FALSE,sub=NULL) {
    pdf(file=lib$plotDir(sub,"/",lib$escapeStringInFileName(name),".pdf"), height=height, width=width, pointsize=pointsize,useDingbats=useDingbats)
  },
  
  pngPlot = function(name,width,height,pointsize=12,res=150,sub=NULL) {
    if(require("Cairo")) {
      CairoPNG(lib$plotDir(sub,"/",lib$escapeStringInFileName(name),".png"), height=height, width=width, pointsize=pointsize, bg="white", res=res, units="in")		
    }
    else {
      png(lib$plotDir(sub,"/",lib$escapeStringInFileName(name),".png"), height=height, width=width, pointsize=pointsize, bg="white", res=res, units="in")
    }
  },
  
  ggsave = function(p, name, width=NA, height=NA, scale=1, res=150, type="png", noWarn=TRUE, expand0=FALSE, addBase=FALSE, sub=NULL) {	
    fname <- lib$plotDir(sub, "/", name, ".", gsub("cairo-","",type))		
    
    if(addBase) p <- lib$ggAddPanelBase(p)
    if(expand0) p <- p + scale_y_continuous(expand=c(0,0))
    #if((type=="png"||type=="cairo-png") && require("Cairo")) type <-  function(...) Cairo(file=fname, width=width, height=height, dpi=res, pointsize=10, type="png", units="in")
    if(noWarn) {
      suppressWarnings(
        ggsave(filename=fname, device = type, plot=p, dpi=res, scale=scale, width=width, height=height, units="in")
      )
    }
    else {
      ggsave(filename=fname, device = type, plot=p, dpi=res, scale=scale, width=width, height=height, units="in")
    }
  },
  
  init = function() {
    lib$loadLibraries(c("RColorBrewer","ggplot2","data.table","fastcluster"))
  },
  
  loadLibraries = function(packages) {
    allOK <- TRUE
    for(package in packages) {
      allOK <- allOK & lib$loadLibrary(package)
    }
    return(allOK)
  },
  
  loadLibrary = function(package, defaultMirror = "http://cran.at.r-project.org") {
    # largely taken from http://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
    
    # first make sure a CRAN mirror has been chosen:
    if(getOption("repos")["CRAN"] == "@CRAN@") {
      options(repos = c(CRAN = defaultMirror))
    }		
    
    if(!eval(parse(text=paste("require(",package,",quietly=TRUE)")))) {
      lib$statusMessage("Installing missing package:", package)
      
      if (!requireNamespace("BiocManager")) install.packages("BiocManager")
      
      # package not installed, try to install from bioconductor:
      tryCatch(BiocManager::install(package),error=function(e) {}) # just ignore error and proceed if this fails
      if(!eval(parse(text=paste("require(",package,",quietly=TRUE)")))) {
        # not a bioconductor package, try to install from CRAN
        
        # then try to install:
        install.packages(package)
        
        return(eval(parse(text=paste("require(",package,",quietly=TRUE)"))))
      }
    }
    
    return(TRUE)
  },
  
  curTimeStr = function() {
    format(Sys.time(), "%Y%m%d_%H%M")
  },
  
  # Transform a label to camelCase by replacing underscores and dots
  # and capitalizing the beginnings of words:
  toCamelCase = function(x){
    x <- tolower(x)
    # inspired by http://stackoverflow.com/questions/11672050/how-to-convert-not-camel-case-to-camelcase-in-r
    capit = function(x) paste0(toupper(substring(x, 1, 1)), substring(x, 2, nchar(x)))
    x <- sapply(strsplit(x, "\\."), function(x) paste(capit(x), collapse=""))
    x <- sapply(strsplit(x, "\\_"), function(x) paste(capit(x), collapse=""))
    x <- paste0(tolower(substring(x, 1, 1)), substring(x, 2, nchar(x)))
    return(x)
  },
  
  padjustMatrix = function(ps, meth="fdr") {		
    matrix(p.adjust(as.numeric(ps),method=meth),ncol=ncol(ps),dimnames=dimnames(ps))
  },
  
  stringReplacements = function(subs, txt) {
    for(n in names(subs)) {
      txt <- gsub(n, subs[n], txt)
    }
    txt
  },
  
  toUnderscoreCase = function(x) {
    x <- gsub("(.)([A-Z][a-z]+)", "\\1_\\2", x)
    x <- tolower(gsub("([a-z0-9])([A-Z])", "\\1_\\2", x))
    return(x)
  },
  
  regsToDt = function(regs) {
    data.table(chrom=as.character(seqnames(regs)), start=start(regs), end=end(regs))
  },
  
  dtToGr = function(dt, chrCol="chrom", startCol="start", endCol="end", metaCols=c()) {
    loadLibrary("GenomicRanges")
    
    argList <- list()
    for(n in metaCols) {
      argList[[n]] <- dt[,get(n)]
    }
    
    argList$ranges <- IRanges(dt[,get(startCol)],dt[,get(endCol)])
    argList$seqnames <- dt[,get(chrCol)]
    
    do.call(GRanges, args=argList)
    
    #do.call(with(dt, GRanges(seqnames=get(chrCol), ranges=IRanges(get(startCol),get(endCol)), arg)))
  },
  
  dfToGr = function(df, chrCol="chrom", startCol="start", endCol="end") {
    loadLibrary("GenomicRanges")
    with(df, GRanges(get(chrCol), IRanges(get(startCol),get(endCol))))
  },
  
  grToDt = function(gr, chrCol="chrom", startCol="start", endCol="end") {		
    dt <- data.table(chrom=as.character(seqnames(gr)), start=start(gr), end=end(gr))
    setnames(dt, c(chrCol, startCol, endCol))
    dt
  },
  
  dt2namedVec = function(d, nameCol="ID", valCol=setdiff(colnames(d),nameCol)) {
    v <- d[,get(valCol)]
    names(v) <- d[,get(nameCol)]
    v	
  },
  
  df2namedVec = function(d, nameCol="ID", valCol=setdiff(colnames(d),nameCol)) {
    v <- d[,valCol]
    names(v) <- d[,nameCol]
    v	
  },
  
  splitBigRegions = function(regs, maxSize=1000) {
    regLen <- regs[,end-start]
    i <- regLen>maxSize
    if(sum(i)>0) {
      msg("split: ", sum(i))
      # split, then call recursively
      regsA <- regs[which(!i),]
      regsB <- regs[which(i),]
      regsB <- rbind(
        regsB[,.(chrom, start, end=round((start+end)/2) )], # first half of regions
        regsB[,.(chrom, start=round((start+end)/2)+1, end=end )] # ... and the second
      )
      rbind(regsA, lib$splitBigRegions(regsB, maxSize=maxSize))
    }
    else {
      # no more regions to split, exit recursion
      return(regs)
    }
  },
  
  mergeAll = function(x) {
    if(!is.list(x)) lib$statusMessage("use with lists only!")
    else {
      commonIds <- rownames(x[[1]])
      for(i in 2:length(x)) {
        commonIds <- intersect(commonIds,rownames(x[[i]]))
      }
      
      mergedD <- x[[1]][commonIds,]
      for(i in 2:length(x)) {
        mergedD <- data.frame(mergedD,x[[i]][commonIds,])
      }
      
      mergedD
    }
  },
  
  escapeStringInFileName = function(str) {
    gsub("[\\+\\(\\)]",".",gsub("/","-", str),perl=TRUE)
  },
  
  corPanelFun = function(x, y, digits=2, prefix='', cex.cor) { 
    usr <- par('usr'); 
    on.exit(par(usr)); 
    par(usr = c(0, 1, 0, 1)); 
    r = (cor(x, y,use='pairwise')); 
    txt <- format(c(round(r,digits), 0.123456789), digits=digits)[1]; 
    txt <- paste(prefix, txt, sep=''); 
    if(missing(cex.cor)) cex <- 0.6/strwidth(txt); 
    text(0.5, 0.5, txt, cex = cex*.6 );
  },
  
  corPanelFunSpearman = function(x, y, digits=2, prefix='', cex.cor) { 
    usr <- par('usr'); on.exit(par(usr)); par(usr = c(0, 1, 0, 1)); r = (cor(x, y,use='pairwise',method='spearman')); txt <- format(c(round(r,digits), 0.123456789), digits=digits)[1]; txt <- paste(prefix, txt, sep=''); if(missing(cex.cor)) cex <- 0.6/strwidth(txt); text(0.5, 0.5, txt, cex = cex*.6 );
  },
  
  getCategoryColors = function(categories, pal=NULL) {
    lib$loadLibrary("RColorBrewer")
    n <- length(categories)
    cols <- NA
    if(!is.null(pal)) {
      cols <- brewer.pal(max(n,3),pal)
    }
    else if(n == 1) cols <- "black"
    else if(n <= 8) {
      cols <- brewer.pal(max(n,3),"Set2")
      if(n < 3) cols <- cols[1:n]
    }
    else if(n <= 9) {
      cols <- brewer.pal(n,"Set1")
    }
    else if(n <= 12) {
      cols <- brewer.pal(n,"Set3")
    }
    else cols <- rainbow(n)
    return(structure(cols,names=as.character(categories)))
  },
  
  determineGroupColours = function(d.annots,col.by=NULL) {	
    lib$loadLibrary("RColorBrewer")
    
    if(!is.null(col.by)) {
      d.annots <- data.frame("condition"=d.annots[,col.by],"group"=d.annots[,col.by])		
    }
    
    col.palette <- c("Greys","Reds","Blues","Greens","Oranges","Purples","YlOrBr","RdPu","PuBuGn")
    for(grpType in intersect(c("condition","condition.short","group"),colnames(d.annots))) {
      t <- table(unique(d.annots[,c("group",grpType)])[,grpType]) #
      if(nrow(t)>0) break
    }	
    if(nrow(t)==0) {
      g <- unique(d.annots$group)
      t <- rainbow(length(g))
      names(t) <- g
      return(t)
    }
    group.cols <- c()
    i <- 1
    for(grp in names(t)) {
      tmp <- rainbow(t[[grp]])
      if(t[[grp]]==1) {
        tmp <- brewer.pal(3,col.palette[i])[3]
      }
      else if(t[[grp]]<8) {
        tmp <- brewer.pal(t[[grp]]+1,col.palette[i])[-1]
      }
      
      names(tmp) <- unique(d.annots[d.annots[,grpType]==grp,"group"])
      group.cols <- c(group.cols,tmp)
      i <- i+1
    }
    group.cols
  },
  
  listCompare = function(l1,l2) {
    return(list("common"=intersect(l1,l2),"l1only"=setdiff(l1,l2),"l2only"=setdiff(l2,l1)))
  },
  
  tic = function() {
    lib$ticTime <- proc.time()
  },
  
  toc = function() {
    deltaTime <- proc.time() - lib$ticTime
    lib$statusMessage(sprintf("Time elapsed: user:%.3f system:%.3f elapsed:%.3f", deltaTime[1], deltaTime[2], deltaTime[3]))
  },
  
  ticTime = proc.time(),
  
  dist = function(x,y=NULL,method="euclidean",...) {
    d <- NULL
    if(method %in% c("euclidean","maximum","manhattan","canberra","binary","minkowski")) {
      lib$loadLibrary("fastcluster") # optional speed-up for all standard dist functions
      d <- dist(x, method=method, ...)
    }
    else if(method %in% c("bray-curtis","jaccard","difference","sorensen","gower","mod-gower10","modgower2")) { #"mahalanobis",--> inefficient, use other implementation from "vegan" package
      lib$loadLibrary("ecodist")
      d <- ecodist::distance(x, method=method, ...)
    }
    else if(method %in% c("cor","pearson","spearman","kendall","cosine","mcd","ogk")) {
      lib$loadLibrary("MKmisc")
      if(method=="cor") method = "pearson"
      d <- MKmisc::corDist(x, method=method, ...)
    }
    else if(method %in% c("dtw")) {
      lib$loadLibrary("dtw")
      if(method=="cor") method = "pearson"
      d <- dtw::dtwDist(x, ...)
    }
    else if(method %in% c("kulczynski","morisita","horn","mountford","raup","binomial","chao","cao","mahalanobis")) {
      lib$loadLibrary("vegan")
      d <- vegan::vegdist(x, method=method, ...)
    }
    else if(method %in% c("mi","kld")) {
      lib$loadLibrary("bioDist")
      if(method=="mi") d <- MIdist(x)
      else if(method=="kld") d <- KLD.matrix(x)
    }
    else if(method %in% c("krls")) {
      lib$loadLibrary("KRLS")
      d <- gausskernel(x, ...)
    }
    else if(method %in% c("mse","mae","rmse","rmsle","ae","auc","ce")) {
      lib$loadLibrary("Metrics")
      if(is.null(y)) {
        y <- x[2,]
        x <- x[1,]
      }
      if(method=="mse") {
        d <- mse(x,y)
      }
      else if(method=="mae") {
        d <- mae(x,y)
      }
      else if(method=="rmse") {
        d <- rmse(x,y)
      }
      else if(method=="rmsle") {
        d <- rmse(x,y)
      }
      else if(method=="ae") {
        d <- ae(x,y)
      }
      else if(method=="auc") {
        d <- auc(x,y)
      }
      else if(method=="ce") {
        d <- ce(x,y)
      }
    }
    else lib$statusMessage("Warning: Unknown distance metric:",method)
    
    return(d)
  },
  
  prettyIntegers = function(x) {
    formatC(x, big.mark = ",", format = "d")
  },
  
  prettyDoubles = function(x,digits=2) {
    sprintf(paste0("%.",digits,"f"),x)
  },
  
  curateRegionDb = function(regionDB) {
    capFirst <- function(str) paste0(toupper(substr(str, 1, 1)), substr(str, 2, nchar(str)))
    regionDB$regionAnno[,cellType:=capFirst(gsub("(cell|progenitor|precursor|phage|cyte|blast)s","\\1", tolower(cellType), perl=TRUE))]
    regionDB
  },
  
  augmentLolaRes = function(lolaResAll, qThresh=0.05, qvalCol="qvalMean", pvalCol="pval", orCol="logOddsRatio") {
    lolaResAll[, padj:=p.adjust(get(pvalCol), method="fdr")]
    lolaResAll[,log2odds:=log2(get(orCol))]
    
    lolaResAll[antibody%in%c("GR","NR3C1"), antibody:="NR3C1"]
    lolaResAll[antibody%in%c("GABP"), antibody:="GABPA"]
    
    lolaResAll[,term:=antibody]
    lolaResAll[is.na(term),term:=gsub("_\\(.+\\)$","",gsub("GSM\\d+_","",gsub("Human_","",gsub("wgEncode.wg","",gsub(".(bed|narrowPeak)","",filename)))))]
    lolaResAll[collection=="sheffield_dnase", term:=paste0("DNase #",gsub(".bed","",filename), " (", sapply(strsplit(description,";"),function(x) paste(substr(x,1,3),collapse=";")), ")")]
    lolaResAll[,term:=gsub("^(EGFP|C)\\W","",gsub("_\\(.+\\)$","",toupper(term)))]
    
    lolaResAll[,term:=gsub("E?P300", "EP300", term, ignore.case=T)]
    lolaResAll[,term:=gsub("PU\\.?1", "SPI1", term, ignore.case=T)]
    lolaResAll[,term:=gsub("[--]", "", term)]
    lolaResAll[,term:=gsub("POL(II|LL)", "POL2", term, ignore.case=T)]
    lolaResAll[,term:=gsub("POLIII", "POL3", term, ignore.case=T)]
    lolaResAll[term=="ERALPHA_A",term:="ESR1"]
    lolaResAll[term=="TCF7L2_C9B9",term:="TCF7L2"]
    
    lolaResAll
  },
  
  getContrastingColours = function(n,theme="Set1") {
    lib$loadLibrary("RColorBrewer")
    
    if(n <= 8 & n >= 3)	{
      brewer.pal(n,theme)
    }
    else {
      rainbow(n)
    }
  },
  
  getDivergingColourScale = function(values,r=values,colours=c(rgb(0.230, 0.299, 0.754),rgb(0.865, 0.865, 0.865),rgb(0.706, 0.016, 0.150)),interval=0.8) {
    extr <- max(abs(r))*interval
    v <- (values + extr)/(2*extr)
    v <- pmax(0,pmin(1,v))
    x <- colorRamp(colours)(v)
    grDevices::rgb(x[,1], x[,2], x[,3], maxColorValue = 255)
  },
  
  # Multiple plot function
  #
  # ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
  # - cols:   Number of columns in layout
  # - layout: A matrix specifying the layout. If present, 'cols' is ignored.
  #
  # If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
  # then plot 1 will go in the upper left, 2 will go in the upper right, and
  # 3 will go all the way across the bottom.
  # http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
  multiplot = function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    lib$loadLibrary("grid")
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
      # Make the panel
      # ncol: Number of columns of plots
      # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                       ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
      print(plots[[1]])
    } 
    else {
      # Set up the page
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
      
      # Make each plot, in the correct location
      for (i in 1:numPlots) {
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
        
        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                        layout.pos.col = matchidx$col))
      }
    }
  },
  
  ggTheme = function(p,flipX = FALSE, boldTitles=FALSE, withAxisLines=FALSE, withGrid = FALSE, topLegend=FALSE, noLegendTitle=FALSE, baseSize = 12, baseFamily = "") {
    lib$loadLibrary("ggplot2")
    
    ggTheme <- theme_bw(base_size = baseSize, base_family = baseFamily) %+replace% 
      theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), 
            panel.grid.major = element_blank(), axis.text=element_text(size=baseSize-1), strip.text = element_text(face=ifelse(boldTitles,"bold","plain"),size=max(baseSize-2,baseSize*0.9)), axis.title = element_text(face=ifelse(boldTitles,"bold","plain"),size=max(baseSize-2,baseSize*0.9)),
            strip.background = element_blank(), legend.key = element_blank())
    if(!withGrid) {
      ggTheme <- ggTheme %+replace% theme(panel.grid.minor = element_blank())
    }
    if(flipX) {
      ggTheme <- ggTheme %+replace% theme(axis.text.x = element_text(angle=90,hjust=1))
    }
    if(withAxisLines) {
      ggTheme <- ggTheme %+replace% theme(axis.line=element_line(), axis.line.x=element_line(), axis.line.y=element_line())
    }
    
    if(topLegend) {
      ggTheme <- ggTheme %+replace% theme(legend.position="top")
    }
    if(noLegendTitle) {
      ggTheme <- ggTheme %+replace% theme(legend.title=element_blank())
    }
    
    return(ggTheme)
  },
  
  ggAddPanelBase = function(p, expand0 = FALSE) {
    p <- p +  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)
    if(expand0) p <- p + scale_y_continuous(expand=c(0,0))
    return(p)
    
    #return(p + expand_limits(y=0) + geom_hline(yintercept=0))
  },
  
  paste_ = function(...) {
    paste(..., sep="_")
  },
  
  mixedSort = function(x, charFirst=TRUE, ...) {
    charX <- as.character(gsub("[\\d\\.]","",x,perl=TRUE)) 
    numX <- as.numeric(gsub("[^\\d\\.]","",x,perl=TRUE))
    if(charFirst) x[order(charX, numX, ...)]
    else x[order(numX, charX, ...)]
  },
  
  # maps new data into a predefined PCA space
  projectPCA = function(pcaResult, newData) {
    # see http://stats.stackexchange.com/questions/2592/how-to-project-a-new-vector-onto-pca-space
    scale(newData, pcaResult$center, pcaResult$scale) %*% pcaResult$rotation
  },
  
  findHull = function(df, dimName1="PC1", dimName2="PC2") {
    #http://stats.stackexchange.com/questions/22805/how-to-draw-neat-polygons-around-scatterplot-regions-in-ggplot2
    df[chull(df[,dimName1], df[,dimName2]), ] 
  },
  
  # Nathan's function:
  # function to return matrix of memory consumption
  objectSizes = function() {
    gc();
    message("Objects in MB:")
    objects = ls(envir=.GlobalEnv);
    classes = sapply(objects, function(object) { class(get(object))[1] });
    sizes =  sapply(objects, function (object) { object.size(get(object)) } )
    a = data.frame(MB=sizes/1e6, class=classes)
    ord = order(sizes, decreasing=TRUE)
    a2 = a[ord,];
    a2 = a2[! a2$class == "function",]
    print(head(a2, 30))
    message("Sum: ", signif(sum(sizes/1e6),4), " MB (", signif(sum(sizes/1e9),4), " GB)")
  },	
  
  forEach = function(...) {
    void <- sapply(...)
  },
  
  nullToNA = function(x) {
    ifelse(is.null(x), NA, x)
  },
  
  makeRanges = function(dtRegs) {
    with(dtRegs, GRanges(seqnames = chrom, ranges = IRanges(start = start, end = end)))
  },
  
  aggregateSignal = function(sig, regs, signalName="count", index=TRUE, len=50) {
    #statusMessage("Aggregate signal..")
    sigGr <- NULL
    if("end"%in%colnames(sig)) {
      sigGr <- with(as.data.frame(sig), GRanges(seqnames=chrom, ranges=IRanges(start=start, end=end), strand=rep("*", nrow(sig)), intensity=get(signalName)))
    }
    else {
      sigGr <- with(as.data.frame(sig), GRanges(seqnames=chrom, ranges=IRanges(start=start, end=start+50), strand=rep("*", nrow(sig)), intensity=get(signalName)))
    }
    
    refGr <- with(as.data.frame(regs), GRanges(seqnames=chrom, ranges=IRanges(start=start, end=end), strand=rep("*", nrow(regs))))
    
    #statusMessage("Finding overlaps...")
    fo <- findOverlaps(sigGr, refGr)
    setkey(sig, chrom, start)
    
    if(index) {
      #statusMessage("Indexing signal...")
      sig[,inSet:=FALSE] 
      sig[unique(queryHits(fo)),inSet:=TRUE] 
      return(sig)
    }
    else {
      #statusMessage("Aggregating signal...")
      sig <- sig[queryHits(fo),] 
      sig[,regionID:=subjectHits(fo)] 
      maxC <- sig[, list(medianX=median(get(signalName),na.rm=T),meanX=mean(get(signalName),na.rm=T),maxX=max(get(signalName),na.rm=T),minX=min(get(signalName),na.rm=T)), by=list(sampleName,regionID)]
      #maxC <- maxC[order(regionID),]
      return(maxC)
    }
  },
  
  aggregateSignalCustom = function(sig, regs, aggFun=function(x) mean(x, na.rm=TRUE), signalName="count", index=TRUE, len=50) {
    #statusMessage("Aggregate signal..")
    sigGr <- NULL
    if("end"%in%colnames(sig)) {
      sigGr <- with(as.data.frame(sig), GRanges(seqnames=chrom, ranges=IRanges(start=start, end=end), strand=rep("*", nrow(sig)), intensity=get(signalName)))
    }
    else {
      sigGr <- with(as.data.frame(sig), GRanges(seqnames=chrom, ranges=IRanges(start=start, end=start+50), strand=rep("*", nrow(sig)), intensity=get(signalName)))
    }
    
    refGr <- with(as.data.frame(regs), GRanges(seqnames=chrom, ranges=IRanges(start=start, end=end), strand=rep("*", nrow(regs))))
    
    #statusMessage("Finding overlaps...")
    fo <- findOverlaps(sigGr, refGr)
    setkey(sig, chrom, start)
    
    msg("sigGr ", length(sigGr))
    msg("refGr ", length(refGr))
    msg("q hits ", length(queryHits(fo)))
    msg("sub hits ", length(subjectHits(fo)))
    
    if(index) {
      #statusMessage("Indexing signal...")
      sig[,inSet:=FALSE] 
      sig[unique(queryHits(fo)),inSet:=TRUE] 
      return(sig)
    }
    else {
      #statusMessage("Aggregating signal...")
      sig <- sig[queryHits(fo),] 
      sig[,regionID:=subjectHits(fo)] 
      maxC <- sig[, list(aggVal=aggFun(get(signalName))), by=list(sampleName,regionID)]
      #maxC <- maxC[order(regionID),]
      return(maxC)
    }
  },
  
  capFirst = function(str) paste0(toupper(substr(str, 1, 1)), substr(str, 2, nchar(str))),
  
  pickColors = function(dA) {
    dfA <- as.data.frame(dA)
    sapply(colnames(dfA),function(curCol) {
      lib$getCategoryColors(unique(dfA[,curCol]))
    }, simplify=FALSE)
  },
  
  dtToDf = function(dt, rownameCol=1) {
    df <- as.data.frame(dt)
    rownames(df) <- df[,rownameCol]
    if(is.numeric(rownameCol)) {
      df <- df[,-rownameCol,drop=FALSE] 
    }
    else {
      df <- df[,setdiff(colnames(df),rownameCol),drop=FALSE] 
    }
    df
  },
  
  blendColors = function(cols, alpha=100) {
    blended <- round(rowSums(sapply(cols, col2rgb))/length(cols))
    rgb(blended[1], blended[2], blended[3], alpha=alpha, maxColorValue=255)
  },	
  
  naToZero = function(d) {
    sapply(d, function(x) {
      if(is.null(x) | is.na(x) | length(x)==0) return(0)
      else return(x)
    })
  },
  
  abscap = function(x, cap=0.99) {
    thresh <- quantile(abs(x), cap, na.rm=T)
    i <- !is.na(x) & abs(x)>=thresh
    x[i] <- sign(x[i]) * thresh
    x
  },
  
  absmax = function(x) {
    x[which.max(abs(x))]
  },
  
  absmin = function(x) {
    x[which.min(abs(x))]
  },
  
  pToSig = function(p, ns="", threshs=c(0.05, 0.01, 0.001)) {
    sapply(p, function(pval) ifelse(pval>threshs[1], ns, ifelse(pval<=threshs[3], "***", ifelse(pval<=threshs[2], "**", "*"))) )
  },
  
  runEnrichr = function(geneLists, maxSize=2000, databases=c("GO_Biological_Process_2015")) {		
    rbindlist(sapply(names(geneLists), function(listName){
      gs <- geneLists[[listName]]
      n <- length(gs)
      
      msgF("%s (%s genes)", listName, lib$prettyIntegers(n))
      
      if(n > maxSize) {
        msgF("	--> sample to %s", lib$prettyIntegers(maxSize))
        gs <- gs[sample(n,maxSize)]
      }	
      
      
      dtTmp=try(as.data.table(lib$runEnrichrOnGeneList(gs,databases = databases)),silent = FALSE)
      if(!any(grepl("Error",dtTmp)) && nrow(dtTmp) > 0){
        return(dtTmp)
      }
      print(dtTmp)
      
      return(NULL)
      
    }, simplify=FALSE), idcol="grp")
  },
  
  # modified version of original Enrichr function to return additional fields:
  runEnrichrOnGeneList = function (gene.list, databases = db.list, fdr.cutoff = NULL) {
    req.body <- list(list = paste(gene.list, collapse = "\n"))
    post.req <- httr::POST("http://amp.pharm.mssm.edu/Enrichr/enrich", 
                           encode = "multipart", body = I(req.body))
    if (!grepl("success", httr::http_status(post.req)$category, 
               ignore.case = T)) 
      stop("Posting gene list to EnrichR failed")
    database.enrichments <- list()
    for (idx in 1:length(databases)) {
      database <- databases[idx]
      get.req <- httr::GET(paste("http://amp.pharm.mssm.edu/Enrichr/enrich?backgroundType=", 
                                 database, sep = ""))
      if (!grepl("success", httr::http_status(get.req)$category, 
                 ignore.case = T)) 
        stop("Retrieving results from EnrichR failed")
      response.content <- enrichR:::mungeResponseContent(httr::content(get.req)[[database]])
      if (length(response.content) > 1) {
        database.res <- data.table::rbindlist(response.content)
        database.res[, 1] <- rep(database, nrow(database.res))
        database.enrichments[[idx]] <- database.res[, paste("V", 
                                                            c(1, 2, 3, 4, 5, 7, 6), sep = ""), with = F] # NEEDED TO MODIFY THIS SO THAT ALL COLUMNS OF THE RESPONSE ARE STORED!
      }
    }
    query.results <- as.data.frame(data.table::rbindlist(database.enrichments))
    colnames(query.results) <- c("database", "category", "pval", 
                                 "zscore", "combscore", "qval", "genes")
    if (!is.null(fdr.cutoff)) {
      query.results <- query.results[query.results$qval < fdr.cutoff, ]
    }
    return(query.results)
  },
  
  sanitizeEnrichrLabels = function(txt) {
    sapply(txt, function(x) {
      x <-  gsub(" (Mus musculus|Mouse)", " (mouse)", gsub(" (Homo sapiens|Human)", " (human)", gsub("_", " ", x), ignore.case=TRUE), ignore.case=TRUE)
      x <- gsub(" \\(NAFLD\\)", "", x)
      x <- gsub("^MP\\d+\\s+", "", x, perl=TRUE)
      x <- gsub("\\s+(hsa|WP)\\d+$", "", x, perl=TRUE)
      x <- gsub("\\s+\\(GO:.+\\)$", "", x, perl=TRUE)
      x <- gsub(" \\d+ ChIP-.+ ", " ", x, perl=TRUE)
      x <- gsub("positive", "pos.", x, perl=TRUE)
      x <- gsub("negative", "neg.", x, perl=TRUE)
      x <- gsub("regulation of", "reg. of", x, perl=TRUE)
      x <- gsub("involved in ", "in ", x, perl=TRUE)
      x <- gsub("ligase activity", "ligase", x, perl=TRUE)		
      x <- gsub("(GSE\\d+) sample \\d+", "\\1", x, perl=TRUE)
      x <- gsub("UMLS ", "", x)
      x <- gsub("\\s+", " ", x, perl=TRUE)
      x <- gsub("in DMSO-Rat-Primary rat ", "DMSO-Rat ", x)
      x <- gsub("\\_(\\w|\\d){8}-(\\w|\\d){4}-(\\w|\\d){4}-(\\w|\\d){4}-(\\w|\\d){12}", "",x)
      x <- gsub("\\U3e30613c", " ", x)
      x <- gsub("organization.from", "organization from", x) # replace a weird character
      x <- gsub("_", " ", x)
      x <- gsub("protein[- ]*protein interactions*", "PPI", x, perl=TRUE, ignore.case=TRUE)
      x <- gsub("Human", "", x, perl=TRUE, ignore.case=TRUE)
      x <- gsub(" sample \\d+", "", x, perl=TRUE, ignore.case=TRUE)
      x <- gsub("transcription factor", "TF", x, perl=TRUE, ignore.case=TRUE)
      x <- gsub(" pathway", "", x, perl=TRUE, ignore.case=TRUE)
      x <- gsub("hsa\\d+$", "", x, perl=TRUE, ignore.case=TRUE)
      x <- gsub("Homo sapiens", "", x, perl=TRUE, ignore.case=TRUE)
      x <- gsub("mus musculus", "MM", x, perl=TRUE, ignore.case=TRUE)
      x <- gsub(".bed", "", x)
      x <- gsub("in absence of", "without", x)
      x <- gsub("treatment", "treat", x, ignore.case=TRUE)
      x <- gsub("shRNA targeting ", " ", x, ignore.case=TRUE)
      x <- gsub("transactivating \\(TA\\) domain", "", x, ignore.case=TRUE)
      x <- gsub(" (or|and) ", "/", x, ignore.case=TRUE)
      x <- gsub("(Neg|pos).tive regulation", "\\1. regulation", x, ignore.case=TRUE)
      x <- gsub("using ", "", x, ignore.case=TRUE)
      x <- gsub("knock-*down", "KD", x, ignore.case=TRUE)
      x <- gsub("knock-*out", "KO", x, ignore.case=TRUE)
      x <- gsub("antibody", "ab", x, ignore.case=TRUE)
      x <- gsub("GSM\\d+ ", "", x, ignore.case=TRUE,perl=TRUE)
      x <- gsub("hesc", "", x)
      x <- gsub("  +", " ", x)
      lib$capFirst(substr(x, 1, 64))
    })
  },
  
  augmentEnrichrResults = function(enrichrRes, adjMeth="bonferroni") {
    enrichrRes$n <- sapply(strsplit(enrichrRes$genes,","), length)
    enrichrRes[, padj:=p.adjust(pval, method=adjMeth)]
    enrichrRes[, term:=lib$sanitizeEnrichrLabels(category)]
    setnames(enrichrRes, lib$toCamelCase(gsub("[ \\-]","_",names(enrichrRes))))
    enrichrRes
  },
  
  loadTranscriptAnnotation = function(genomeBuild, nameCol="Gene name", typeCol="Transcript type", trIdCol="Transcript stable ID", gIdCol="Gene stable ID", chromCol="Chromosome/scaffold name", startCol="Gene start (bp)", endCol="Gene end (bp)", strandCol="Strand") {
    transcriptAnnotation <- unique(fread(tail(system(paste0("ls ", lib$resourceDir("genomes/",genomeBuild,"_cdna/transcripts_*")),intern=TRUE),n=1)))[get(chromCol)%in%c(1:100,"X","Y","M","MT")]
    
    transcriptAnnotation[,chrom:=paste0("chr",gsub("MT","M",get(chromCol)))]
    transcriptAnnotation[,start:=get(startCol)]
    transcriptAnnotation[,end:=get(endCol)]				
    transcriptAnnotation[,tss:=start]
    transcriptAnnotation[,strand:=get(strandCol)]
    transcriptAnnotation[strand==-1,tss:=end]
    transcriptAnnotation[,ensT:=get(trIdCol)]
    transcriptAnnotation[,ensG:=get(gIdCol)]
    transcriptAnnotation[,type:=get(typeCol)]
    transcriptAnnotation[,geneSymbol:=gsub("\\.\\d+", "", get(nameCol))]
    
    setkey(transcriptAnnotation,ensT)
    
    transcriptAnnotation
  }
  
)
message("... done!")




# quick aliases:
naToZero <- lib$naToZero
absmax <- lib$absmax
absmin <- lib$absmin
abscap <- lib$abscap
statusMessage <- lib$statusMessage
msg <- statusMessage
msgF <- lib$statusMessageFormatted
paste_ <- lib$paste_
os <- lib$objectSizes
loadLibrary <- lib$loadLibrary
loadLibraries <- lib$loadLibraries
setCurrentAnalysis <- lib$setCurrentAnalysis
lib$defaultPlotCmd <- lib$svgPlot
lib$defaultPlot <- lib$defaultPlotCmd
defPlot <- lib$defaultPlot
defTheme <- lib$ggTheme
forEach <- lib$forEach
gg <- lib$ggsave
analyze <- lib$analyze
rblapply <- lib$rblapply
splitApply <- lib$splitApply
rowApply <- lib$rowApply
run <- lib$run
runBox <- lib$runBox
scriptPath <- lib$scriptPath
dependOn <- lib$dependOn


resultsDir <- plotDir <- lib$plotDir <- lib$resultsDir
baseResultsDir <- lib$baseResultsDir
metaDir <- lib$metaDir
dataDir <- lib$dataDir
resourceDir <- lib$resourceDir
toolsDir <- lib$toolsDir

tryCatch({ 
  lib$init() 
}, error=function(e) {
  print(e)
  msg("Continuing despite error in init()!")
})


# SimpleCache wrapper that generates a lock file when it starts the process and
# deletes it afterwards. Meant as a simple way to distribute the generation of many
# caches across multiple computers (a poor man's parallel processing across machines)
createCacheWithLock <- function(n, fun) {
  if(!file.exists(paste0(simpleCache::getCacheDir(),"/",n,".RData"))) {
    lockFile <- paste0("lock__", Sys.getenv("COMPUTER"), "__", n, ".lockfile")
    msgF("Create cache '%s' with lockfile %s...", n, lockFile)
    fs <- list.files(simpleCache::getCacheDir(), pattern=paste0("lock__.+.",n,".lockfile"))
    #print(fs)
    if(length(fs) > 0) {
      msgF("Lock files exist: %s", paste(fs, collapse=", "))
      if(length(fs)>1 || fs[1]!=lockFile) {
        msgF("--> cache creation locked on other computer, skipped.")
        return(NULL)
      } 
    }
    if(length(fs) > 0) {
      msgF("--> cache creation locked on same computer, canceled process? Replacing lock...")
    }
    lockFilePath <- paste0(simpleCache::getCacheDir(),"/", lockFile)
    #print(lockFilePath)
    write.csv(Sys.time(), file=lockFilePath)
    tryCatch({
      simpleCache::simpleCache(cacheName=n, noload=T, instruction=fun)
      msgF("Cache created: %s", n)
    }, 
    error=function(e) {
      print(paste("Error:  ",e))
      NULL
    },
    finally= {
      file.remove(lockFilePath)
      msgF("Lockfile %s deleted.\n", lockFile)
    }
    )
  }
}
fwrite <- function(x, ...) data.table::fwrite(tryCatch(as.data.table(x, keep.rownames=T), error=function(e) as.data.table(as.matrix(x), keep.rownames=T)), ...)
forEachDT <- function(dt, fun) {
  invisible({ void <- FALSE; if(nrow(dt)>0) void <- apply(dt, 1, function(x) fun(as.list(x))); void <- TRUE})
}






selectFromList <- function(l, att) {
  sapply(l, function(x) x[[att]], simplify=F)
}
hclustImpute <- function(d, ...) {
  d[is.na(d)] <- max(d, na.rm=T)
  hclust(d, ...)
}



plotGridPaged <- function(pl, nc=6, nr=3) {
  nPlotsPerPage <- nc*nr
  pageNums <- floor(((1:length(pl))-1) / nPlotsPerPage)
  lapply(unique(pageNums), function(pageNum) {
    cowplot::plot_grid(plotlist=pl[pageNums==pageNum], ncol=nc, nrow=nr)
  })
}

plotMetaDataDots <- function(pData, clusterBounds, m, plotTitle=m, filtNA=T, withText=F, withGuide="legend") {
  loadLibrary("ggrastr")
  
  pDataTmp <- data.table(pData)
  if(filtNA) pDataTmp <- pDataTmp[!is.na(get(m)) & get(m)!="NA",]
  pDataTmp[, val:=as.character(get(m))]
  
  uniq <-  pDataTmp[,unique(val)]
  cols <- getColors(m, uniq, na_col=alpha("grey",0))
  if(length(uniq)>length(cols)) stop("length(uniq)>length(cols)")
  
  topN <- 1
  topPerc <- pDataTmp[, .(X=mean(X),Y=mean(Y), .N), by=.(v=val,cluster_id)][,.(X,Y, perc=N/sum(N)*100, lbl=sprintf("%s = %.1f%%", v, N/sum(N)*100), rnk=rank(-N, ties.method="random")), by=cluster_id][rnk<=topN,][order(rnk), .(X,Y,perc=max(perc),lbl=paste(lbl, collapse="\n")), by=cluster_id]
  
  #withGuide <- ifelse((length(uniq)<11 & topPerc[,sum(perc<95)]!=0) | forceLegend, "legend", FALSE)
  
  pDots <- ggplot(pDataTmp[sample(nrow(pDataTmp),min(MAX_DOTS_UMAP,nrow(pDataTmp))),], aes(x=X, y=Y, color=val)) + geom_point_rast(size=0.5) 
  pDots <- pDots + geom_polygon(aes(group=cluster_id), size=0.1, color="darkgrey", fill=NA, data=clusterBounds)
  pDots <- pDots + scale_color_manual(values=cols, guide=withGuide) 
  if(withText) pDots <- pDots + geom_text_repel(aes(label=lbl), size=2, color="darkgrey", data=topPerc)
  stylePlot(pDots, plotTitle)
}

plotDataDots <- function(pData, clusterBounds, m, plotTitle=m, colScaleDef=c("#DDDDDD",brewer.pal(5,"YlGnBu")[2:5]), fixColorBreaks=NULL) {
  loadLibrary("ggrastr")
  
  pDataTmp <- data.table(pData)[!is.na(get(m)) & get(m)!="NA",]
  pDataTmp[, val:=as.numeric(get(m))]
  
  cols <- colScaleDef
  if(m%in%names(colorPalettes)) cols <- colorPalettes[[m]]
  cols <- alpha(cols,0.75)
  
  if(!is.null(fixColorBreaks)) {
    pDataTmp[val<min(fixColorBreaks), val:=min(fixColorBreaks)]
    pDataTmp[val>max(fixColorBreaks), val:=max(fixColorBreaks)]
    
    cols <- colorRampPalette(cols)(length(fixColorBreaks)-1)
    for(i in 1:(length(fixColorBreaks)-1)) {
      pDataTmp[val>=fixColorBreaks[i] & val<fixColorBreaks[i+1], mbin:=i]
      names(cols)[i] <- sprintf("[%.1f,%.1f[", fixColorBreaks[i], fixColorBreaks[i+1])
    }
    pDataTmp[val==fixColorBreaks[i+1], mbin:=i]
    names(cols)[i] <- sprintf("[%.1f,%.1f]", fixColorBreaks[i], fixColorBreaks[i+1])
    
    pDataTmp[, val:=factor(mbin, levels=1:(length(fixColorBreaks)-1), labels=names(cols))]
    
    cols <- c(cols, "NA"=gray)
  }
  
  pDots <- ggplot(pDataTmp[sample(nrow(pDataTmp),min(MAX_DOTS_UMAP,nrow(pDataTmp))),], aes(x=X, y=Y, color=val)) + geom_point_rast(size=0.5) 
  pDots <- pDots + geom_polygon(aes(group=cluster_id), color="darkgrey", fill=NA, size=0.1, data=clusterBounds)
  if(is.null(fixColorBreaks)) {
    pDots <- pDots + scale_color_gradientn(colors=cols, guide=guide_colorbar(ticks.color="black", barwidth=0.5, barheight=4)) 
  } else {	
    pDots <- pDots + scale_color_manual(values=cols)
  }
  
  stylePlot(pDots, plotTitle, border=F)
}

MAX_DOTS_UMAP <- 50000

legendPlot <- function(colorContinuous=NULL, colorCategorical=NULL) {
  p <- ggplot()+theme_void()
  if(!is.null(colorContinuous)) {
    colorData <- data.table(c(colorContinuous$from, colorContinuous$to))
    setnames(colorData, colorContinuous$title)
    p <- p + geom_point(x=NA,y=NA,aes_string(color=colorContinuous$title), size=3, data=colorData) + scale_color_gradientn(colors=colorContinuous$cols, guide=guide_colorbar())
  }	
  if(!is.null(colorCategorical)) {
    colorData <- data.table(names(colorCategorical$values))
    setnames(colorData, colorCategorical$title)
    p <- p + geom_point(x=NA,y=NA,aes_string(color=colorCategorical$title), size=3, data=colorData) + scale_color_manual(values=colorCategorical$values, guide=guide_legend(ncol=ceiling(length(colorCategorical$values)/9)))
  }	
  p
}
textPlot <- function(txt, flip=F) {
  ggplot()+theme_void()+ggtitle(txt)+theme(plot.title=element_text(angle=ifelse(flip,90,0)))
}
emptyPlot <- function(txt) {
  ggplot()+theme_void()
}
stylePlot <- function(p, plotTitle, border=T) {	
  p + ggtitle(plotTitle) +
    defTheme() +
    theme(
      legend.position = c(1, 0),
      legend.title=element_blank(),
      legend.justification = c("right", "bottom"),
      legend.box.just = "right",
      legend.box.background = element_rect(color=ifelse(border,"black",NA), size=0.5),
      legend.box.spacing = unit(1,"pt"),
      legend.box.margin = ggplot2::margin(2, 2, 2, 2),
      legend.key.size = unit(3,"pt"),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      axis.text.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.ticks.y=element_blank(), 
      validate = TRUE
    ) +
    scale_y_continuous(expand=expansion(mult=0.1)) +
    scale_x_continuous(expand=expansion(mult=0.1))
}

# dataDir <- function(...) {
# paste0(Sys.getenv("DATA"),"/os_tumors/",...)
# }
getColors <- function(type, values, sorted=TRUE, na_col=alpha("grey",0.25)) {
  values <- unique(values)
  if(sorted) values <- sort(values)
  if(type%in%names(colorPalettes)) {
    cols <- colorPalettes[[type]]
  } else {
    cols <- lib$getCategoryColors(values)
  }
  cols[["NA"]] <- na_col
  cols[values]
}
getGeoSupp <- function(gsm, baseDir=resultsDir("geo"), filt=NULL) {
  loadLibrary("GEOquery")
  dir.create(baseDir, showWarnings=F)	
  supp <- as.data.table(getGEOSuppFiles(gsm, fetch=F, filter_regex=filt))
  supp[, fpath:=paste0(baseDir,"/",fname)]
  forEachDT(supp[!file.exists(fpath),], function(x) {
    msg("download ", x$url)
    download.file(x$url, destfile=x$fpath)
  })	
  return(supp)
}

