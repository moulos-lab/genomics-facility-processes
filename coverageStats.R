# This script calculates, for each provided region target:
# - Minimum coverage over target
# - Maximum coverage over target
# - Mean coverage over target
# - Median coverage over target
# - Total stats (but this requires more thinking, maybe a different sheet)
#   as we need to add min, max, etc apart from mean
#
# seqlevelStyles must also be taken into account
#
# This script should could be used along with HybridStat's exon coverage script
# for targetted sequencing and bamstats.
#
coverageStats <- function(
    bams,
    targets=NULL,
    controls=NULL,
    fragLen=300,
    asObject=FALSE,
    pairedOpts=list(isPaired=NULL,asPaired=FALSE,ifPaired="list"),
    outFormat=c("xlsx","txt"),
    avgSamples=FALSE,
    restricBamToTargets=FALSE,
    rc=NULL
) {
    if (missing(bams))
        stop("One or more BAM files must be provided!")
    # TODO: Add check if any file does not exist
    
    if (is.null(targets))
        warning("Targets not provided! Total coverage will be calculated",
            immediate.=TRUE)
    else if (!(is.character(targets) && file.exists(targets)))
        stop("The targets file must be a valid file!")
        
    if (!is.null(controls)) {
        if (!(is.character(controls) && file.exists(controls))) {
            warning("Control regions must be a valid file! Ignoring...",
                immediate.=TRUE)
            controls <- NULL
        }
    }
        
    if (!require(GenomicAlignments))
        stop("Bioconductor package GenomicAlignments is required!")
    if (!require(GenomicRanges))
        stop("Bioconductor package GenomicAlignments is required!")
        
    outFormat <- outFormat[1]
    if (!is.character(outFormat) || !(outFormat %in% c("txt","xlsx")))
        stop("outFormat must be one of \"txt\", \"xlsx\"")
    if (outFormat == "xlsx" && !require(openxlsx))
        stop("R package openxlsx is required for Excel output")
        
    if (!is.numeric(fragLen) || fragLen < 0)
        stop("fragLen must be an integer > 0!")
        
    if (!is.null(rc) && rc==0) # Convenience for accepting from command line
        rc <- NULL # Use 1 core if rc=0
       
    # Read targets - must be a valid bed file!
    if (!is.null(targets)) {
        message("Reading targets file ",targets)
        preBed <- read.table(targets,header=FALSE)
        if (is.null(preBed$V4)) {
            names(preBed)[1:3] <- c("chromosome","start","end")
            bed <- GRanges(preBed)
            names(bed) <- paste(seqnames(bed),":",start(bed),"-",end(bed),
                sep="")
        }
        else {
            names(preBed)[1:4] <- c("chromosome","start","end","name")
            bed <- GRanges(preBed)
            names(bed) <- bed$name
        }
    }
    
    # Read controls - a simple 1-column text file without header!
    if (!is.null(controls)) {
        message("Reading controls file ",controls)
        ctrls <- as.character(read.table(controls,header=FALSE)[,1])
    }
    else
        ctrls <- NULL
    
    # If less than 2 files, sample averaging is meaningless...
    if (avgSamples && length(bams)==1) {
        warning("Sample averaging is not meaningful for one sample! ",
            "Ignoring...",immediate.=TRUE)
        avgSamples <- FALSE
    }
    
    # Import reads supporting the targets
    bamRanges <- NULL
    if (!is.null(targets) && restricBamToTargets) {
        if (!is.null(fragLen)) {
            message("Expanding target locations to enclose reads")
            w <- width(bed)
            bamRanges <- suppressWarnings(promoters(bed,upstream=fragLen,
                downstream=0))
            bamRanges <- resize(bamRanges,width=w+2*fragLen)
        }
        else
            bamRanges <- bed
    }
    
    att <- NULL
    if (!is.null(targets))
        att <- preBed[,seq_len(3)]
    
    #TODO: Pass all internal parameters to parallel function, now reading
    #      from upper level which is not very safe.
    results <- cmclapply(basename(bams),function(b) {
        message("Sample ",b)
        message("  reading BAM file")
        if (is.null(targets))
            reads <- readBamAlignments(b)
        else {
            if (restricBamToTargets)
                reads <- readBamAlignments(b,gr=bamRanges)
            else {
                allReads <- readBamAlignments(b,asRanges=TRUE)
                ov <- findOverlaps(allReads,bed,ignore.strand=TRUE)
                reads <- allReads[queryHits(ov)]
            }
        }
        # TODO: Pass the other supported parameters to readBamAlignments
        
        # Coverage operations
        message("  calculating coverage stats")
        co <- calcCoverage(reads,bed,rc=rc)
        
        # Min coverage
        minCo <- vapply(co,min,numeric(1))
        # Max coverage
        maxCo <- vapply(co,max,numeric(1))
        # Mean coverage
        meanCo <- vapply(co,mean,numeric(1))
        # Median coverage
        medianCo <- vapply(co,median,numeric(1))
        # First base coverage
        firstCo <- vapply(co,function(x) {
            return(as.numeric(x[1]))
        },numeric(1))
        # Last base coverage
        lastCo <- vapply(co,function(x) {
            return(as.numeric(x[length(x)]))
        },numeric(1))
        
        # Assemble all stats
        D <- data.frame(
            min_coverage=minCo,
            max_coverage=maxCo,
            mean_coverage=meanCo,
            median_coverage=medianCo,
            first_base_coverage=firstCo,
            last_base_coverage=lastCo,
            row.names=names(co)
        )
        
        # Explicit names
        D <- cbind(name=rownames(D),D)
        
        if (!is.null(att))
            D <- cbind(att,D)
        
        # If controls, then separate them
        if (!is.null(ctrls)) {
            # Safety...
            if (!all(ctrls %in% rownames(D)))
                ctrls <- intersect(rownames(D),ctrls)
            ampl <- setdiff(rownames(D),ctrls)
            Dd <- D[ampl,]
            Dc <- D[ctrls,]
        }
        else {
            Dd <- D
            Dc <- NULL
        }
        
        # Create total stats
        totalD <- data.frame(
            name="total_average_targets",
            average_min_coverage=mean(min(Dd$min_coverage)),
            average_max_coverage=mean(max(Dd$max_coverage)),
            average_mean_coverage=mean(Dd$mean_coverage),
            average_median_coverage=mean(Dd$median_coverage),
            average_first_base_coverage=mean(Dd$first_base_coverage),
            average_last_base_coverage=mean(Dd$last_base_coverage)
        )
        if (!is.null(Dc))
            totalC <- data.frame(
                name="total_average_controls",
                average_min_coverage=mean(min(Dc$min_coverage)),
                average_max_coverage=mean(max(Dc$max_coverage)),
                average_mean_coverage=mean(Dc$mean_coverage),
                average_median_coverage=mean(Dc$median_coverage),
                average_first_base_coverage=mean(Dc$first_base_coverage),
                average_last_base_coverage=mean(Dc$last_base_coverage)
            )
        else
            totalC <- NULL
        Dt <- rbind(totalD,totalC)

        return(list(targets=Dd,controls=Dc,totals=Dt))
    },rc=rc)
    #names(results) <- gsub(".bam","",basename(bams),ignore.case=TRUE)
    names(results) <- basename(bams)
       
    if (avgSamples) {
        avgStatsTargets <- .averageAcrossSamples(results,"targets")
        rownames(avgStatsTargets) <- rownames(results[[1]]$targets)
        if (!is.null(ctrls)) {
            avgStatsControls <- .averageAcrossSamples(results,"controls")
            rownames(avgStatsControls) <- rownames(results[[1]]$controls)
        }
        avgStatsTotals <- .averageAcrossSamples(results,"totals")
        rownames(avgStatsTotals) <- rownames(results[[1]]$totals)
        
        # Attach coordinates etc
        avgStatsTargets <- 
            cbind(results[[1]]$targets[,c(1,2,3)],avgStatsTargets)
        if (!is.null(ctrls)) {
            avgStatsControls <- 
                cbind(results[[1]]$controls[,c(1,2,3)],avgStatsControls)
    }
    
    if (outFormat == "txt") {
        message(paste0("Writing summary output to separate text files, one ",
            "group of files per sample"))
        for (b in bams) {
            base <- basename(b)
            parts <- strsplit(base,"\\.")[[1]]
            if (length(parts) == 1) { # No extension
                outTargets <- paste(parts,"_covstats_TARGETS.txt",sep="")
                if (!is.null(ctrls))
                    outControls <- paste(parts,"_covstats_CONTROLS.txt",sep="")
                outTotals <- paste(parts,"_covstats_TOTALS.txt",sep="")
            }
            else {
                outTargets <- paste(paste(parts[1:(length(parts)-1)],
                    collapse="."),"_covstats_TARGETS.txt",sep="")
                if (!is.null(ctrls))
                    outControls <- paste(paste(parts[1:(length(parts)-1)],
                        collapse="."),"_covstats_CONTROLS.txt",sep="")
                outTotals <- paste(paste(parts[1:(length(parts)-1)],
                    collapse="."),"_covstats_TOTALS.txt",sep="")
            }
            
            write.table(results[[base]]$targets,
                file=file.path(dirname(b),outTargets),sep="\t",quote=FALSE,
                row.names=FALSE)
            if (!is.null(ctrls))
                write.table(results[[base]]$controls,
                    file=file.path(dirname(b),outControls),sep="\t",quote=FALSE,
                    row.names=FALSE)
            write.table(results[[base]]$totals,
                file=file.path(dirname(b),outTotals),sep="\t",quote=FALSE,
                row.names=FALSE)
        }
        
        if (avgSample) {
            outTargets <- file.path(dirname(bams[1]),
                paste("average_covstats_summary_",
                format(Sys.time(),"%Y-%m-%d-%H-%M-%S"),"_TARGETS.txt",sep=""))
            if (!is.null(ctrls))
                outControls <- file.path(dirname(bams[1]),
                    paste("average_covstats_summary_",
                    format(Sys.time(),"%Y-%m-%d-%H-%M-%S"),"_CONTROLS.txt",
                    sep=""))
            outTotals <- file.path(dirname(bams[1]),
                paste("average_covstats_summary_",
                format(Sys.time(),"%Y-%m-%d-%H-%M-%S"),"_TOTALS.txt",sep=""))
                
            write.table(avgStatsTargets,file=outTargets,sep="\t",quote=FALSE,
                row.names=FALSE)
            if (!is.null(ctrls))
                write.table(avgStatsControls,file=outControls,sep="\t",
                    quote=FALSE,row.names=FALSE)
            write.table(avgStatsTotals,file=outTotals,sep="\t",quote=FALSE,
                row.names=FALSE)
        }
    }
    else if (outFormat == "xlsx") {
        message("Writing summary output to one spreadsheet file per sample")
        for (b in bams) {
            base <- basename(b)
            parts <- strsplit(base,"\\.")[[1]]
            if (length(parts) == 1) # No extension
                out <- paste(parts,"_covstats.xlsx",sep="")
            else
                out <- paste(paste(parts[1:(length(parts)-1)],collapse="."),
                    "_covstats.xlsx",sep="")
            out <- file.path(dirname(b),out)
            if (!is.null(ctrls))
                xlsx <- list(
                    targets=results[[base]]$targets,
                    controls=results[[base]]$controls,
                    totals=results[[base]]$totals
                )
            else
                xlsx <- list(
                    targets=results[[base]]$targets,
                    totals=results[[base]]$totals
                )
            write.xlsx(xlsx,file=out,keepNA=TRUE)
        }
        
        if (avgSample) {
            out <- file.path(dirname(bams[1]),paste("average_covstats_summary_",
                format(Sys.time(),"%Y-%m-%d-%H-%M-%S"),".xlsx",sep=""))
            out <- file.path(dirname(b),out)
            if (!is.null(ctrls))
                xlsx <- list(
                    targets=avgStatsTargets,
                    controls=avgStatsControls,
                    totals=avgStatsTotals
                )
            else
                xlsx <- list(
                    targets=avgStatsTargets,
                    totals=avgStatsTotals
                )
            write.xlsx(xlsx,file=out,keepNA=TRUE)
        }
    }
    # TODO: Define output dir as argument
    
    if (asObject)
        return(results)
}

readBamAlignments <- function(bam,gr=NULL,params=NULL,isPaired=NULL,
    asPaired=FALSE,ifPaired=c("pairs","list"),seq=FALSE,asRanges=FALSE) {
    ifPaired <- tolower(ifPaired[1])
    if (!(ifPaired %in% c("pairs","list")))
        stop("ifPaired must be one of \"pairs\" or \"pairs\"")
    
    bamIndex <- paste(bam,"bai",sep=".")
    if (!file.exists(paste0(bam,".bai")))
        prepareBam(bam)
    if (is.null(isPaired) && asPaired)
        isPaired <- testPairedEndBam(bam)
    else
        isPaired <- FALSE
    
    if (is.null(params)) {
        if (!is.null(gr)) {
            if (seq)
                params <- ScanBamParam(
                    flag=scanBamFlag(
                        isSecondaryAlignment=FALSE,
                        isSupplementaryAlignment=FALSE
                    ),
                    which=gr,
                    what="seq"
                )
            else
                params <- ScanBamParam(
                    flag=scanBamFlag(
                        isSecondaryAlignment=FALSE,
                        isSupplementaryAlignment=FALSE
                    ),
                    which=gr
                )
        }
        else {
            if (seq)
                params <- ScanBamParam(
                    flag=scanBamFlag(
                        isSecondaryAlignment=FALSE,
                        isSupplementaryAlignment=FALSE
                    ),
                    what="seq"
                )
        }
    }
    else {
        if (!is.null(gr)) {
            if (seq) {
                bamWhich(params) <- gr
                bamWhat(params) <- "seq"
            }
            else
                bamWhich(params) <- gr
        }
        else {
            if (seq)
                bamWhat(params) <- "seq"
        }
    }
    
    if (isPaired && asPaired) {
        bamCon <- BamFile(file=bam,asMates=TRUE)
        if (ifPaired == "pairs")
            # TODO: Control strandMode
            reads <- readGAlignmentPairs(file=bam,index=bamIndex,
                param=params,with.which_label=TRUE,strandMode=1)
        else if (ifPaired == "list")
            reads <- readGAlignmentsList(file=bamCon,index=bamIndex,
                param=params,with.which_label=TRUE)
    }
    else {
        bamCon <- BamFile(file=bam)
        #reads <- readGAlignments(file=bamCon,index=bamIndex,
        #    param=params,with.which_label=TRUE)
        reads <- readGAlignments(file=bamCon,index=bamIndex,param=params)
    }
    
    #Error in `colnames<-`(`*tmp*`, value = colnames) :
    #  more column names than columns
    #In addition: Warning message:
    #In Rsamtools:::.load_bamcols_from_scanBam_res(res, param, with.which_label = #with.which_label) :
    #  'which_label' is ignored when 'param' is missing or doesn't have a 'which' component

    if (asRanges)
        return(as(reads,"GRanges"))
    else
        return(reads)
}

calcCoverage <- function(input,mask,strand=NULL,ignore.strand=TRUE,rc=NULL) {
    if (!is(input,"GRanges") && !is(input,"GAlignments") && !is.list(input))
        stop("The input argument must be a GenomicRanges object, a ",
            "GenomicAlignments object or a list of either")
    if (!is(mask,"GRanges") && !is(mask,"GRangesList"))
        stop("The mask argument must be a GRanges or GRangesList object")
    if (!is.null(strand) && !is.list(strand) && !isBam && !isBigWig) {
        message("Retrieving ",strand," reads...")
        input <- input[strand(input)==strand]
    }
    cov <- coverageFromRanges(input,mask,ignore.strand,rc=rc)
    gc(verbose=FALSE)
    return(cov) # Rle
}

coverageFromRanges <- function(input,mask,ignore.strand,rc=NULL) {
    allChrs <- as.character(seqlevels(input))
    inChrs <- as.character(unique(seqnames(input)))
    maskChrs <- as.character(unique(seqnames(mask)))
    preCov <- coverage(input)
    preCov <- preCov[allChrs]
    chrs <- Reduce("intersect",list(allChrs,inChrs,maskChrs))
    maskList <- split(mask,seqnames(mask))
    maskList <- maskList[chrs]
    covs <- cmclapply(names(maskList),function(x,maskList,preCov) {
        return(lazyRangesCoverage(x,maskList,preCov))
    },maskList,preCov,rc=rc)
    covs <- unlist(covs)
    if (is.null(names(covs)))
        names(covs) <- names(mask)
    return(covs)
}

lazyRangesCoverage <- function(x,maskList,preCov) {
    message("    processing ",x)
    m <- maskList[[x]]
    pre <- preCov[[x]]
    if (!is.null(m) && !is.null(pre)) { # Sanity...
        co <- pre
        V <- Views(co,ranges(m))
        cot <- unlist(viewApply(V,function(x) x))
        names(cot) <- names(m)
        inv <- which(strand(m)=="-")
        if (length(inv))
            cot[inv] <- lapply(cot[inv],rev)
        return(cot)
    }
    else {
        message("    ",x," not found!")
        return(Rle(NA))
    }
}

prepareBam <- function(b) {
    tryCatch({
        # Will fail if BAM unsorted
        message("Indexing BAM file ",b)
        indexBam(b)
    },error=function(e) {
        warning("Caught error ",e," while indexing BAM file ",b,
            "! Will try to sort now...",immediate.=TRUE)
        message("Sorting BAM file ",b)
        file.rename(b,paste0(b,".uns"))
        ff <- sub(pattern="(.*)\\..*$",replacement="\\1",b)
        sortBam(paste0(b,".uns"),ff)
        file.remove(paste0(b,".uns"))
        message("Indexing BAM file ",b)
        indexBam(b)
    },finally="")
}

.averageAcrossSamples <- function(results,what) {
    if (what != "totals") {
        minCov <- as.matrix(do.call("cbind",lapply(results,function(x) {
            return(x[[what]][,"min_coverage"])
        })))
        
        maxCov <- as.matrix(do.call("cbind",lapply(results,function(x) {
            return(x[[what]][,"max_coverage"])
        })))
        
        meanCov <- as.matrix(do.call("cbind",lapply(results,function(x) {
            return(x[[what]][,"mean_coverage"])
        })))
        
        medianCov <- as.matrix(do.call("cbind",lapply(results,function(x) {
            return(x[[what]][,"median_coverage"])
        })))
        
        fbCov <- as.matrix(do.call("cbind",lapply(results,function(x) {
            return(x[[what]][,"first_base_coverage"])
        })))
        
        lbCov <- as.matrix(do.call("cbind",lapply(results,function(x) {
            return(x[[what]][,"last_base_coverage"])
        })))
    }
    else {
        minCov <- as.matrix(do.call("cbind",lapply(results,function(x) {
            return(x[[what]][,"average_min_coverage"])
        })))
        
        maxCov <- as.matrix(do.call("cbind",lapply(results,function(x) {
            return(x[[what]][,"average_max_coverage"])
        })))
        
        meanCov <- as.matrix(do.call("cbind",lapply(results,function(x) {
            return(x[[what]][,"average_mean_coverage"])
        })))
        
        medianCov <- as.matrix(do.call("cbind",lapply(results,function(x) {
            return(x[[what]][,"average_median_coverage"])
        })))
        
        fbCov <- as.matrix(do.call("cbind",lapply(results,function(x) {
            return(x[[what]][,"average_first_base_coverage"])
        })))
        
        lbCov <- as.matrix(do.call("cbind",lapply(results,function(x) {
            return(x[[what]][,"average_last_base_coverage"])
        })))
    }
    
    
    return(data.frame(
        average_min_coverage=apply(minCov,1,mean),
        average_max_coverage=apply(maxCov,1,mean),
        average_mean_coverage=apply(meanCov,1,mean),
        average_median_coverage=apply(medianCov,1,mean),
        average_first_base_coverage=apply(fbCov,1,mean),
        average_last_base_coverage=apply(lbCov,1,mean)
    ))
}

cmclapply <- function(...,rc) {
    if (suppressWarnings(!requireNamespace("parallel")) 
        || .Platform$OS.type!="unix")
        m <- FALSE
    else {
        m <- TRUE
        ncores <- parallel::detectCores()
        if (ncores==1) 
            m <- FALSE
        else {
            if (!missing(rc) && !is.null(rc))
                ncores <- ceiling(rc*ncores)
            else 
                m <- FALSE
        }
    }
    if (m)
        return(mclapply(...,mc.cores=ncores,mc.set.seed=FALSE))
    else
        return(lapply(...))
}


#bams <- c("IonXpress_027_R_2019_04_05_11_42_09_user_IONAS-415-GK_PHlab_PMN_190405_GK3R227-240_PH3R59-64_PMN1-6_Auto_user_IONAS-415-GK_PHlab_PMN_190405_GK3R227-240_PH3R59-64_PMN1-6_581.bam")
#targets <- "../../temp/CancerHotSpot-v2.dna_manifest.20180509_chr.bed"

# Test
testIt <- function() {
    library(GenomicAlignments)
    library(GenomicRanges)
    library(Rsamtools)
    library(openxlsx)
    
    bams <- dir("/media/raid/data/fleming/precmed/emqn2/smallbam",
        pattern=".bam$")
    targets <- "/media/raid/data/fleming/precmed/pmn1/solid_custom_merged.bed"
    controls <- "/media/raid/data/fleming/precmed/pmn1/control_amplicons.txt"
    fragLen <- 300
    asObject <- FALSE
    pairedOpts <- list(isPaired=NULL,asPaired=FALSE,ifPaired="list")
    outFormat <- c("xlsx","txt")
    avgSamples <- TRUE
    restricBamToTargets <- FALSE
    rc <- 0.75
    
    
}




