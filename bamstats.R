getBamStats <- function(targets,mq=20,ref=NULL,splicing=NULL,
    outFormat=c("txt","xlsx"),outDir=getwd(),outBase=NULL,rc=NULL,
    .verbose=ifelse(is.null(rc),TRUE,FALSE)) {
    if (!requireNamespace("Rsamtools"))
        stop("Bioconductor package Rsamtools is required!")
    if (!requireNamespace("GenomicAlignments"))
        stop("Bioconductor package GenomicAlignments is required!")
        
    if ("xlsx" %in% outFormat && !requireNamespace("openxlsx")) {
        warning("The required CRAN package openxlsx is not installed! Only ",
            "text output will be produced...",immediate.=TRUE)
        outFormat <- "txt"
    }
    
    if (!is.list(targets)) {
        if (file.exists(targets))
            targets <- readTargets(targets)
        else
            stopwrap("targets must be either the result of readTargets ",
                "function or a valid targets file!")
    }
    
    # The BAM files
    bams <- unlist(targets$files,use.names=FALSE)
    nams <- unlist(targets$samples,use.names=FALSE)
    
    # Single or paired?
    if (!is.null(targets$paired))
        paired <- ifelse(targets$paired[[1]][1]=="single",FALSE,TRUE)
    else { # Pairing not provided with targets, must guess - BAM index required
        bi <- vapply(bams,function(x) {
            return(file.exists(paste0(x,".bai")))
        },logical(1))
        if (!all(bi)) {
            nbi <- which(!bi)
            cmclapply(bams[nbi],.prepareBam,rc=rc)
        }
        paired <- testPairedEndBam(bams[1])
        message("")
    }
    
    # Splicing?
    if (!is.logical(splicing)) {
        message("Checking for spliced reads")
        splicing <- .detectSplicing(bams[1])
    }
    
    refRegions <- NULL
    if (!is.null(ref)) {
        if (is.list(ref)) { # metaseqR2 takes care of reference regions
            if (is.null(ref$org))
                stop("A metaseqR2 supported or imported organism name must be ",
                    " provided when ref is a list!")
            if (is.null(ref$refdb)) {
                warning("A metaseqR2 supported or imported reference source ",
                    "must be provided when ref is a list! Using ensembl...",
                    immediate.=TRUE)
                ref$refdb <- "ensembl"
            }
            if (is.null(ref$version)) {
                warning("A version should be provided with a metaseqR2 ",
                    "supported or imported organism! Will use the latest...",
                    immediate.=TRUE)
                ref$version <- "auto"
            }
            if (is.null(ref$localDb))
                # No warning as loadAnnotation will try to fetch on-the-fly
                # but let's try default location
                ref$localDb <- file.path(system.file(package="metaseqR2"),
                    "annotation.sqlite")
            
            if (is.null(ref$countType))
                stop("The countType must be provided when ref is a list!")
            else {
                if (ref$countType == "utr" && is.null(ref$opts)) {
                    warning("3' UTR reference region options not provided, ",
                        "will use the defaults...")
                    ref$opts <- list(frac=1,minLength=300,downstream=50)
                }
                if (ref$countType == "exon" && is.null(ref$opts)) {
                    warning("Exon reference region options not provided, ",
                        "will use the defaults...")
                    ref$opts <- 50
                }
            }
            
            refRegions  <- loadAnnotation(genome=ref$org,refdb=ref$refdb,
                level="gene",type=ref$countType,version=ref$version,
                db=ref$localDb,summarized=TRUE)
            
            # Now, redefine reference regions according to options
            switch(ref$countType,
                gene = {},
                exon = {
                    refRegions <- resize(refRegions,
                        width=width(refRegions) + ref$opts,fix="start")
                    refRegions <- resize(refRegions,
                        width=width(refRegions) + ref$opts,fix="end")
                },
                utr = {
                    nw <- ref$opts$frac*width(refRegions)
                    nw <- ifelse(nw < ref$opts$minLength,ref$opts$minLength,nw)
                    refRegions <- resize(refRegions,width=nw,fix="end")
                    refRegions <- resize(refRegions,
                        width=width(refRegions) + ref$opts$downstream,
                        fix="start")
                }
            )
        }
        
        # Is it a BED file? No additional options then
        if (is.character(ref) && file.exists(ref)) {
            tmp <- read.delim(ref,header=FALSE)
            # FIXME: Name according to bed columns
            names(tmp) <- c("chromosome","start","end","name","score","strand",
                "x","y")
            refRegions <- tryCatch(GRanges(tmp),error=function(e) {
                message("Caught error ",e)
                stop(ref," does not seem to be a valid BED file!")
            },finally="")
        }
    }

    # Now, the main job
    if (paired)
        statsList <- cmclapply(bams,getPairedBamStats,refRegions,mq,
            splicing,.verbose,rc=rc)
    else
        statsList <- cmclapply(bams,getSingleBamStats,refRegions,mq,
            splicing,.verbose,rc=rc)
    
    # Collect results
    readsList <- lapply(statsList,function(x) {
        return(unlist(x$reads))
    })
    pctList <- lapply(statsList,function(x) {
        return(unlist(x$pct))
    })
    hybList <- lapply(statsList,function(x) {
        return(unlist(x$hyb))
    })
    names(readsList) <- names(pctList) <- names(hybList) <- nams
    
    readsDf <- data.frame(sample_name=nams,do.call("rbind",readsList))
    pctDf <- data.frame(sample_name=nams,do.call("rbind",pctList))
    hybDf <- data.frame(sample_name=nams,do.call("rbind",hybList))
    
    if (is.null(outBase))
        outBase <- paste0("alignment_statistics_",
            format(Sys.time(),"%Y-%m-%d-%H-%M-%S"))
    if ("txt" %in% outFormat) {
        readsOutput <- file.path(outDir,paste0(outBase,"_READS.txt"))
        pctOutput <- file.path(outDir,paste0(outBase,"_PCT.txt"))
        hybOutput <- file.path(outDir,paste0(outBase,"_HYBRID.txt"))
        write.table(readsDf,file=readsOutput,sep="\t",quote=FALSE,
            row.names=FALSE)
        write.table(pctDf,file=pctOutput,sep="\t",quote=FALSE,row.names=FALSE)
        write.table(hybDf,file=hybOutput,sep="\t",quote=FALSE,row.names=FALSE)
    }
    if ("xlsx" %in% outFormat) {
        xlsList <- list(reads=readsDf,percentages=pctDf,hybrid=hybDf)
        output <- file.path(outDir,paste0(outBase,".xlsx"))
        write.xlsx(xlsList,file=output,keepNA=TRUE)
    }
    
    return(list(reads=readsDf,pct=pctDf,hyb=hybDf))
}

getSingleBamStats <- function(bam,targets=NULL,mq=20,splicing=FALSE,
    .verbose=TRUE) {
    message("Getting alignment statistics for file ",bam)
    
    # Total reads and bases
    if (.verbose) message("  total reads and bases")
    stru <- .getTotalReads(bam)
    totalReads <- stru$reads
    totalBases <- stru$bases
    totalReadsPct <- "100%"
    totalBasesPct <- "100%"
    totalReadsHyb <- paste0(totalReads," (100%)")
    totalBasesHyb <- paste0(totalBases," (100%)")

    # Aligned reads and bases
    if (.verbose) message("  aligned reads and bases")
    stru <- .getAlignedReads(bam)
    alignedReads <- stru$reads
    alignedBases <- stru$bases
    alignedReadsPct <- paste0(round(100*alignedReads/totalReads,2),"%")
    alignedBasesPct <- paste0(round(100*alignedBases/totalBases,2),"%")
    alignedReadsHyb <- paste0(alignedReads," (",alignedReadsPct,")")
    alignedBasesHyb <- paste0(alignedBases," (",alignedBasesPct,")")

    # Uniquely aligned reads and bases
    if (.verbose) message("  uniquely aligned reads and bases")
    stru <- .getAlignedReads(bam,mq=mq)
    uniqAlignedReads <- stru$reads
    uniqAlignedBases <- stru$bases
    uniqAlignedReadsPct <- paste0(round(100*uniqAlignedReads/totalReads,2),"%")
    uniqAlignedBasesPct <- paste0(round(100*uniqAlignedBases/totalBases,2),"%")
    uniqAlignedReadsHyb <- paste0(uniqAlignedReads," (",uniqAlignedReadsPct,")")
    uniqAlignedBasesHyb <- paste0(uniqAlignedBases," (",uniqAlignedBasesPct,")")
    
    if (splicing) {
        # Aligned and spliced reads
        if (.verbose) message("  aligned and spliced reads")
        alignedSplicedReads <- .getAlignedSplicedReads(bam)
        alignedSplicedReadsPct <- 
            paste0(round(100*alignedSplicedReads/totalReads,2),"%")
        alignedSplicedReadsHyb <- 
            paste0(alignedSplicedReads," (",alignedSplicedReadsPct,")")
        
        # Uniquely aligned and spliced reads
        if (.verbose) message("  uniquely aligned and spliced reads")
        uniqAlignedSplicedReads <- .getAlignedSplicedReads(bam,mq=mq)
        uniqAlignedSplicedReadsPct <- 
            paste0(round(100*uniqAlignedSplicedReads/totalReads,2),"%")
        uniqAlignedSplicedReadsHyb <- 
            paste0(uniqAlignedSplicedReads," (",uniqAlignedSplicedReadsPct,")")
    }
    else
        alignedSplicedReads <- alignedSplicedReadsPct <- 
            alignedSplicedReadsHyb <- uniqAlignedSplicedReads <-
            uniqAlignedSplicedReadsPct <- uniqAlignedSplicedReadsHyb <- NULL
    
    # Reads on target regions
    onTargetReads <- onTargetBases <- onTargetUniqReads <- 
        onTargetUniqBases <- 0
    onTargetReadsPct <- onTargetBasesPct <- onTargetUniqReadsPct <- 
        onTargetUniqBasesPct <- "0%"
    onTargetReadsHyb <- onTargetBasesHyb <- onTargetUniqReadsHyb <- 
        onTargetUniqBasesHyb <- "0 (0%)"
    if (!is.null(targets)) {
        targets <- reduce(targets)
        
        if (.verbose) message("  on target aligned reads and bases")
        stru <- .getOnTargetReads(bam,targets,splicing=splicing)
        onTargetReads <- stru$reads
        onTargetBases <- stru$bases
        onTargetReadsPct <- paste0(round(100*onTargetReads/totalReads,2),"%")
        onTargetBasesPct <- paste0(round(100*onTargetBases/totalBases,2),"%")
        onTargetReadsHyb <- paste0(onTargetReads," (",onTargetReadsPct,")")
        onTargetBasesHyb <- paste0(onTargetBases," (",onTargetBasesPct,")")
        
        if (.verbose) message("  on target uniquely aligned reads and bases")
        stru <- .getOnTargetReads(bam,targets,splicing=splicing,mq=mq)
        onTargetUniqReads <- stru$reads
        onTargetUniqBases <- stru$bases
        onTargetUniqReadsPct <- 
            paste0(round(100*onTargetUniqReads/totalReads,2),"%")
        onTargetUniqBasesPct <- 
            paste0(round(100*onTargetUniqBases/totalBases,2),"%")
        onTargetUniqReadsHyb <- 
            paste0(onTargetUniqReads," (",onTargetUniqReadsPct,")")
        onTargetUniqBasesHyb <- 
            paste0(onTargetUniqBases," (",onTargetUniqBasesPct,")")
    }
    
    return(list(
        reads=list(
            total_reads=totalReads,
            aligned_reads=alignedReads,
            uniquely_aligned_reads=uniqAlignedReads,
            aligned_and_spliced_reads=alignedSplicedReads,
            uniquely_aligned_and_spliced_reads=uniqAlignedSplicedReads,
            on_target_reads=onTargetReads,
            on_target_uniquely_aligned_reads=onTargetUniqReads,
            total_bases=totalBases,
            aligned_bases=alignedBases,
            uniquely_aligned_bases=uniqAlignedBases,
            on_target_bases=onTargetBases,
            on_target_uniquely_aligned_bases=onTargetUniqBases
        ),
        pct=list(
            total_reads=totalReadsPct,
            aligned_reads=alignedReadsPct,
            uniquely_aligned_reads=uniqAlignedReadsPct,
            aligned_and_spliced_reads=alignedSplicedReadsPct,
            uniquely_aligned_and_spliced_reads=uniqAlignedSplicedReadsPct,
            on_target_reads=onTargetReadsPct,
            on_target_uniquely_aligned_reads=onTargetUniqReadsPct,
            total_bases=totalBasesPct,
            aligned_bases=alignedBasesPct,
            uniquely_aligned_bases=uniqAlignedBasesPct,
            on_target_bases=onTargetBasesPct,
            on_target_uniquely_aligned_bases=onTargetUniqBasesPct
        ),
        hyb=list(
            total_reads=totalReadsHyb,
            aligned_reads=alignedReadsHyb,
            uniquely_aligned_reads=uniqAlignedReadsHyb,
            aligned_and_spliced_reads=alignedSplicedReadsHyb,
            uniquely_aligned_and_spliced_reads=uniqAlignedSplicedReadsHyb,
            on_target_reads=onTargetReadsHyb,
            on_target_uniquely_aligned_reads=onTargetUniqReadsHyb,
            total_bases=totalBasesHyb,
            aligned_bases=alignedBasesHyb,
            uniquely_aligned_bases=uniqAlignedBasesHyb,
            on_target_bases=onTargetBasesHyb,
            on_target_uniquely_aligned_bases=onTargetUniqBasesHyb
        )
    ))
}

getPairedBamStats <- function(bam,targets=NULL,mq=20,splicing=FALSE,
    .verbose=TRUE) {
    message("Getting alignment statistics for file ",bam)
    
    # Total reads and bases
    if (.verbose) message("  total reads and bases")
    stru <- .getTotalReads(bam)
    totalReads <- stru$reads
    totalBases <- stru$bases
    totalReadsPct <- "100%"
    totalBasesPct <- "100%"
    totalReadsHyb <- paste0(totalReads," (100%)")
    totalBasesHyb <- paste0(totalBases," (100%)")
    
    # Total paired reads
    if (.verbose) message("  total paired reads and bases")
    stru <- .getPairedTotalPairedReads(bam)
    totalPairedReads <- stru$reads
    totalPairedBases <- stru$bases
    totalPairedReadsPct <- paste0(round(100*totalPairedReads/totalReads,2),"%")
    totalPairedBasesPct <- paste0(round(100*totalPairedBases/totalBases,2),"%")
    totalPairedReadsHyb <- paste0(totalPairedReads," (",totalPairedReadsPct,")")
    totalPairedBasesHyb <- paste0(totalPairedBases," (",totalPairedBasesPct,")")
    
    # Total read-pairs
    if (.verbose) message("  total read pairs")
    totalReadPairs <- totalPairedReads/2
    totalReadPairsPct <- "100%"
    totalReadPairsHyb <- paste0(totalReadPairs," (100%)")

    # Aligned reads and bases
    if (.verbose) message("  aligned reads and bases")
    stru <- .getAlignedReads(bam)
    alignedReads <- stru$reads
    alignedBases <- stru$bases
    alignedReadsPct <- paste0(round(100*alignedReads/totalReads,2),"%")
    alignedBasesPct <- paste0(round(100*alignedBases/totalBases,2),"%")
    alignedReadsHyb <- paste0(alignedReads," (",alignedReadsPct,")")
    alignedBasesHyb <- paste0(alignedBases," (",alignedBasesPct,")")
    
    # Properly paired aligned reads
    if (.verbose) message("  properly paired aligned reads and bases")
    stru <- .getPairedProperPairedAlignedReads(bam)
    properPairAlignedReads <- stru$reads
    properPairAlignedBases <- stru$bases
    properPairAlignedReadsPct <- 
        paste0(round(100*properPairAlignedReads/totalReads,2),"%")
    properPairAlignedBasesPct <- 
        paste0(round(100*properPairAlignedBases/totalBases,2),"%")
    properPairAlignedReadsHyb <- 
        paste0(properPairAlignedReads," (",properPairAlignedReadsPct,")")
    properPairAlignedBasesHyb <- 
        paste0(properPairAlignedBases," (",properPairAlignedBasesPct,")")
        
    # Properly paired aligned read-pairs
    if (.verbose) message("  properly paired aligned read pairs")
    properPairAlignedReadPairs <- properPairAlignedReads/2
    properPairAlignedReadPairsPct <- 
        paste0(round(100*properPairAlignedReadPairs/totalReadPairs,2),"%")
    properPairAlignedReadPairsHyb <- paste0(properPairAlignedReadPairs,
        " (",properPairAlignedReadPairsPct,")")

    # Uniquely aligned reads and bases
    if (.verbose) message("  uniquely aligned reads and bases")
    stru <- .getAlignedReads(bam,mq=mq)
    uniqAlignedReads <- stru$reads
    uniqAlignedBases <- stru$bases
    uniqAlignedReadsPct <- paste0(round(100*uniqAlignedReads/totalReads,2),"%")
    uniqAlignedBasesPct <- paste0(round(100*uniqAlignedBases/totalBases,2),"%")
    uniqAlignedReadsHyb <- paste0(uniqAlignedReads," (",uniqAlignedReadsPct,")")
    uniqAlignedBasesHyb <- paste0(uniqAlignedBases," (",uniqAlignedBasesPct,")")
    
    # Properly paired and uniquely aligned reads
    if (.verbose) message("  properly paired uniquely aligned reads and bases")
    stru <- .getPairedProperPairedAlignedReads(bam,mq=mq)
    uniqProperPairAlignedReads <- stru$reads
    uniqProperPairAlignedBases <- stru$bases
    uniqProperPairAlignedReadsPct <- 
        paste0(round(100*uniqProperPairAlignedReads/totalReads,2),"%")
    uniqProperPairAlignedBasesPct <- 
        paste0(round(100*uniqProperPairAlignedBases/totalBases,2),"%")
    uniqProperPairAlignedReadsHyb <- paste0(uniqProperPairAlignedReads,
        " (",uniqProperPairAlignedReadsPct,")")
    uniqProperPairAlignedBasesHyb <- paste0(uniqProperPairAlignedBases,
        " (",uniqProperPairAlignedBasesPct,")")
    
    # Properly paired and uniquely aligned read-pairs
    if (.verbose) message("  properly paired uniquely aligned read pairs")
    uniqProperPairAlignedReadPairs <- uniqProperPairAlignedReads/2
    uniqProperPairAlignedReadPairsPct <- 
        paste0(round(100*uniqProperPairAlignedReadPairs/totalReadPairs,2),"%")
    uniqProperPairAlignedReadPairsHyb <- paste0(uniqProperPairAlignedReadPairs,
        " (",uniqProperPairAlignedReadPairsPct,")")
    
    if (splicing) {
        # Aligned and spliced reads
        if (.verbose) message("  aligned and spliced reads")
        alignedSplicedReads <- .getAlignedSplicedReads(bam)
        alignedSplicedReadsPct <- 
            paste0(round(100*alignedSplicedReads/totalReads,2),"%")
        alignedSplicedReadsHyb <- 
            paste0(alignedSplicedReads," (",alignedSplicedReadsPct,")")
        
        # Uniquely aligned and spliced reads
        if (.verbose) message("  uniquely aligned and spliced reads")
        uniqAlignedSplicedReads <- .getAlignedSplicedReads(bam,mq=mq)
        uniqAlignedSplicedReadsPct <- 
            paste0(round(100*uniqAlignedSplicedReads/totalReads,2),"%")
        uniqAlignedSplicedReadsHyb <- 
            paste0(uniqAlignedSplicedReads," (",uniqAlignedSplicedReadsPct,")")
    }
    else
        alignedSplicedReads <- alignedSplicedReadsPct <- 
            alignedSplicedReadsHyb <- uniqAlignedSplicedReads <-
            uniqAlignedSplicedReadsPct <- uniqAlignedSplicedReadsHyb <- NULL
    
    # Chimeric reads
    if (.verbose) message("  chimeric reads")
    chimericReads <- .getChimericReads(bam)
    chimericReadsPct <- paste0(round(100*chimericReads/totalReads,2),"%")
    chimericReadsHyb <- paste0(round(100*chimericReads/totalReads,2),"%")
    
    # Uniquely aligned chimeric reads
    if (.verbose) message("  uniquely aligned chimeric reads")
    uniqChimericReads <- .getChimericReads(bam,mq=mq)
    uniqChimericReadsPct <- 
        paste0(round(100*uniqChimericReads/totalReads,2),"%")
    uniqChimericReadsHyb <- 
        paste0(uniqChimericReads," (",uniqChimericReadsPct,")")

    # Reads on target regions
    onTargetReads <- onTargetBases <- onTargetUniqReads <- 
        onTargetUniqBases <- 0
    onTargetReadsPct <- onTargetBasesPct <- onTargetUniqReadsPct <- 
        onTargetUniqBasesPct <- "0%"
    onTargetReadsHyb <- onTargetBasesHyb <- onTargetUniqReadsHyb <- 
        onTargetUniqBasesHyb <- "0 (0%)"
    if (!is.null(targets)) {
        targets <- reduce(targets)
        
        if (.verbose) message("  on target aligned reads and bases")
        stru <- .getOnTargetReads(bam,targets,splicing=splicing)
        onTargetReads <- stru$reads
        onTargetBases <- stru$bases
        onTargetReadsPct <- paste0(round(100*onTargetReads/totalReads,2),"%")
        onTargetBasesPct <- paste0(round(100*onTargetBases/totalBases,2),"%")
        onTargetReadsHyb <- paste0(onTargetReads," (",onTargetReadsPct,")")
        onTargetBasesHyb <- paste0(onTargetBases," (",onTargetBasesPct,")")
        
        if (.verbose) message("  on target uniquely aligned reads and bases")
        stru <- .getOnTargetReads(bam,targets,splicing=splicing,mq=mq)
        onTargetUniqReads <- stru$reads
        onTargetUniqBases <- stru$bases
        onTargetUniqReadsPct <- 
            paste0(round(100*onTargetUniqReads/totalReads,2),"%")
        onTargetUniqBasesPct <- 
            paste0(round(100*onTargetUniqBases/totalBases,2),"%")
        onTargetUniqReadsHyb <- 
            paste0(onTargetUniqReads," (",onTargetUniqReadsPct,")")
        onTargetUniqBasesHyb <- 
            paste0(onTargetUniqBases," (",onTargetUniqBasesPct,")")
    }
    
    return(list(
        reads=list(
            total_reads=totalReads,
            total_paired_reads=totalPairedReads,
            total_read_pairs=totalReadPairs,
            aligned_reads=alignedReads,
            proper_paired_aligned_reads=properPairAlignedReads,
            proper_paired_aligned_read_pairs=properPairAlignedReadPairs,
            uniquely_aligned_reads=uniqAlignedReads,
            proper_paired_uniquely_aligned_reads=uniqProperPairAlignedReads,
            proper_paired_uniquely_aligned_read_pairs=
                uniqProperPairAlignedReadPairs,
            aligned_and_spliced_reads=alignedSplicedReads,
            uniquely_aligned_and_spliced_reads=uniqAlignedSplicedReads,
            chimeric_reads=chimericReads,
            uniquely_aligned_chimeric_reads=uniqChimericReads,
            on_target_reads=onTargetReads,
            on_target_uniquely_aligned_reads=onTargetUniqReads,
            total_bases=totalBases,
            total_paired_bases=totalPairedBases,
            aligned_bases=alignedBases,
            proper_paired_aligned_bases=properPairAlignedBases,
            uniquely_aligned_bases=uniqAlignedBases,
            proper_paired_uniquely_aligned_bases=uniqProperPairAlignedBases,
            on_target_bases=onTargetBases,
            on_target_uniquely_aligned_bases=onTargetUniqBases
        ),
        pct=list(
            total_reads=totalReadsPct,
            total_paired_reads=totalPairedReadsPct,
            total_read_pairs=totalReadPairsPct,
            aligned_reads=alignedReadsPct,
            proper_paired_aligned_reads=properPairAlignedReadsPct,
            proper_paired_aligned_read_pairs=properPairAlignedReadPairsPct,
            uniquely_aligned_reads=uniqAlignedReadsPct,
            proper_paired_uniquely_aligned_reads=uniqProperPairAlignedReadsPct,
            proper_paired_uniquely_aligned_read_pairs=
                uniqProperPairAlignedReadPairsPct,
            aligned_and_spliced_reads=alignedSplicedReadsPct,
            uniquely_aligned_and_spliced_reads=uniqAlignedSplicedReadsPct,
            chimeric_reads=chimericReadsPct,
            uniquely_aligned_chimeric_reads=uniqChimericReadsPct,
            on_target_reads=onTargetReadsPct,
            on_target_uniquely_aligned_reads=onTargetUniqReadsPct,
            total_bases=totalBasesPct,
            total_paired_bases=totalPairedBasesPct,
            aligned_bases=alignedBasesPct,
            proper_paired_aligned_bases=properPairAlignedBasesPct,
            uniquely_aligned_bases=uniqAlignedBasesPct,
            proper_paired_uniquely_aligned_bases=uniqProperPairAlignedBasesPct,
            on_target_bases=onTargetBasesPct,
            on_target_uniquely_aligned_bases=onTargetUniqBasesPct
        ),
        hyb=list(
            total_reads=totalReadsHyb,
            total_paired_reads=totalPairedReadsHyb,
            total_read_pairs=totalReadPairsHyb,
            aligned_reads=alignedReadsHyb,
            proper_paired_aligned_reads=properPairAlignedReadsHyb,
            proper_paired_aligned_read_pairs=properPairAlignedReadPairsHyb,
            uniquely_aligned_reads=uniqAlignedReadsHyb,
            proper_paired_uniquely_aligned_reads=uniqProperPairAlignedReadsHyb,
            proper_paired_uniquely_aligned_read_pairs=
                uniqProperPairAlignedReadPairsHyb,
            aligned_and_spliced_reads=alignedSplicedReadsHyb,
            uniquely_aligned_and_spliced_reads=uniqAlignedSplicedReadsHyb,
            chimeric_reads=chimericReadsHyb,
            uniquely_aligned_chimeric_reads=uniqChimericReadsHyb,
            on_target_reads=onTargetReadsHyb,
            on_target_uniquely_aligned_reads=onTargetUniqReadsHyb,
            total_bases=totalBasesHyb,
            total_paired_bases=totalPairedBasesHyb,
            aligned_bases=alignedBasesHyb,
            proper_paired_aligned_bases=properPairAlignedBasesHyb,
            uniquely_aligned_bases=uniqAlignedBasesHyb,
            proper_paired_uniquely_aligned_bases=uniqProperPairAlignedBasesHyb,
            on_target_bases=onTargetBasesHyb,
            on_target_uniquely_aligned_bases=onTargetUniqBasesHyb
        )
    ))
}

.getTotalReads <- function(bam) {
    if (is.character(bam))
        return(.getTotalReadsFile(bam))
    else if (is.list(bam))
        return(.getTotalReadsObj(bam))
}

.getTotalReadsFile <- function(bam) {
    params <- ScanBamParam(
        flag=scanBamFlag(
            isSecondaryAlignment=FALSE,
            isSupplementaryAlignment=FALSE
        )
    )
    countObj <- countBam(bam,param=params)
    nReads <- countObj$records
    nBases <- countObj$nucleotides

    return(list(
        reads=nReads,
        bases=nBases
    ))
}

.getTotalReadsObj <- function(obj) {
    return(list(
        reads=length(obj[[1]]$pos),
        bases=sum(width(obj[[1]]$seq))
    ))
}

.getAlignedReads <- function(bam,mq=NA_integer_) {
    if (is.character(bam))
        return(.getAlignedReadsFile(bam,mq))
    else if (is.list(bam))
        return(.getAlignedReadsObj(bam,mq))
}

.getAlignedReadsFile <- function(bam,mq=NA_integer_) {
    params <- ScanBamParam(
        flag=scanBamFlag(
            isUnmappedQuery=FALSE,
            isSecondaryAlignment=FALSE,
            isSupplementaryAlignment=FALSE
        ),
        mapqFilter=mq
    )
    countObj <- countBam(bam,param=params)
    nReads <- countObj$records
    nBases <- countObj$nucleotides
    
    return(list(
        reads=nReads,
        bases=nBases
    ))
}

.getAlignedReadsObj <- function(obj,mq=NA_integer_) {
    if (is.na(mq))
        ind <- which(!is.na(obj[[1]]$pos))
    else    
        ind <- which(!is.na(obj[[1]]$pos) & obj[[1]]$mapq>mq)
        
    return(list(
        reads=length(ind),
        bases=sum(width(obj[[1]]$seq[ind]))
    ))
}

.getAlignedSplicedReads <- function(bam,mq=NA_integer_) {
    if (is.character(bam))
        return(.getAlignedSplicedReadsFile(bam,mq))
    else if (is.list(bam))
        return(.getAlignedSplicedReadsObj(bam,mq))
}

.getAlignedSplicedReadsFile <- function(bam,mq=NA_integer_) {
    params <- ScanBamParam(
        flag=scanBamFlag(
            isUnmappedQuery=FALSE,
            isSecondaryAlignment=FALSE,
            isSupplementaryAlignment=FALSE
        ),
        what="cigar",
        mapqFilter=mq
    )
    bamObj <- scanBam(bam,param=params)
    nReads <- length(grep("N",bamObj[[1]]$cigar))

    return(nReads)
}

.getAlignedSplicedReadsObj <- function(obj,mq=NA_integer_) {
    if (is.na(mq))
        ind <- which(!is.na(obj[[1]]$pos))
    else    
        ind <- which(!is.na(obj[[1]]$pos) & obj[[1]]$mapq>mq)
    return(length(grep("N",obj[[1]]$cigar[ind])))
}

.getPairedTotalPairedReads <- function(bam) {
    if (is.character(bam))
        return(.getPairedTotalPairedReadsFile(bam))
    else if (is.list(bam))
        return(.getPairedTotalPairedReadsObj(bam))
}

.getPairedTotalPairedReadsFile <- function(bam) {
    params <- ScanBamParam(
        flag=scanBamFlag(
            isPaired=TRUE,
            isSecondaryAlignment=FALSE,
            isSupplementaryAlignment=FALSE
        )
    )
    countObj <- countBam(bam,param=params)
    nReads <- countObj$records
    nBases <- countObj$nucleotides

    return(list(
        reads=nReads,
        bases=nBases
    ))
}

.getPairedTotalPairedReadsObj <- function(obj) {
    return(list(
        reads=length(obj[[1]]$pos),
        bases=sum(width(obj[[1]]$seq))
    ))
}

.getPairedProperPairedAlignedReads <- function(bam,mq=NA_integer_) {
    if (is.character(bam))
        return(.getPairedProperPairedAlignedReadsFile(bam,mq))
    else if (is.list(bam))
        return(.getPairedProperPairedAlignedReadsObj(bam,mq))
}

.getPairedProperPairedAlignedReadsFile <- function(bam,mq=NA_integer_) {
    params <- ScanBamParam(
        flag=scanBamFlag(
            isProperPair=TRUE,
            isSecondaryAlignment=FALSE,
            isSupplementaryAlignment=FALSE
        ),
        mapqFilter=mq
    )
    countObj <- countBam(bam,param=params)
    nReads <- countObj$records
    nBases <- countObj$nucleotides
    
    return(list(
        reads=nReads,
        bases=nBases
    ))
}

.getPairedProperPairedAlignedReadsObj <- function(obj,mq=NA_integer_) {
    # 83 : read paired, read mapped in proper pair, read reverse strand,
    #      first in pair
    # 99 : read paired, read mapped in proper pair, mate reverse strand,
    #      first in pair
    # 147: read paired, read mapped in proper pair, read reverse strand,
    #      second in pair
    # 163: read paired, read mapped in proper pair, mate reverse strand,
    #      second in pair
    if (is.na(mq))
        ind <- which(obj[[1]]$flag %in% c(83,99,147,163))
    else
        ind <- which(obj[[1]]$flag %in% c(83,99,147,163) & obj[[1]]$mapq>mq)
    
    return(list(
        reads=length(ind),
        bases=sum(width(obj[[1]]$seq[ind]))
    ))
}

.getPairedProperPairedAlignedReadPairs <- function(bam) {
    if (is.character(bam))
        return(.getPairedProperPairedAlignedReadPairsFile(bam))
    else if (is.list(bam))
        return(.getPairedProperPairedAlignedReadPairsObj(bam))
}

.getPairedProperPairedAlignedReadPairsFile <- function(bam) {
    params <- ScanBamParam(
        flag=scanBamFlag(
            isProperPair=TRUE,
            isFirstMateRead=TRUE,
            isSecondaryAlignment=FALSE,
            isSupplementaryAlignment=FALSE
        )
    )
    countObj <- countBam(bam,param=params)
    
    return(countObj$records)
}

.getPairedProperPairedAlignedReadPairsObj <- function(obj) {
    return(length(which(objp[[1]]$flag %in% c(83,99,147,163)))/2)
}

.getChimericReads <- function(bam,mq=NA_integer_) {
     if (is.character(bam))
        return(.getChimericReadsFile(bam,mq))
    else if (is.list(bam))
        return(.getChimericReadsObj(bam,mq))
}

.getChimericReadsFile <- function(bam,mq=NA_integer_) {
    params <- ScanBamParam(
        flag=scanBamFlag(
            isProperPair=FALSE,
            isUnmappedQuery=FALSE,
            hasUnmappedMate=FALSE,
            isSecondaryAlignment=FALSE,
            isSupplementaryAlignment=FALSE
        ),
        what=c("rname","mrnm"),
        mapqFilter=mq
    )
    bamObj <- scanBam(bam,param=params)
    
    return(length(which(bamObj[[1]]$rname != bamObj[[1]]$mrnm)))
}

.getChimericReadsObj <- function(obj,mq=NA_integer_) {
    params <- ScanBamParam(
        flag=scanBamFlag(
            isProperPair=FALSE,
            isUnmappedQuery=FALSE,
            hasUnmappedMate=FALSE,
            isSecondaryAlignment=FALSE,
            isSupplementaryAlignment=FALSE
        ),
        what=c("rname","mrnm"),
        mapqFilter=mq
    )
    bamObj <- scanBam(bam,param=params)
    
    # Exclude any properly paired mapped read by flags
    if (is.na(mq))
        ind <- which(!(obj[[1]]$flag %in% c(83,99,147,163)))
    else
        ind <- which(!(obj[[1]]$flag %in% c(83,99,147,163)) & obj[[1]]$mapq>mq)
    
    return(length(which(!(obj[[1]]$flag %in% c(83,99,147,163)) 
        & obj[[1]]$rname != obj[[1]]$mrnm)))
}

.getOnTargetReads <- function(bam,targets,splicing=FALSE,mq=NA_integer_) {
    params <- ScanBamParam(
        flag=scanBamFlag(
            isUnmappedQuery=FALSE,
            isSecondaryAlignment=FALSE,
            isSupplementaryAlignment=FALSE
        ),
        #what=c("seq","cigar"),
        what="seq",
        mapqFilter=mq
    )
    
    reads <- as(readGAlignments(file=bam,param=params),"GRanges")
    ov <- findOverlaps(reads,targets,ignore.strand=TRUE)
    ind <- queryHits(ov)
    if (splicing)
        ind <- unique(ind)
    seqs <- reads$seq[ind]
    
    nReads <- length(ind)
    nBases <- sum(lengths(seqs))
    
    return(list(
        reads=nReads,
        bases=nBases
    ))
}

.getOnTargetReadsFile <- function(bam,targets,splicing=FALSE,mq=NA_integer_) {
    params <- ScanBamParam(
        flag=scanBamFlag(
            isUnmappedQuery=FALSE,
            isSecondaryAlignment=FALSE,
            isSupplementaryAlignment=FALSE
        ),
        #what=c("seq","cigar"),
        what="seq",
        mapqFilter=mq
    )
    
    reads <- as(readGAlignments(file=bam,param=params),"GRanges")
    ov <- findOverlaps(reads,targets,ignore.strand=TRUE)
    ind <- queryHits(ov)
    if (splicing)
        ind <- unique(ind)
    seqs <- reads$seq[ind]
    
    nReads <- length(ind)
    nBases <- sum(lengths(seqs))
    
    return(list(
        reads=nReads,
        bases=nBases
    ))
}

.getOnTargetReadsObj <- function(obj,targets,splicing=FALSE,mq=NA_integer_) {
    if (is.na(mq))
        ind <- which(!is.na(obj[[1]]$pos))
    else    
        ind <- which(!is.na(obj[[1]]$pos) & obj[[1]]$mapq>mq)
    
    # seqlengths are not really required for this operation
    reads <- as(GAlignments(
        seqnames=obj[[1]]$rname[ind],
        pos=obj[[1]]$pos[ind],
        cigar=obj[[1]]$cigar[ind],
        strand=obj[[1]]$strand[ind],
        seqlengths=
    ),"GRanges")
    
    ov <- findOverlaps(reads,targets,ignore.strand=TRUE)
    ind <- queryHits(ov)
    if (splicing)
        ind <- unique(ind)
    seqs <- reads$seq[ind]
    
    return(list(
        reads=length(ind),
        bases=sum(lengths(seqs))
    ))
}

.detectSplicing <- function(bam) {
    params <- ScanBamParam(
        flag=scanBamFlag(
            isUnmappedQuery=FALSE,
            isSecondaryAlignment=FALSE,
            isSupplementaryAlignment=FALSE
        ),
        what="cigar"
    )
    bamFile <- BamFile(bam)
    yieldSize(bamFile) <- 100000
    bamObj <- scanBam(bamFile,param=params)
    return(length(grep("N",bamObj[[1]]$cigar))>0)
}

.prepareBam <- function(bam) {
    tryCatch({
        # Will fail if BAM unsorted
        message("Indexing BAM file ",bam)
        indexBam(bam)
    },error=function(e) {
        warning("Caught error ",e," while indexing BAM file ",bam,"! ",
            "Will try to sort now...",immediate.=TRUE)
        message("Sorting BAM file ",bam)
        file.rename(bam,paste0(bam,".uns"))
        ff <- sub(pattern="(.*)\\..*$",replacement="\\1",bam)
        sortBam(paste0(bam,".uns"),ff)
        file.remove(paste0(bam,".uns"))
        message("Indexing BAM file ",bam)
        indexBam(bam)
    },finally="")
}

.readBamObj <- function(bam,paired) {
    # Pairing info not required if not paired-end
    if (paired)
        what <- scanBamWhat()[-1]
    else
        what <- scanBamWhat()[-c(1,9,10,11)]
    params <- ScanBamParam(
        flag=scanBamFlag(
            isSecondaryAlignment=FALSE,
            isSupplementaryAlignment=FALSE
        ),
        what=what
    )
    return(scanBam(bam,param=params))
}

cmclapply <- function (..., rc) {
    if (suppressWarnings(!requireNamespace("parallel")) 
        || .Platform$OS.type != "unix")
        m <- FALSE
    else {
        m <- TRUE
        ncores <- parallel::detectCores()
        if (ncores == 1)
            m <- FALSE
        else {
            if (!missing(rc) && !is.null(rc))
                ncores <- ceiling(rc * ncores)
            else m <- FALSE
        }
    }
    if (m)
        return(mclapply(..., mc.cores = ncores, mc.set.seed = FALSE))
    else return(lapply(...))
}

################################################################################
#
#.getOnTargetReadsSlow <- function(bam,targets,mq=NA_integer_) {
#    params <- ScanBamParam(
#        flag=scanBamFlag(
#            isUnmappedQuery=FALSE,
#            isSecondaryAlignment=FALSE,
#            isSupplementaryAlignment=FALSE
#        ),
#        which=targets,
#        mapqFilter=mq
#    )
#    countObj <- countBam(bam,param=params)
#    nReads <- sum(countObj$records)
#    nBases <- sum(countObj$nucleotides)
#    
#    return(list(
#        reads=nReads,
#        bases=nBases
#    ))
#}
#
