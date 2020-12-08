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
# for targetted sequencing.
#
coverageStats <- function(bams,targets=NULL,fragLen=300,asObject=FALSE,
    pairedOpts=list(isPaired=NULL,asPaired=FALSE,ifPaired="list"),
    outFormat=c("xlsx","txt"),avgSamples=FALSE,rc=NULL) {
    if (missing(bams))
        stop("One or more BAM files must be provided!")
    # TODO: Add check if any file does not exist
    
    if (is.null(targets))
        warning("Targets not provided! Total coverage will be calculated",
            immediate.=TRUE)
    # TODO: Add check if targets file does not exist if provided
    
    if (!require(GenomicAlignments))
        stop("Bioconductor package GenomicAlignments is required!")
    if (!require(GenomicRanges))
        stop("Bioconductor package GenomicAlignments is required!")
        
    outFormat <- outFormat[1]
    if (!is.character(outFormat) || !(outFormat %in% c("txt","xlsx")))
        stop("outFormat must be one of \"txt\", \"xlsx\"")
    if (outFormat == "xlsx" && !require(openxlsx))
        stop("R package openxlsx is required for Excel output")
        
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
    
    # Import reads supporting the targets
    bamRanges <- NULL
    if (!is.null(targets)) {
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
    
    if (!is.null(rc) && length(bams) > 1) {
        results <- cmclapply(basename(bams),function(b,a) {
            message("Sample ",b)
            message("  reading BAM file")
            if (is.null(targets))
                reads <- readBamAlignments(b)
            else
                # This is not working very well... It is better to read the
                # total BAM file and calculate overlaps... Also much slower.
                reads <- readBamAlignments(b,gr=bamRanges)
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
            
            # Assemble
            bs <- basename(b)
            D <- data.frame(
                min_coverage=minCo,
                max_coverage=maxCo,
                mean_coverage=meanCo,
                median_coverage=medianCo,
                first_base_coverage=firstCo,
                last_base_coverage=lastCo,
                row.names=names(co)
            )
            ## Add total stats
            #fRow <- c(
            #    mean(min(D$min_coverage)),
            #    mean(max(D$max_coverage)),
            #    mean(D$mean_coverage),
            #    mean(D$median_coverage),
            #    mean(D$first_base_coverage),
            #    mean(D$last_base_coverage)
            #)
            #D <- rbind(D,fRow)
            #rownames(D)[nrow(D)] <- "TOTAL_AVERAGE"
            # Explicit names
            D <- cbind(name=rownames(D),D)
            
            if (!is.null(a))
                D <- cbind(a,D)
            
            return(D)
        },att,rc=rc)
        names(results) <- basename(bams)
    }
    else {
        results <- vector("list",length(bams))
        names(results) <- basename(bams)
        
        for (b in bams) {
            message("Sample ",b)
            message("  reading BAM file")
            if (is.null(targets))
                reads <- readBamAlignments(b)
            else {
                #reads <- readBamAlignments(b,gr=bamRanges)
                reads <- readBamAlignments(b,asRanges=TRUE,seq=FALSE)
                tmp <- reduce(bamRanges)
                ov <- findOverlaps(reads,tmp)
                ind <- unique(queryHits(ov))
                reads <- reads[ind]
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
            
            # Assemble
            bs <- basename(b)
            results[[bs]] <- data.frame(
                min_coverage=minCo,
                max_coverage=maxCo,
                mean_coverage=meanCo,
                median_coverage=medianCo,
                first_base_coverage=firstCo,
                last_base_coverage=lastCo,
                row.names=names(co)
            )
            ## Add total stats
            #fRow <- c(
            #    mean(min(results[[bs]]$min_coverage)),
            #    mean(max(results[[bs]]$max_coverage)),
            #    mean(results[[bs]]$mean_coverage),
            #    mean(results[[bs]]$median_coverage),
            #    mean(results[[bs]]$first_base_coverage),
            #    mean(results[[bs]]$last_base_coverage)
            #)
            #results[[bs]] <- rbind(results[[bs]],fRow)
            #rownames(results[[bs]])[nrow(results[[bs]])] <- "TOTAL_AVERAGE"
            # Explicit names
            results[[bs]] <- cbind(name=rownames(results[[bs]]),results[[bs]])
            
            if (!is.null(att))
                results[[bs]] <- cbind(att,results[[bs]])
        }
    }
       
    #if (avgSamples && length(bams) > 1) {
    #   minco <- sapply(results,function(x) {
    #       x[,"min_coverage"]
    #   })
    #   maxco <- sapply(results,function(x) {
    #       x[,"max_coverage"]
    #   })
    #   meanco <- sapply(results,function(x) {
    #       x[,"mean_coverage"]
    #   })
    #   medianco <- sapply(results,function(x) {
    #       x[,"median_coverage"]
    #   })
    #}
    
    if (outFormat == "txt") {
        message("Writing summary output to separate text files")
        for (b in bams) {
            base <- basename(b)
            parts <- strsplit(base,"\\.")[[1]]
            if (length(parts) == 1) # No extension
                out <- file.path(dirname(b),paste(parts,".covstats.txt",sep=""))
            else
                out <- paste(paste(parts[1:(length(parts)-1)],collapse="."),
                    "covstats",parts[length(parts)],sep=".")
            out <- file.path(dirname(b),out)
            write.table(results[[base]],file=out,sep="\t",quote=FALSE,
                row.names=FALSE)
        }
    }
    else if (outFormat == "xlsx") {
        message("Writing summary output to a single xls file")
        out <- file.path(dirname(bams[1]),paste("covstats_summary_",
            format(Sys.time(),"%Y-%m-%d-%H-%M-%S"),".",outFormat,sep=""))
        write.xlsx(results,file=out,keepNA=TRUE)
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
                params <- ScanBamParam(which=gr,what="seq")
            else
                params <- ScanBamParam(which=gr)
        }
        else {
            if (seq)
                params <- ScanBamParam(what="seq")
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
