#' Peak union calculation
#' 
#' 
#' The function goes over each BAM file in the directory and finds the expression peaks that satisfy the coverage boundary and length criteria in each file. Then it unifies the peak information to obtain a single set of peak genomic coordinates.
#' 
#' @param bam_location The directory containing BAM files.
#' @param bam_txt_list Optional newline separated text file of filenames of bam files. File must be located in bam_location directory. 
#' @param low_coverage_cutoff An integer indicating the low coverage threshold value.
#' @param high_coverage_cutoff An integer indicating the high coverage threshold value.
#' @param peak_width An integer indicating the minimum peak width.
#' @param paired_end_data A boolean indicating if the reads are paired-end.
#' @param strandedness A string outlining the type of the sequencing library: stranded, or reversely stranded. Defaults to "stranded"; "unstranded" is rejected with an error.
#' @param scanbamparam An optional \code{Rsamtools::ScanBamParam} object giving
#'   full control of the BAM read filter. When \code{NULL} (the default) an
#'   internal filter is built from \code{mapqFilter} plus alignment-flag
#'   exclusions: unmapped, QC-failing, secondary and supplementary alignments
#'   are dropped, and for paired-end data the read must additionally be properly
#'   paired with a mapped mate. A supplied \code{scanbamparam} takes precedence
#'   and \code{mapqFilter} is then ignored.
#' @param mapqFilter Integer. Minimum mapping quality (MAPQ) a read must have to
#'   be retained during coverage construction, used only when
#'   \code{scanbamparam = NULL}. Default 10, which keeps uniquely-mapped reads
#'   on the common bacterial aligners while discarding the ambiguous and
#'   multi-mapping reads that cluster at MAPQ 0 to 3. MAPQ is aligner-specific,
#'   so do NOT set \code{mapqFilter} above the MAPQ your aligner assigns a
#'   uniquely-mapped read, or all reads are discarded and coverage is empty.
#'   Maximum recommended values by aligner: BWA-MEM 60, minimap2 60, Rsubread
#'   align/subjunc 40, BWA aln (backtrack) 37, Bowtie2 42 (but see note). Note:
#'   Bowtie2 uses a non-monotonic MAPQ in which some uniquely-mapped reads that
#'   carry mismatches score only 3 or 8, below the default of 10, so Bowtie2
#'   users should set \code{mapqFilter = 1}. Set \code{mapqFilter = NA} to
#'   disable mapping-quality filtering while keeping the flag exclusions. A
#'   read-quality filter that retains no reads triggers a warning.
#' @param coverage_model Character, one of "fragment" or "footprint", setting how
#'   paired-end reads contribute to coverage. It applies only to paired-end data
#'   and is ignored for single-end. "fragment" (the default) treats each read
#'   pair as one span from its leftmost to rightmost aligned base, so the
#'   unsequenced insert between the mates is counted as covered. "footprint"
#'   instead counts only the aligned blocks of the two mates, so the insert
#'   between them is left uncovered and the coverage reflects only the bases
#'   that were actually sequenced.
#' 
#' @return A named list with two IRanges objects, `plus` and `minus`, holding the
#'   unified peak coordinates for each strand.
#' 
#' @import IRanges
#' @import GenomicAlignments
#' @import Rsamtools
#' @importFrom utils capture.output read.csv read.delim write.table
#'
#' @export
peak_union_calc <- function(bam_location = ".", bam_txt_list = "", low_coverage_cutoff, high_coverage_cutoff, peak_width, paired_end_data = FALSE, strandedness = "stranded", scanbamparam = NULL, mapqFilter = 10, coverage_model = c("fragment", "footprint")) {
  coverage_model <- match.arg(coverage_model)
  ## Validate strandedness at entry, before any BAM is read.
  valid_strandedness <- c("stranded", "reversely_stranded")
  if (!strandedness %in% valid_strandedness) {
    stop("Invalid 'strandedness' value: '", strandedness,
         "'. Feature prediction requires stranded data; must be one of: ",
         paste(valid_strandedness, collapse = ", "),
         ". ('unstranded' is valid for read counting but not for prediction.)",
         call. = FALSE)
  }
  ## Find all BAM files in the directory.
  if (bam_txt_list != ""){
    bam_files <- readLines(bam_txt_list)
    bam_files <- lapply(bam_files, function(x) paste(bam_location, x, sep = "/"))
  } else {
    bam_files <- list.files(path = bam_location, pattern = "\\.BAM$", full.names = TRUE, ignore.case = TRUE)
  }

  ## Helper: selected peak IRanges for one strand's reads.
  compute_strand_peaks <- function(reads, f) {
    strand_cvg <- coverage(reads)
    list_components <- names(strand_cvg)
    if (length(list_components) != 1) {
      stop(paste("Invalid BAM file:", f, sep = " "))
    }
    target <- list_components
    ## Cut the coverage vector to obtain the expression peaks above the low cut-off.
    strand_rle <- strand_cvg[[target]]
    peaks <- IRanges::slice(strand_rle, lower = low_coverage_cutoff, includeLower = TRUE)
    ## Count, per peak, the positions above the high cut-off, then keep the peaks
    ## whose count exceeds the minimum width. This reproduces the count-based test.
    above_high_rle <- strand_rle > high_coverage_cutoff
    above_high_views <- Views(above_high_rle, ranges(peaks))
    positions_above <- viewSums(above_high_views)
    selected_peaks <- peaks[positions_above > peak_width]
    IRanges(start = start(selected_peaks), end = end(selected_peaks))
  }

  plus_union  <- IRanges()
  minus_union <- IRanges()
  empty_bams  <- character(0)
  ## Build the default read-quality filter when none is supplied.
  ## mapqFilter (default 10) sets the minimum MAPQ and is used only on this
  ## default path; a supplied scanbamparam takes precedence and mapqFilter is
  ## ignored. Default 10 keeps uniquely-mapped reads on the common bacterial
  ## aligners while dropping the ambiguous and multi-mapping reads that cluster
  ## at MAPQ 0 to 3. MAPQ is aligner-specific: never set it above the MAPQ your
  ## aligner gives a uniquely-mapped read, or coverage empties. Unique-read
  ## ceilings: BWA-MEM 60, minimap2 60, Rsubread 40, bwa aln 37, bowtie2 42.
  ## Bowtie2 is the exception: its non-monotonic MAPQ scores some unique reads
  ## carrying mismatches 3 or 8 (below 10), so Bowtie2 users should pass
  ## mapqFilter = 1. Set mapqFilter = NA to disable MAPQ filtering while keeping
  ## the flag exclusions; for total control pass a full ScanBamParam.
  if (is.null(scanbamparam)) {
    if (paired_end_data) {
      scanbamflag <- scanBamFlag(isUnmappedQuery = FALSE,
                                 isPaired = TRUE,
                                 isProperPair = TRUE,
                                 isNotPassingQualityControls = FALSE,
                                 hasUnmappedMate = FALSE,
                                 isSecondaryAlignment = FALSE,
                                 isSupplementaryAlignment = FALSE)
    } else {
      scanbamflag <- scanBamFlag(isUnmappedQuery = FALSE,
                                 isNotPassingQualityControls = FALSE,
                                 isSecondaryAlignment = FALSE,
                                 isSupplementaryAlignment = FALSE)
    }
    scanbamparam <- ScanBamParam(flag = scanbamflag, mapqFilter = mapqFilter)
  }
  ## Read each BAM once, split by strand, and accumulate both peak unions.
  for (f in bam_files) {
    if (paired_end_data) {
      strand_mode <- if (strandedness == "reversely_stranded") 2 else 1
      read_pairs <- readGAlignmentPairs(f, strandMode = strand_mode, param = scanbamparam)
      ## coverage_model (paired-end only): "fragment" (default) counts the whole
      ## pair from leftmost to rightmost base; "footprint" counts only the aligned
      ## blocks of each mate, leaving the gap between mates uncovered.
      if (coverage_model == "footprint") {
        file_alignment <- unlist(grglist(read_pairs))
      } else {
        file_alignment <- granges(read_pairs)
      }
    } else {
      file_alignment <- readGAlignments(f, param = scanbamparam)
    }
    if (length(file_alignment) == 0L) empty_bams <- c(empty_bams, f)
    reads_plus  <- file_alignment[strand(file_alignment) == "+"]
    reads_minus <- file_alignment[strand(file_alignment) == "-"]
    ## Single-end reversely-stranded data interprets the read strand in reverse.
    if (!paired_end_data & strandedness == "reversely_stranded") {
      plus_reads  <- reads_minus
      minus_reads <- reads_plus
    } else {
      plus_reads  <- reads_plus
      minus_reads <- reads_minus
    }
    plus_union  <- IRanges::union(plus_union,  compute_strand_peaks(plus_reads,  f))
    minus_union <- IRanges::union(minus_union, compute_strand_peaks(minus_reads, f))
  }
  ## Warn once if an active mapqFilter left any BAM with no reads.
  if (length(empty_bams) > 0L && !is.null(scanbamparam) && !is.na(bamMapqFilter(scanbamparam))) {
    warning(paste0(length(empty_bams), " BAM file(s) yielded no reads after filtering with ",
                   "mapqFilter = ", bamMapqFilter(scanbamparam), " (e.g. ", basename(empty_bams[1]),
                   "). This usually means the threshold exceeds the aligner's maximum MAPQ ",
                   "(bwa aln ~37, bowtie2 ~42); lower mapqFilter or set it to NA."),
            call. = FALSE, immediate. = TRUE)
  }
  return(list(plus = plus_union, minus = minus_union))
}

#' Peak checking for the second coverage threshold and width.
#' 
#' This is a helper function that is used to examine if the peak had a continuous stretch of a given width that has coverage above the high cut-off value.
#' 
#' @param View_line A line from a RleViews object.
#' @param high_cutoff An integer indicating the high coverage threshold value.
#' @param min_sRNA_length An integer indicating the minimum sRNA length (peak width).
#' 
#' @return Returns a RleViews line if it satisfies conditions.
#' 
#' @export
#' 
peak_analysis <- function(View_line, high_cutoff, min_sRNA_length) {
  ## This is a helper function that is used to examine if the peak had a continuous stretch of a given width that has coverage above the high cut-off value.
  cvg_string <- as.vector(View_line)
  target_peak <- which(cvg_string>high_cutoff)
  if (length(target_peak)>min_sRNA_length) {
    return(View_line)
  }
}


#' Extract major features from the annotation file
#' 
#' The function extracts parent features, plus any tRNA/rRNA rows (retained for masking even when child features); it also excludes non-coding RNAs already annotated in the file.
#' 
#' @param annotation_file  GFF3 genome annotation file.
#' @param annot_file_directory The directory path for the annotation file (default is '.')
#' @param target_strand A character string indicating the strand. Supports two valies; '+' and '-'.
#' @param original_sRNA_annotation A string indicating how the biotype of pre-annotated ncRNA, which can be found in the attribute column.In case if the user does not know how the sRNA is annotated, it can be set as "unknown". In this case, all RNAs apart from tRNAs and rRNAs will be removed from the selection.
#' 
#' @return A dataframe with the major features for a set strand, plus any tRNA/rRNA rows retained for masking.
#' 
#' @importFrom utils read.delim
#' @export
major_features <- function(annotation_file, annot_file_directory = ".", target_strand, original_sRNA_annotation) {
  gff_cache <- .resolve_gff_cache(annotation_file, annot_file_directory)
  gff <- gff_cache$parsed
  ## Pre-annotated sRNAs always have a defined biotype, which can be found in the attribute column.
  ## The following code creates a regex that will recognise pre-annotated sRNAs
  ori_sRNA_biotype <- c()
  ## In case if the user does not know how the sRNA is annotated, it can be set as "unknown". In this case, all RNAs apart from tRNAs and rRNAs will be removed from the selection.
  if (original_sRNA_annotation=="unknown") {
    ori_sRNA_biotype <- "biotype=.*?[^tr]RNA;"
  } else {
    ori_sRNA_biotype <- paste0("biotype=", original_sRNA_annotation)
  }
  
  ## Select the major genomic features: remove child features (CDS, mRNA etc.), previously annotated sRNAs and extra features, but retain tRNA/rRNA for masking.
  is_trna_rrna <- gff[,3] %in% c("tRNA", "rRNA")
  major_f <- gff[(!grepl("Parent", gff[,9], ignore.case = TRUE) & gff[,3]!='chromosome' & gff[,3]!='biological_region' & !grepl(ori_sRNA_biotype, gff[,9], ignore.case = TRUE) & gff[,3]!='region' & gff[,3]!='sequence_feature') | is_trna_rrna,]
  ## Select only major features for the target strand.
  m_strand_features <- data.frame()
  if (target_strand=="+") {
    m_strand_features <- major_f[major_f[,7]=="+",]
  } else if (target_strand=="-") {
    m_strand_features <- major_f[major_f[,7]=="-",]
  } else {
    return(major_f)
  }
  
  return(m_strand_features)
  
}



#' sRNA prediction
#' 
#' The function for prediction and annotation of sRNA.
#' 
#' @param major_strand_features A dataframe containing the major features for a particular strand.
#' @param target_strand A character string indicating the strand. Supports two values; '+' and '-'.
#' @param union_peak_ranges An IRanges object containing genomic coordinated for all peaks detected on the target strand.
#' 
#' @return An IRanges object containing coordinates and names of the predicted sRNA.
#' 
#' @export
sRNA_calc <- function(major_strand_features, target_strand, union_peak_ranges) {
  ## This function predicts sRNAs.
  ## Convert strand feature coordinates into IRanges.
  strand_IRange <- IRanges(start = major_strand_features[,4], end = major_strand_features[,5])
  ## Select only the ranges that do not overlap the annotated features
  ## Also, disregard the ranges that finish/start 1 position before the genomic feature, because they should be considered as UTRs.
  IGR_sRNAs <- union_peak_ranges[IRanges::match(union_peak_ranges, subsetByOverlaps(union_peak_ranges, strand_IRange, maxgap = 1L), nomatch = 0) == 0,]
  ## Construct the IDs for the new sRNAs to be added into the attribute column of the annotation.
  if (target_strand=="+") {
    names(IGR_sRNAs) <- apply(as.data.frame(IGR_sRNAs),1, function(x) paste0("ID=putative_sRNA:p", x[1], "_", x[2], ";"))
  } else if (target_strand== "-") {
    names(IGR_sRNAs) <- apply(as.data.frame(IGR_sRNAs),1, function(x) paste0("ID=putative_sRNA:m", x[1], "_", x[2], ";"))
  } else {
    stop("Select strand")
  }
  return(IGR_sRNAs)
}



#' UTR prediction
#' 
#' Function for prediction and annotation of UTRs.
#' 
#' @param major_strand_features A dataframe containing the major features for a particular strand.
#' @param target_strand A character string indicating the strand. Supports two valies; '+' and '-'.
#' @param union_peak_ranges An IRanges object containing genomic coordinated for all peaks detected on the target strand.
#' @param min_UTR_length An integer indicating the minimum UTR length.
#' 
#' @return An IRanges object containing coordinates and names of the predicted UTRs.
#' 
#' @export
UTR_calc <- function(major_strand_features, target_strand, union_peak_ranges, min_UTR_length) {
  ## This function predicts UTRs.
  ## Convert strand feature coordinates into IRanges.
  strand_IRange <- IRanges(start = major_strand_features[,4], end = major_strand_features[,5])
  ## Find the peak union ranges that overlap with genomic features. Also, include the ranges that do not overlap the features but start/finish 1 position away from it.
  overapping_features <- subsetByOverlaps(union_peak_ranges, strand_IRange, maxgap = 1L)
  ## Join the overlapping features with the genomic ones for further cutting.
  overapping_features <- c(overapping_features, strand_IRange)
  ## Cut the overlapping features on teh border where they overlap with the genomic features.
  split_features <- disjoin(overapping_features)
  ## Now select only the UTR "overhangs" that are created by cutting overlapping features on the border.
  UTRs <- split_features[IRanges::match(split_features, subsetByOverlaps(split_features, strand_IRange), nomatch = 0) == 0]
  ## Select only UTRs that satisfy the minimum length condition.
  UTRs <- UTRs[width(UTRs)>=min_UTR_length,]
  ## Construct the IDs for the new UTRs to be added into the attribute column of the annotation.
  if (target_strand=="+") {
    names(UTRs) <- apply(as.data.frame(UTRs),1, function(x) paste0("ID=putative_UTR:p", x[1], "_", x[2],";"))
  } else if (target_strand== "-") {
    names(UTRs) <- apply(as.data.frame(UTRs),1, function(x) paste0("ID=putative_UTR:m", x[1], "_", x[2],";"))
  } else {
    stop("Select strand")
  }
  return(UTRs)
  
}



#' Strand annotation
#' 
#' This function constructs the full annotation for the strand.
#' 
#' @param target_strand A character string indicating the strand. Supports two values; '+' and '-'.
#' @param sRNA_IRanges An IRanges object containing coordinates and names of the predicted sRNAs.
#' @param UTR_IRanges An IRanges object containing coordinates and names of the predicted UTRs.
#' @param major_strand_features A dataframe containing the major features for a particular strand.
#' 
#' @return A dataframe containing strand annotation populated with the prediction features, build in accordance with GFF3 file format.
#' 
#' @export
strand_feature_editor <- function(target_strand, sRNA_IRanges, UTR_IRanges, major_strand_features) {
  ## Collect GFF attributes that carry no parseable ID, and any neighbour name that
  ## comes back purely numeric (a wrong-match signature), to report each once after
  ## the loop.
  unparsed_attrs <- character(0)
  numeric_ids    <- character(0)

  ## Isolate the ID= value, bounded to the first ';', so the capture can never run
  ## on into a later attribute such as Dbxref=GeneID:885041. Records the attribute
  ## when no ID= token is present.
  parse_id_token <- function(attr) {
    token <- sub("ID=([^;]*).*", "\\1", attr)
    if (token == attr) {
      unparsed_attrs <<- c(unparsed_attrs, attr)
    }
    token
  }
  ## Feature name: the part after the first ':' for an Ensembl type:name token, or
  ## the whole token for a hyphen or bare ID (RefSeq gene-Rv0001, Prokka PROKKA_0001).
  parse_name <- function(attr) {
    name <- sub("^[^:]*:", "", parse_id_token(attr))
    if (grepl("^[0-9]+$", name)) {
      numeric_ids <<- c(numeric_ids, name)
    }
    name
  }
  ## Feature type: the part before the first ':' within the ID token.
  parse_type <- function(attr) sub(":.*", "", parse_id_token(attr))
  
  ## Join the sRNA and UTR ranges together.
  sRNA_UTR <- c(sRNA_IRanges, UTR_IRanges)
  n_features <- length(sRNA_UTR)
  
  ## Create information to go into corresponding columns.
  seqid <- rep(major_strand_features[1,1], n_features)
  empty_col <- rep(".", n_features)
  ## Create a dataframe for sRNAs and UTRs with all GFF3 file columns.
  cmp_strand <- data.frame()
  if (target_strand=="+") {
    cmp_strand <- data.frame(seqid, empty_col, empty_col, as.integer(start(sRNA_UTR)), as.integer(end(sRNA_UTR)), empty_col, rep("+", n_features), empty_col, names(sRNA_UTR), stringsAsFactors = FALSE)
  } else if (target_strand== "-") {
    cmp_strand <- data.frame(seqid, empty_col, empty_col, as.integer(start(sRNA_UTR)), as.integer(end(sRNA_UTR)), empty_col, rep("-", n_features), empty_col, names(sRNA_UTR), stringsAsFactors = FALSE)
  } else {
    stop("Select strand")
  }
  
  names(cmp_strand) <-names(major_strand_features)
  ## Join the sRNA/UTR dataframe with the dataframe for the major features.
  cmp_strand <- rbind(cmp_strand, major_strand_features)
  ## Order the dataframe by the feature start position.
  cmp_strand <- cmp_strand[order(cmp_strand[,4]),]
  
  ## Set the previous feature name to be the ID of the last feature in the chromosome, accounting for the fact that bacterial genomes are circular.
  previous_feature_name <- parse_name(cmp_strand[nrow(cmp_strand),9])
  
  for (i in 1:nrow(cmp_strand)) {
    ## Determine feature type from the attribute column if the third column is empty.
    feature_name <- parse_name(cmp_strand[i,9])
    if (cmp_strand[i,3]==".") {
      feature_type <- parse_type(cmp_strand[i,9])
      cmp_strand[i,3] <- feature_type
      ## Find the name of the next feature in the annotation.
      next_feature_name <- c()
      if (i+1 <= nrow(cmp_strand)) {
        next_feature_name <- parse_name(cmp_strand[i+1,9])
      } else {
        next_feature_name <- parse_name(cmp_strand[1,9])
      }
      
      ## Build feature attribute column information for sRNAs and UTRs, including upstream and downstream features (with reagards to the strand).
      feature_attribute <- c()
      
      if (target_strand=="+") {
        feature_attribute <- paste0(cmp_strand[i,9],"upstream_feature=", previous_feature_name, ";downstream_feature=", next_feature_name)
      } else if (target_strand== "-") {
        feature_attribute <- paste0(cmp_strand[i,9],"upstream_feature=", next_feature_name, ";downstream_feature=", previous_feature_name)
      } else {
        stop("Select strand")
      }
      cmp_strand[i,9] <- feature_attribute
      
    }
    previous_feature_name <- feature_name
  }
  
  ## Report rows with no parseable ID once, as a single deduplicated summary.
  if (length(unparsed_attrs) > 0L) {
    n_failed <- length(unique(unparsed_attrs))
    warning(paste0(n_failed, " of ", nrow(cmp_strand),
                   " feature rows had no parseable ID= attribute (e.g. ",
                   unparsed_attrs[1], "); each feature row should carry an ID= tag."),
            call. = FALSE, immediate. = TRUE)
  }
  ## Report purely-numeric neighbour names once: a parse wrong-match signature.
  if (length(numeric_ids) > 0L) {
    warning(paste0(length(unique(numeric_ids)),
                   " feature name(s) parsed to a purely-numeric value (e.g. ",
                   numeric_ids[1], "); this is a signature of an ID-parsing wrong-match."),
            call. = FALSE, immediate. = TRUE)
  }
  return(cmp_strand)
  
}





#' Prediction and annotation of sRNAs and UTRs from RNA-seq data
#' 
#' A wrapper function that executes all prediction steps for each strand and builds the final GFF3 annotation.
#' 
#' @param bam_directory The directory containing BAM files.
#' @param bam_list Optional newline separated text file of filenames of bam files. File must be located in bam_location directory. 
#' @param original_annotation_file GFF3 genome annotation file.
#' @param annot_file_dir The directory containing the GFF3 annotation file.
#' @param output_file A string containing the name of an output file.
#' @param original_sRNA_annotation A string indicating how the biotype of pre-annotated ncRNA, which can be found in the attribute column.In case if the user does not know how the sRNA is annotated, it can be set as "unknown". In this case, all RNAs apart from tRNAs and rRNAs will be removed from the selection.
#' @param low_coverage_cutoff An integer indicating the low coverage threshold value.
#' @param high_coverage_cutoff An integer indicating the high coverage threshold value.
#' @param min_sRNA_length An integer indicating the minimum peak width/sRNA length.
#' @param min_UTR_length An integer indicating the minimum UTR length.
#' @param paired_end_data A boolean indicating if the reads are paired-end.
#' @param strandedness A string outlining the type of the sequencing library: stranded, or reversely stranded. Defaults to "stranded"; "unstranded" is rejected with an error.
#' @param scanbamparam An optional \code{Rsamtools::ScanBamParam} object giving
#'   full control of the BAM read filter. When \code{NULL} (the default) an
#'   internal filter is built from \code{mapqFilter} plus alignment-flag
#'   exclusions: unmapped, QC-failing, secondary and supplementary alignments
#'   are dropped, and for paired-end data the read must additionally be properly
#'   paired with a mapped mate. A supplied \code{scanbamparam} takes precedence
#'   and \code{mapqFilter} is then ignored.
#' @param mapqFilter Integer. Minimum mapping quality (MAPQ) a read must have to
#'   be retained during coverage construction, used only when
#'   \code{scanbamparam = NULL}. Default 10, which keeps uniquely-mapped reads
#'   on the common bacterial aligners while discarding the ambiguous and
#'   multi-mapping reads that cluster at MAPQ 0 to 3. MAPQ is aligner-specific,
#'   so do NOT set \code{mapqFilter} above the MAPQ your aligner assigns a
#'   uniquely-mapped read, or all reads are discarded and coverage is empty.
#'   Maximum recommended values by aligner: BWA-MEM 60, minimap2 60, Rsubread
#'   align/subjunc 40, BWA aln (backtrack) 37, Bowtie2 42 (but see note). Note:
#'   Bowtie2 uses a non-monotonic MAPQ in which some uniquely-mapped reads that
#'   carry mismatches score only 3 or 8, below the default of 10, so Bowtie2
#'   users should set \code{mapqFilter = 1}. Set \code{mapqFilter = NA} to
#'   disable mapping-quality filtering while keeping the flag exclusions. A
#'   read-quality filter that retains no reads triggers a warning.
#' @param coverage_model Character, one of "fragment" or "footprint", setting how
#'   paired-end reads contribute to coverage. It applies only to paired-end data
#'   and is ignored for single-end. "fragment" (the default) treats each read
#'   pair as one span from its leftmost to rightmost aligned base, so the
#'   unsequenced insert between the mates is counted as covered. "footprint"
#'   instead counts only the aligned blocks of the two mates, so the insert
#'   between them is left uncovered and the coverage reflects only the bases
#'   that were actually sequenced.
#' 
#' @return Outputs a new GFF3 file populated with predicted sRNAs and UTRs.
#'
#' 
#' @export
feature_file_editor <- function(bam_directory = ".", bam_list = "", original_annotation_file, annot_file_dir = ".", output_file, original_sRNA_annotation, low_coverage_cutoff, high_coverage_cutoff, min_sRNA_length, min_UTR_length, paired_end_data = FALSE, strandedness  = "stranded", scanbamparam = NULL, mapqFilter = 10, coverage_model = c("fragment", "footprint")) {
  ## Validate strandedness at entry, before any filesystem work. peak_union_calc()
  ## validates too, so a direct caller of that function is guarded as well.
  valid_strandedness <- c("stranded", "reversely_stranded")
  if (!strandedness %in% valid_strandedness) {
    stop("Invalid 'strandedness' value: '", strandedness,
         "'. Feature prediction requires stranded data; must be one of: ",
         paste(valid_strandedness, collapse = ", "),
         ". ('unstranded' is valid for read counting but not for prediction.)",
         call. = FALSE)
  }
  test <- list.files(path = bam_directory, pattern = "\\.BAM$", full.names = TRUE, ignore.case = TRUE)
  if (length(test) > 0){
    ## Load the original GFF once for the whole wrapper.
    gff_cache <- load_gff_cache(original_annotation_file, annot_file_dir)

    ## Plus strand
    peak_sets <- peak_union_calc(bam_location = bam_directory, bam_txt_list = bam_list, low_coverage_cutoff, high_coverage_cutoff, min_sRNA_length, paired_end_data, strandedness, scanbamparam = scanbamparam, mapqFilter = mapqFilter, coverage_model = coverage_model)
    plus_strand_peaks  <- peak_sets$plus
    message("Extracted plus strand data from BAM files")
    maj_plus_features <- major_features(gff_cache, annot_file_directory = annot_file_dir, "+", original_sRNA_annotation)
    plus_sRNA <- sRNA_calc(maj_plus_features, "+", plus_strand_peaks)
    plus_UTR <- UTR_calc(maj_plus_features, "+", plus_strand_peaks, min_UTR_length)
    ## tRNA/rRNA are masked (above) but dropped from the neighbour set, so neighbour fields are unchanged.
    plus_neighbour_features <- maj_plus_features[!(maj_plus_features[,3] %in% c("tRNA", "rRNA")), ]
    plus_annot_dataframe <- strand_feature_editor("+", plus_sRNA, plus_UTR, plus_neighbour_features)
    message("Built plus strand annotation dataframe")
    ## Minus strand
    minus_strand_peaks <- peak_sets$minus
    message("Extracted minus strand data from BAM files")
    maj_minus_features <- major_features(gff_cache, annot_file_directory = annot_file_dir, "-", original_sRNA_annotation)
    minus_sRNA <- sRNA_calc(maj_minus_features, "-", minus_strand_peaks)
    minus_UTR <- UTR_calc(maj_minus_features, "-", minus_strand_peaks, min_UTR_length)
    minus_neighbour_features <- maj_minus_features[!(maj_minus_features[,3] %in% c("tRNA", "rRNA")), ]
    minus_annot_dataframe <- strand_feature_editor("-", minus_sRNA, minus_UTR, minus_neighbour_features)
    message("Built minus strand annotation dataframe")
  
    ## Creating the final annotation dataframe by combining both strand dataframe and adding missing information like child features from the original GFF3 file.
    annotation_dataframe <- rbind(gff_cache$parsed, plus_annot_dataframe, minus_annot_dataframe)
    ## Remove all teh repeating information.
    annotation_dataframe <- unique(annotation_dataframe)
    ## Order the dataframe by feature start coordinates.
    annotation_dataframe <- annotation_dataframe[order(annotation_dataframe[,4]),]
    message("Prepared complete annotation dataframe")
  
    ## Restore the original header from the cache.
    f <- gff_cache$raw_lines
    header <- c()
    i <- 1
    while (grepl("#", f[i], fixed = TRUE)) {
      f_line <- f[i]
      header <- c(header,f_line)
      i <- i+1
    }
    # add a line to indicate the origin of the file (single # commas should be ignored by programs)
    header <- c(header, "# produced by baerhunter")
  
    message("Building output file now")
  
    ## Create the final GFF3 file.
    write.table(header, output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(annotation_dataframe, output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
  
    invisible(output_file)
  }else{
    stop("No BAMs in bam directory!")
  }
}

#commit1 completed
#commit2 completed
#commit3 completed
#commit4 completed
#commit5 completed
#commit6 completed
#commit7 completed
#commit8 completed
#commit9 completed
#commit10 completed
