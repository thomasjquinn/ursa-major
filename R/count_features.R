#' read_annotation_file function
#'
#' This function pastes together the path and filename of annotation file (if not in current directory) and tests file for existence.
#'
#' @param annot_dir Directory where genome annotation file is located
#' @param annot_file  GFF3 or GTF genome annotation file
#'
#' @return annot_file_loc Complete path of existing annotation file
#'
#' @export
read_annotation_file <- function(annot_dir, annot_file){
  annot_file_loc <- file.path(annot_dir, annot_file)
  stopifnot("Annotation file not found" = file.exists(annot_file_loc))
  return(annot_file_loc)
}


#' Load a GFF annotation into a reusable cache
#'
#' Reads a GFF3 annotation file once and returns a list holding both the
#' raw lines and the parsed dataframe, so downstream baerhunter functions
#' can reuse a single parse instead of re-reading the file from disk. The
#' path is resolved the same way the prediction and counting stages resolve
#' annotation paths: a complete path is used as given, and a directory is
#' prepended only when \code{annot_file_directory} is not ".". The resolved
#' path is checked for existence before it is read.
#'
#' @param annotation_file A GFF3 genome annotation file, given either as a
#'   complete path or as a filename to be resolved against
#'   \code{annot_file_directory}.
#' @param annot_file_directory The directory containing the annotation file
#'   (default "."). When ".", \code{annotation_file} is used as given.
#'
#' @return A list with three named slots: \code{path} (the resolved file
#'   path), \code{raw_lines} (the annotation as a character vector, read by
#'   \code{readLines}) and \code{parsed} (the tab-delimited dataframe, read
#'   by \code{read.delim}).
#'
#' @importFrom utils read.delim
#' @export
load_gff_cache <- function(annotation_file, annot_file_directory = ".") {
  annot_file_loc <- if (annot_file_directory == ".") {
    annotation_file
  } else {
    file.path(annot_file_directory, annotation_file)
  }
  stopifnot("Annotation file not found" = file.exists(annot_file_loc))
  list(
    path      = annot_file_loc,
    raw_lines = readLines(annot_file_loc),
    parsed    = read.delim(annot_file_loc, header = FALSE,
                           comment.char = "#")
  )
}


#' Resolve a GFF argument to a cache
#'
#' Internal dispatch helper. Returns its argument unchanged if it is already
#' a GFF cache (a list carrying the \code{path}, \code{raw_lines} and
#' \code{parsed} slots); otherwise treats it as a file path and builds a
#' cache via \code{load_gff_cache}. This lets every public function that
#' accepts a GFF path also accept a pre-built cache without a signature
#' change.
#'
#' @param x Either a GFF cache or a path to a GFF3 annotation file.
#' @param annot_file_directory The directory containing the annotation file (default ".").
#'
#' @return A GFF cache list (see \code{load_gff_cache}).
#'
#' @keywords internal
.resolve_gff_cache <- function(x, annot_file_directory = ".") {
  if (is.list(x) &&
      all(c("path", "raw_lines", "parsed") %in% names(x))) {
    return(x)
  }
  load_gff_cache(x, annot_file_directory)
}


#' find_strandedness function
#'
#' This function translates the user-inputted strandedness parameter to the required integer input for strandSpecific arg of featureCounts.
#'
#' @param strand_param user input 'stranded' or 'reversely-stranded'
#'
#' @return strand_sp integer value
#'
#' @export
find_strandedness <- function(strand_param){
  strand_sp <- integer()
  strand_strings <- c("unstranded", "reversely_stranded", "stranded")
  stopifnot("Invalid strandedness parameter: must be either 'unstranded', 'stranded' or 'reversely_stranded'" = strand_param %in% strand_strings)
  if(strand_param=="unstranded"){
    strand_sp <- 0
  } else if(strand_param=="stranded"){
    strand_sp <- 1
  } else if(strand_param=="reversely_stranded"){
    strand_sp <- 2
  }
  return(strand_sp)
}


#' SAF converter
#'
#' This function converts gff file into simplified annotation format with 5 columns, appropriate for use in featureCounts/tpm_norm_flagging.
#'
#' @param ann_file A GFF3 annotation file.
#' @param exclude A boolean to indicate whether or not to include rRNA/tRNA features
#'
#' @return A dataframe in Simplified Annotation Format
#'
#' @importFrom stringr str_match
#' @export
make_saf <- function(ann_file, exclude=FALSE){
  gff_cache <- .resolve_gff_cache(ann_file)
  gff <- gff_cache$parsed
  ## check that file is correct format (9 cols)
  stopifnot("annotation file format is invalid" = ncol(gff)==9)
  ## Select only the major genomic features: remove all child features (like CDS, mRNA etc.) and extra features
  if (!exclude){
    major_f <- gff[!grepl("Parent", gff[,9], ignore.case = TRUE) &
                     gff[,3]!='chromosome' & gff[,3]!='biological_region' &
                     gff[,3]!='region' & gff[,3]!='sequence_feature',]
  }else{ # if you want to exclude rRNA and tRNA features)
    major_f <- gff[!grepl("Parent", gff[,9], ignore.case = TRUE) &
                     gff[,3]!='chromosome' & gff[,3]!='biological_region' &
                     gff[,3]!='region' &
                     gff[,3]!='sequence_feature' &
                     gff[,3]!='tRNA_gene' &
                     gff[,3]!='rRNA_gene',]
  }
  saf_df <- as.data.frame(matrix(0, ncol=5, nrow=nrow(major_f)))
  colnames(saf_df) <- c("GeneID", "Chr", "Start", "End", "Strand")
  saf_df$GeneID   <- str_match(major_f$V9, "ID=(.*?);")[,2]
  saf_df$Chr      <- major_f$V1
  saf_df$Start    <- major_f$V4
  saf_df$End      <- major_f$V5
  saf_df$Strand   <- major_f$V7
  return(saf_df)
}

#' count_features function
#'
#' This is a function to employ Rsubread featureCounts to quantify expression of annotated and predicted elements
#'
#' @param bam_dir The directory where bam files located
#' @param annotation_dir The directory where annotation file is located
#' @param annotation_file The complete annotation file: GFF3 or GTF genome annotation file
#' @param output_dir The full directory path for CSV output files to be written
#' @param output_filename The name for the output files--for example dataset name
#' @param chromosome_alias_file A comma-delimited TXT file containing a character string with the chromosome names. This file has to have two columns: first with the chromosome name in the annotation file, second with the chromosome name in the BAM file.
#' @param strandedness A string outlining the type of the sequencing library: unstranded, stranded, or reversely stranded.
#' @param is_paired_end A boolean indicating if the reads are paired-end.
#' @param excl_rna A boolean indicating if misc RNA features (rRNA, tRNA) are excluded from quantification. (Defaults=TRUE)
#' @param largest_overlap A boolean; if TRUE, assigns each read to the feature with the largest number of overlapping bases. Maps to featureCounts largestOverlap. Combined with the package's `fraction = TRUE`, Rsubread >= 2.14.0 is recommended, since earlier versions silently miscount that combination. (Default: FALSE)
#' @param frac_overlap_feature Minimum fraction of a feature that must be overlapped before a read is assigned to it. Maps to featureCounts fracOverlapFeature. (Default: 0)
#' @param read_to_pos Reduce each read to a single base before counting: 5 for the 5' end, 3 for the 3' end, or NULL to count the whole read. Maps to featureCounts read2pos. (Default: NULL)
#' @param count_multi_mapping_reads A boolean; if FALSE, reads mapping to multiple locations are excluded from counts. Maps to featureCounts countMultiMappingReads. (Default: FALSE, reproducing 2019 behaviour)
#' @param count_read_pairs A boolean; for paired-end data, if TRUE each fragment is counted once rather than each mate separately. Ignored for single-end data. Maps to featureCounts countReadPairs. Requires Rsubread >= 2.4.3. (Default: TRUE)
#' @param ... Optional parameters passed on to featureCounts(). Note that allowMultiOverlap and fraction are set internally to TRUE and cannot be overridden via ...
#'
#' @return Count tables for each feature are written into separate files, as well as the result summary.
#'
#' @import Rsubread
#' @importFrom tools file_path_sans_ext
#' @importFrom utils write.table
#' @export
count_features <- function(bam_dir=".",
                           annotation_dir=".",
                           annotation_file,
                           output_dir=".",
                           output_filename="dataset",
                           chromosome_alias_file,
                           strandedness,
                           is_paired_end,
                           excl_rna = TRUE,
                           largest_overlap = FALSE,
                           frac_overlap_feature = 0,
                           read_to_pos = NULL,
                           count_multi_mapping_reads = FALSE,
                           count_read_pairs = TRUE,
                           ...){
  ## function to call rsubread featureCounts

  ## Load the annotation once (path resolved, existence checked, parsed and raw lines cached).
  gff_cache <- load_gff_cache(annotation_file, annotation_dir)

  ## Compile a list of BAM files present in the bam directory
  bam_files <- list.files(path = bam_dir, pattern = "\\.BAM$", full.names = TRUE, ignore.case = TRUE)
  stopifnot("Empty bam directory" = length(bam_files) > 0)

  ## Convert to SAF from the cache
  nsaf_df <- make_saf(ann_file = gff_cache, exclude = excl_rna)

  ## if output directory exists, create filenames and path to output file
  stopifnot("Output directory doesn't exist" = dir.exists(output_dir))
  output_file  <- paste(output_dir, output_filename, sep = "/")
  count_file_name <- paste0(output_file, "_Counts.csv")
  summary_file_name <- paste0(output_file, "_Count_summary.csv")

  ## The strandedness set by the user is translated into the integer for strandSpecific argument of the featureCounts function.
  strand_specific <- find_strandedness(strandedness)

  ## Paired-end is FALSE by default unless specified by the user otherwise.
  paired_end <- isTRUE(is_paired_end)

  ## Extract BAM file names without extension.
  sample_names <- c(file_path_sans_ext(basename(bam_files)))

  ## Execute featureCounts function. Feature counts and summary stats are written into separate TAB delimited files.
  fc <- featureCounts(bam_files,
                      annot.ext = nsaf_df,
                      chrAliases = chromosome_alias_file,
                      strandSpecific = strand_specific,
                      isPairedEnd = paired_end,
                      # countReadPairs requires Rsubread >= 2.4.3 (present in the 2.4.3 reference manual, RELEASE_3_12, 30 March 2021);
                      # forwarded unconditionally, with the matching DESCRIPTION Imports floor added at the documentation commit.
                      countReadPairs = count_read_pairs,
                      countMultiMappingReads = count_multi_mapping_reads,
                      allowMultiOverlap = TRUE,
                      fraction = TRUE,
                      largestOverlap = largest_overlap,
                      fracOverlapFeature = frac_overlap_feature,
                      read2pos = read_to_pos,
                      ...)
  colnames(fc$counts) <- sample_names
  write.table(fc$counts, count_file_name, sep = "\t", quote = FALSE)
  colnames(fc$stat) <- c('Status',sample_names)
  write.table(fc$stat, summary_file_name, sep = "\t", quote = FALSE)

  invisible(NULL)
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
