#' TPM calculator
#'
#' This function uses feature count tables to calculate TPM values for each gene and sample.
#'
#' @param count_table A CSV file containing feature counts for each sample.
#' @param complete_ann A GFF3 annotation file or SAF dataframe
#' @param is_gff A boolean indicating whether annotation is gff file, default=TRUE
#' @param excl_rna A boolean indicating if misc RNA features (rRNA, tRNA) are excluded from normalisation. (Defaults=TRUE)
#' @param output_file A string indicating the name of the output file.
#'
#' @return A dataframe with TPM values for each gene and sample; the same is written into the output file.
#'
#'
#' @importFrom utils read.delim write.table
#'
#' @export
tpm_normalisation <- function(count_table, complete_ann, is_gff = TRUE, output_file = NA, excl_rna = TRUE) {
  ## check output directory exists
  if (!is.na(output_file)) {
    out_dir <- dirname(output_file)
    stopifnot("Output directory doesn't exist." = dir.exists(out_dir))
  }
  ## make saf from gff (uses make_saf function)
  if (is_gff){
    nsaf_df <- make_saf(ann_file=complete_ann, exclude = excl_rna)
  }else{
    nsaf_df <- complete_ann
  }
  feature_names <- nsaf_df$GeneID

  ## Load in the count table
  count_df <- read.delim(count_table)

  ## Calculate the length of the features, if they are in the right order.
  stopifnot("Wrong feature order. Is excl_rna param set the same as for count_features?" = all(rownames(count_df) == feature_names))
  feature_lengths <- c()
  feature_lengths <- (nsaf_df$End-nsaf_df$Start+1)/1000

  ## Calculate RPK by dividing the feature count of each gene (per feature) by its length in kilobases.
  rpk_df <- count_df / feature_lengths
  ## Calculate the scaling factor for each sample by summing up all RPKs per sample and dividing by a million.
  sample_rpk_sum <- colSums(rpk_df)
  scaling_fact <- sample_rpk_sum/1000000

  ## Divide all RPK by the corresponding sampling factor.
  tpm_df <- sweep(rpk_df, 2, scaling_fact, "/")

  colnames(tpm_df) <- colnames(count_df)
  rownames(tpm_df) <- feature_names

  ## If the output file is set by the user, write the TPM dataframe into it.
  if (!is.na(output_file)) {
    write.table(tpm_df, output_file, sep = "\t", quote = FALSE)
  }
  return(tpm_df)
}


#' Flagging features depending on TPM value profile
#'
#' A helper function to analyse each row of the TPM table. Each feature gets allocated a flag depending on the expression profile.
#'
#' @param tpm_data A CSV file containing TPM values for each normalised feature in each sample.
#' @param complete_annotation A GFF3 annotation file.
#' @param output_file A string indicating the name of the output file.
#'
#' @return The path to the output GFF3 file, returned invisibly. The written file is the input annotation with an expression flag added to the attribute column of each flagged feature.
#'
#'
#' @importFrom utils read.delim write.table
#' @export
tpm_flagging <- function(tpm_data, complete_annotation, output_file) {


  tpm_analyser <- function(num_vec) {
    if (all(num_vec<0.5)) {
      return("expression_below_cutoff")
    } else  if (any(num_vec > 1000)) {
      return("high_expression_hit")
    } else if (any((num_vec > 10) & (num_vec <= 1000))) {
      return("medium_expression_hit")
    } else if (any((num_vec >= 0.5) & (num_vec <= 10))) {
      return("low_expression_hit")
    }
  }
  ## Read the TPM table and create a flag vector in the corresponding gene order.
  norm_data <- read.delim(tpm_data, header = TRUE)
  # return if the data frame contains non-numeric data (comparisons will fail otherwise)
  if (!is.numeric(as.matrix(norm_data))) {
    stop("Input file contains unexpected non-numeric data\n")
  }
  flags <- apply(norm_data, 1, function(x) tpm_analyser(as.vector(x)))
  flag_names <- names(flags)
  ann_file <- readLines(complete_annotation)
  ## Add the flag to the corresponding feature's attribute column.
  feature_names <- sub(".*?ID=(.*?);.*", "\\1", ann_file)
  matched <- feature_names %in% flag_names
  if (!any(matched)) {
    warning("No annotation feature IDs matched the TPM table; check that the GFF ID format is consistent with the count table.",
            call. = FALSE, immediate. = TRUE)
  }
  new_annot <- ann_file
  new_annot[matched] <- paste0(ann_file[matched], ";expression_flag=",
                               flags[feature_names[matched]])

  write.table(new_annot, output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  invisible(output_file)
}


#' Filtering features by a target flag.
#'
#' A function to filter the marked features selected by the user by the flag of choice.
#'
#' @param flagged_annotation_file A flagged GFF3 annotation file.
#' @param target_features A string indicating feature type to filter (optional).
#' @param target_flag A string indicating a flag for filtering.
#' @param output_file A string indicating the name of the output file.
#'
#' @return A GFF3 file where a target feature type is filtered by the expression flag of interest.
#'
#'
#' @importFrom utils read.delim write.table
#' @export
tpm_flag_filtering <- function(flagged_annotation_file, target_features = c("putative_sRNA", "putative_UTR"), target_flag, output_file) {

  ##load in annotation data.
  annot_data <- read.delim(flagged_annotation_file, header = FALSE, comment.char = "#")

  is_target  <- annot_data[, 3] %in% target_features
  flag_match <- grepl(target_flag, annot_data[, 9], fixed = TRUE)
  keep       <- !is_target | flag_match
  filtered_selection <- annot_data[keep, ]
  ## Restore the original header.
  f <- readLines(flagged_annotation_file)
  header <- c()
  i <- 1
  while (grepl("#", f[i], fixed = TRUE)) {
    header <- c(header,f[i])
    i <- i+1
  }
  ## Write a new GFF3 file.
  write.table(header, output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(filtered_selection, output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)

}

#commit1 completed
#commit2 completed
#commit3 completed
#commit4 completed
#commit5 completed
#commit6 completed
