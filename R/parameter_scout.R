# =============================================================================
# parameter_scout.R
#
# Advisory tool. Measures per-base coverage in the intergenic regions of each
# BAM, reports its percentiles, draws one figure of the distributions, and
# offers two cutoff pairs for low_coverage_cutoff and high_coverage_cutoff.
# It sets no parameter and nothing else in baerhunter calls it.
#
# Exported: parameter_scout(), plot_scout_distribution(), suggest_cutoffs().
# Set paired_end_data to match the library.
#
# Needs count_features.R for .resolve_gff_cache() and feature_file_editor.R
# for major_features(). The BAM-reading helpers are copies of that file's, so
# the percentiles are on the scale peak_union_calc() later slices.
# =============================================================================


## --- section: internal constants -------------------------------------------

#' Maximum BAM files per scan; the reasoning is with \code{.check_bam_limit()}
#'
#' @noRd
.BH_MAX_BAMS <- 9L


## --- section: figure configuration -----------------------------------------
##
## Named constants for the figure. The low cutoff moves the slice boundaries,
## so raising it can yield more features as well as shorter ones; the high
## cutoff only sets how much evidence a slice needs, so raising it always
## yields a strict subset of the lower setting's output. Only the second is
## monotone, so a reader who takes the stringent pair to mean fewer
## predictions will sometimes be wrong.

## Fill palette, decorative only. One of "Okabe-Ito" (8 hues, colour-vision
## safe), "Dark2", "Set2", "Paired", "viridis" or "grey".
PALETTE <- "Okabe-Ito"

## Fill transparency.
FILL_ALPHA <- 0.75

## "area" gives every violin equal area; "max" gives them equal height.
VIOLIN_SCALE <- "area"

## Kernel bandwidth on log2 coverage; NULL uses Silverman's weighted rule.
BW_OVERRIDE <- NULL

## Device settings.
FIG_WIDTH   <- 9      # inches
FIG_HEIGHT  <- 7      # inches
FIG_RES     <- 300    # dpi, PNG only
TABLE_FRAC  <- 0.32   # fraction of figure height given to the table

## The aggregate row: the median across BAMs of each percentile column.
AGG_ROW_LABEL  <- "Median across BAMs"
AGG_ROW_GAP    <- 0.35     # row-heights of clear space above the aggregate row
AGG_RULE_COL   <- "grey20" # darker than the header rule's grey40
AGG_RULE_WIDTH <- 0.7

## --- the annotation block on the right of the plot -------------------------

## TRUE renders the annotation block with plotmath, which is what allows bold.
BOLD_AGG <- TRUE

## The lead sentence and the width it wraps to.
ANN_LEAD <- paste("Possible Cutoff Parameters: Use the shape of the",
                  "distributions to inform your choice.")
WRAP_CHARS <- 52

## Text size for the block, and one size down for the note.
ANN_SIZE      <- 3.1
ANN_NOTE_SIZE <- 2.7

## Hanging indent in characters, as an x offset: plotmath drops leading space.
ANN_EXDENT <- 0

## NULL leaves the note unwrapped; a number wraps it to that width.
NOTE_WRAP <- NULL
NOTE_MAX_CHARS <- 70

## Horizontal anchor of the block, a fraction of the x range; NULL computes it.
ANN_X_FRAC    <- 0.62
ANN_RIGHT_FRAC <- 0.98

## Vertical start of the block, in line heights. It does not move the tags.
ANN_TOP_OFFSET <- 5.4

## Headroom above the marker tags, as a fraction of the plotted range. The y
## view is pinned by coord_cartesian() and the block is drawn into this.
ANN_HEADROOM_MULT <- 0.26

## Width of one character at ANN_SIZE, as a fraction of panel width.
ANN_CHAR_FRAC <- 0.0091

## Vertical spacing, as a fraction of the plotted y range.
ANN_LINE_FRAC <- 0.050    # between consecutive lines
ANN_GAP_FRAC  <- 0.027    # extra space above each pair line and the note

## Percentiles given a dashed marker; not the offered set.
MARKER_PCTS <- c("P25", "Median", "P75")

## The note below the two pairs. Its wording does not vary with the data.
ANN_NOTE      <- "Note: Try these first. Other cutoffs may suit your data better."
SHOW_ANN_NOTE <- TRUE

## The two offered pairs: Relaxed (Median / P75) and Stringent (P75 / P80).
## P75 is the high of one and the low of the other, so choosing P75 for both
## is not a valid pair. Values are looked up by percentile name, so the labels
## and the numbers cannot drift apart.
PAIR_LABELS <- c(relaxed = "Relaxed", stringent = "Stringent")
PAIR_LOW    <- c(relaxed = "Median",  stringent = "P75")
PAIR_HIGH   <- c(relaxed = "P75",     stringent = "P80")

## Line style for the percentile markers; all markers share one style.
CUT_LINETYPE <- "dashed"

## Geometry of one file's row: violin above the baseline, box below.
VIOLIN_H <- 0.72
BOX_H    <- 0.16
BOX_GAP  <- 0.14


## --- section: internal helpers, BAM reading --------------------------------

## Copies of the BAM reading in peak_union_calc() (feature_file_editor.R).
## Their bodies must stay identical to it.

#' List the BAM files for a run
#'
#' Internal. From a newline-separated list file, or by scanning a directory.
#'
#' @param bam_location The directory containing BAM files.
#' @param bam_txt_list Optional newline separated text file of BAM filenames.
#' @return A character vector of BAM file paths.
#' @noRd
.list_bam_files <- function(bam_location = ".", bam_txt_list = "") {
  if (bam_txt_list != "") {
    bam_files <- readLines(bam_txt_list)
    bam_files <- lapply(bam_files, function(x) paste(bam_location, x, sep = "/"))
  } else {
    bam_files <- list.files(path = bam_location, pattern = "\\.BAM$",
                            full.names = TRUE, ignore.case = TRUE)
  }
  bam_files
}

#' Build the default BAM read filter
#'
#' Internal. Drops unmapped, QC-failing, secondary and supplementary
#' alignments; paired-end reads must also be properly paired with a mapped mate.
#'
#' @param paired_end_data A boolean indicating if the reads are paired-end.
#' @param mapqFilter Integer. Minimum mapping quality; see
#'   \code{\link{peak_union_calc}}.
#' @return A \code{ScanBamParam}.
#' @noRd
.default_scanbamparam <- function(paired_end_data, mapqFilter = 10) {
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
  ScanBamParam(flag = scanbamflag, mapqFilter = mapqFilter)
}

#' Read one BAM and split its reads by strand
#'
#' Internal. Reads one BAM under the supplied filter, applies the paired-end
#' coverage model, and splits reads by strand, inverting for single-end
#' reversely-stranded libraries. Returns reads, not coverage: the caller runs
#' \code{coverage()} itself, which keeps \code{compute_strand_peaks()} inside
#' \code{peak_union_calc()} unchanged.
#'
#' strandMode 1 takes the pair strand from the first mate ("stranded"), 2 from
#' the last mate ("reversely_stranded": dUTP, Illumina TruSeq).
#'
#' @param f Path to a BAM file.
#' @param paired_end_data A boolean indicating if the reads are paired-end.
#' @param strandedness "stranded" or "reversely_stranded".
#' @param scanbamparam A \code{ScanBamParam} object.
#' @param coverage_model "fragment" or "footprint".
#' @return A list with \code{plus_reads}, \code{minus_reads}, \code{n_reads}
#'   (alignments retained by the filter) and \code{seqinfo}.
#' @noRd
.strand_split_reads <- function(f, paired_end_data, strandedness,
                                scanbamparam, coverage_model) {
  if (paired_end_data) {
    strand_mode <- if (strandedness == "reversely_stranded") 2 else 1
    read_pairs <- readGAlignmentPairs(f, strandMode = strand_mode, param = scanbamparam)
    ## "fragment" counts the whole pair, "footprint" only the aligned blocks.
    ## Built as GRanges here because coverage() on GAlignmentPairs is footprint.
    if (coverage_model == "footprint") {
      file_alignment <- unlist(grglist(read_pairs))
    } else {
      file_alignment <- granges(read_pairs)
    }
  } else {
    file_alignment <- readGAlignments(f, param = scanbamparam)
  }
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
  list(plus_reads  = plus_reads,
       minus_reads = minus_reads,
       n_reads     = length(file_alignment),
       seqinfo     = seqinfo(file_alignment))
}


## --- section: internal helpers, module -------------------------------------

#' The offered cutoff pairs, resolved from the constants
#'
#' Internal. The single definition of what the module offers; the figure and
#' the writers both read it.
#'
#' @param percentiles The per-BAM percentile table.
#' @return A data frame of one row per pair: \code{pair}, \code{low},
#'   \code{high}, \code{from}.
#' @noRd
.cutoff_pairs <- function(percentiles) {
  need <- unique(c(PAIR_LOW, PAIR_HIGH))
  stopifnot("the percentile table lacks a column the offered pairs need" =
              all(need %in% names(percentiles)))
  agg <- vapply(need, function(n) .agg_percentile(percentiles[[n]]), integer(1))
  data.frame(pair = unname(PAIR_LABELS),
             low  = unname(agg[PAIR_LOW]),
             high = unname(agg[PAIR_HIGH]),
             from = paste(unname(PAIR_LOW), "/", unname(PAIR_HIGH)),
             stringsAsFactors = FALSE)
}


#' Intergenic regions, by strand
#'
#' Internal. The per-strand complement of the major annotated features. tRNA
#' and rRNA are retained as features by \code{major_features()} and so are
#' masked out of the intergenic space; tmRNA, RNase P RNA and SRP RNA are not.
#' The reference name and length come from the BAM header, and
#' \code{seqlengths} is set so \code{gaps()} runs to the chromosome end.
#'
#' @param gff_cache A GFF cache (see \code{load_gff_cache}).
#' @param original_sRNA_annotation Biotype of pre-annotated ncRNA, or "unknown".
#' @param bam_seqinfo A \code{Seqinfo} from the BAM header (single sequence).
#' @return A list with slots \code{plus} and \code{minus}, each a GRanges.
#' @noRd
.igr_regions <- function(gff_cache, original_sRNA_annotation, bam_seqinfo) {
  stopifnot("parameter scout expects a single reference sequence" =
              length(seqnames(bam_seqinfo)) == 1)
  seqname   <- as.character(seqnames(bam_seqinfo))
  seqlength <- unname(seqlengths(bam_seqinfo))
  ## Without a seqlength gaps() would silently drop the terminal region.
  stopifnot("BAM header carries no sequence length; cannot define intergenic regions safely" =
              !is.na(seqlength))

  ## "." returns the full major-feature set, both strands.
  major_f <- major_features(gff_cache, annot_file_directory = ".", ".",
                            original_sRNA_annotation)
  ## GFF3 permits "." in the strand column; keep only stranded features.
  major_f <- major_f[major_f[, 7] %in% c("+", "-"), ]
  stopifnot("No stranded major features found in the annotation" = nrow(major_f) > 0)

  ## Positions from the GFF, seqname from the BAM, so the gaps key to coverage.
  feat <- GenomicRanges::GRanges(
    seqnames   = seqname,
    ranges     = IRanges(start = as.integer(major_f[, 4]),
                         end   = as.integer(major_f[, 5])),
    strand     = major_f[, 7],
    seqlengths = stats::setNames(seqlength, seqname)
  )
  ## gaps() also returns a "*"-strand complement spanning the whole sequence;
  ## taking only "+" and "-" discards it.
  g <- GenomicRanges::gaps(feat)
  list(plus  = g[strand(g) == "+"],
       minus = g[strand(g) == "-"])
}


#' Intergenic coverage values as a single Rle
#'
#' Internal. Restricts a strand's coverage to the intergenic regions, without
#' expanding to a per-base vector.
#'
#' @param reads A GAlignments (or GRanges) for one strand.
#' @param igr_granges A GRanges of intergenic regions on that strand.
#' @return An Rle of the coverage values inside \code{igr_granges}.
#' @noRd
.igr_coverage_rle <- function(reads, igr_granges) {
  cvg <- coverage(reads)                 # RleList keyed by seqname; BAM extends to seqlength
  if (length(igr_granges) == 0L) return(Rle(integer(0)))
  restricted <- cvg[igr_granges]         # RleList, one element per intergenic range
  unlist(restricted, use.names = FALSE)  # one Rle over all intergenic positions
}


#' Weighted quantiles of a run-length encoded vector
#'
#' Internal. Computes quantiles from run values and run lengths, reproducing
#' \code{stats::quantile(type = 7)} without expanding the vector.
#' Zero-coverage positions are excluded.
#'
#' @param x An Rle of coverage values.
#' @param probs Numeric vector of probabilities.
#' @return A named numeric vector of quantiles, with attribute
#'   \code{n_positions} (the number of non-zero positions).
#' @noRd
.rle_weighted_quantile <- function(x, probs) {
  v <- as.numeric(runValue(x))
  w <- as.numeric(runLength(x))
  keep <- v > 0
  v <- v[keep]; w <- w[keep]
  nm <- paste0(probs * 100, "%")
  if (length(v) == 0L) {
    out <- stats::setNames(rep(NA_real_, length(probs)), nm)
    attr(out, "n_positions") <- 0
    return(out)
  }
  o  <- order(v); v <- v[o]; w <- w[o]
  cw <- cumsum(w)
  n  <- cw[length(cw)]
  ## X[k], the k-th value of the sorted expanded vector, without expanding it.
  pick <- function(k) v[findInterval(k - 1, cw) + 1L]
  ## stats::quantile type 7.
  h    <- (n - 1) * probs + 1
  lo   <- floor(h)
  frac <- h - lo
  hi   <- pmin(lo + 1, n)
  out  <- stats::setNames(pick(lo) + frac * (pick(hi) - pick(lo)), nm)
  attr(out, "n_positions") <- n
  out
}


#' Column names for the percentile table
#' @noRd
.pct_colnames <- function(probs) {
  vapply(probs, function(p) {
    if (isTRUE(all.equal(p, 0.5))) "Median" else paste0("P", round(p * 100))
  }, character(1))
}


#' Run-length summary of a coverage Rle, for plotting
#'
#' Internal. Collapses a coverage Rle to (value, weight) rows, dropping zeros,
#' which is what a weighted violin needs without materialising every position.
#'
#' @noRd
.rle_summary <- function(x, file, strand) {
  v <- as.numeric(runValue(x))
  w <- as.numeric(runLength(x))
  keep <- v > 0
  if (!any(keep)) {
    return(data.frame(File = character(0), Strand = character(0),
                      Coverage = numeric(0), Weight = numeric(0),
                      stringsAsFactors = FALSE))
  }
  agg <- stats::aggregate(list(Weight = w[keep]),
                          by = list(Coverage = v[keep]), FUN = sum)
  data.frame(File = file, Strand = strand,
             Coverage = agg$Coverage, Weight = agg$Weight,
             stringsAsFactors = FALSE)
}

## --- section: .agg_percentile(), the single aggregation rule ---------------

#' Aggregate one percentile column across BAMs
#'
#' Internal. The single cross-BAM aggregation rule: the median across files of
#' one percentile column, taken up to the next whole number. Both
#' \code{suggest_cutoffs()} and the figure call it, so they cannot disagree.
#'
#' \code{ceiling()} rather than \code{round()}: the cutoffs are a floor, so
#' rounding up is the conservative direction, and base \code{round()} is
#' round-half-to-even, which would make the result depend on parity.
#' \code{na.rm = TRUE} so one empty BAM does not propagate NA.
#'
#' @param v A numeric vector: one percentile column, one value per BAM.
#' @return A single integer.
#' @noRd
.agg_percentile <- function(v) {
  as.integer(ceiling(stats::median(as.numeric(v), na.rm = TRUE)))
}

## --- section: .check_bam_limit(), the single BAM-count ceiling -------------

#' Enforce the BAM-count ceiling
#'
#' Internal. Called before any BAM is opened and again before the figure is
#' drawn, since a caller can reach the figure from a saved percentile table.
#' Nine, because past it the annotation block overflows the panel and the
#' violins are too compressed to read; the table and the cutoffs are sound at
#' any count.
#'
#' @param n Integer. Number of BAM files, or rows in the percentile table.
#' @param context "scan" or "figure", naming where the count came from.
#' @return Invisibly TRUE, or stops.
#' @noRd
.check_bam_limit <- function(n, context = c("scan", "figure")) {
  context <- match.arg(context)
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
      n %% 1 != 0 || n < 1) {
    stop("parameter scout: the ", context, " count must be a single whole ",
         "number of at least 1, not ", paste(deparse(n), collapse = " "), ".",
         call. = FALSE)
  }
  if (n <= .BH_MAX_BAMS) return(invisible(TRUE))
  detail <- switch(context,
    scan   = paste0("The BAM set has ", n, " files."),
    figure = paste0("The percentile table has ", n, " rows."))
  stop("parameter scout can only inspect up to ", .BH_MAX_BAMS,
       " BAM files at a time. ", detail, "\n",
       "  Past ", .BH_MAX_BAMS, " the annotation block no longer fits the panel ",
       "and the distributions are too compressed to read, so the figure would ",
       "mislead rather than inform.\n",
       "  Scan in subsets of ", .BH_MAX_BAMS, " or fewer, naming each subset ",
       "with bam_txt_list rather than moving files. The percentile table and ",
       "the suggested cutoffs are unaffected by the limit; it is the figure ",
       "that cannot carry more.", call. = FALSE)
}

## --- section: internal helpers, writing ------------------------------------

#' Write one table with LF endings and no scientific notation
#'
#' Internal. \code{open = "wb"} rather than \code{eol = "\\n"}: on Windows a
#' text-mode connection translates the newline regardless. \code{scipen = 999}
#' keeps large values out of scientific notation, and is restored on exit.
#'
#' @noRd
.write_tsv_lf <- function(x, path) {
  old <- getOption("scipen")
  on.exit(options(scipen = old), add = TRUE)
  options(scipen = 999)
  con <- file(path, open = "wb")
  on.exit(close(con), add = TRUE)
  utils::write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE,
                     eol = "\n")
  invisible(path)
}

#' Write plain lines with LF endings
#' @noRd
.write_txt_lf <- function(lines, path) {
  con <- file(path, open = "wb")
  on.exit(close(con), add = TRUE)
  writeLines(lines, con, sep = "\n")
  invisible(path)
}

#' Write a scout result to a directory
#'
#' Internal. Four artefacts: the two tables, the suggested pairs in a form
#' meant to be copied by hand, and a run log. Paths are recorded as given.
#'
#' @param scout The list returned by \code{parameter_scout()}.
#' @param out_dir Destination directory, created if absent.
#' @param bam_location,annotation_file As passed to the scan, for the log.
#' @param elapsed_min Wall-clock minutes, or NA.
#' @param label Short dataset name, recorded in the log.
#' @return Invisibly, a named character vector of the absolute paths written.
#' @noRd
.write_scout_outputs <- function(scout, out_dir, bam_location = NA_character_,
                                 annotation_file = NA_character_,
                                 elapsed_min = NA_real_, label = NA_character_) {
  p  <- scout$parameters
  pr <- scout$pairs

  files <- c(
    percentiles  = .write_tsv_lf(scout$percentiles,
                                 file.path(out_dir, "scout_percentiles.tsv")),
    distribution = .write_tsv_lf(scout$distribution,
                                 file.path(out_dir, "scout_distribution.tsv")),
    ## Both pairs, from scout$pairs, the same object the figure draws from.
    suggestions  = .write_txt_lf(c(
      "Possible Cutoff Parameters",
      "",
      unlist(lapply(seq_len(nrow(pr)), function(i) c(
        sprintf("%s (%s):", pr$pair[i], gsub(" / ", "/", pr$from[i])),
        paste("low_coverage_cutoff  =", pr$low[i]),
        paste("high_coverage_cutoff =", pr$high[i]),
        ""))),
      "You will need to edit feature_file_editor() as Parameter Scout does not set them for you.",
      "",
      "Other coverage cutoffs may suit your data better; however, try these first."
    ), file.path(out_dir, "scout_suggestions.txt")),
    log          = .write_txt_lf(c(
      paste("parameter scout run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      paste("R", getRversion()),
      paste("bam_location    =", bam_location),
      paste("annotation      =", p$annotation),
      paste("annotation_file =", annotation_file),
      paste("bam_files       =", paste(p$bam_files, collapse = ", ")),
      paste("n_bam           =", p$n_bam),
      paste("paired_end      =", p$paired_end_data),
      paste("strandedness    =", p$strandedness),
      paste("coverage_model  =", p$coverage_model),
      paste("mapqFilter      =", p$mapqFilter),
      paste("elapsed_min     =", elapsed_min),
      paste("label           =", label)
    ), file.path(out_dir, "scout_run_log.txt"))
  )
  ## Absolute, so a returned path survives a change of working directory.
  files[] <- normalizePath(files, winslash = "/")
  invisible(files)
}

## --- section: parameter_scout(), exported ----------------------------------

#' Scout coverage cutoffs from intergenic coverage
#'
#' Reports, for each BAM, the distribution of per-base coverage in the
#' intergenic regions, so \code{low_coverage_cutoff} and
#' \code{high_coverage_cutoff} can be chosen from the data, and returns a
#' suggested pair. Advisory: it sets no parameter.
#'
#' Coverage is built by the same machinery \code{peak_union_calc()} uses, so
#' \code{paired_end_data}, \code{strandedness}, \code{scanbamparam},
#' \code{mapqFilter} and \code{coverage_model} must be given the same values
#' here as in \code{feature_file_editor()}, or the cutoffs will not transfer.
#'
#' Percentiles are pooled across strands and exclude zero-coverage positions,
#' so each describes expressed intergenic positions. tRNA and rRNA are masked
#' out of the intergenic space; tmRNA, RNase P RNA and SRP RNA are not.
#'
#' Reads every BAM once, at roughly the cost of one
#' \code{feature_file_editor()} pass.
#'
#' @param bam_location The directory containing BAM files.
#' @param bam_txt_list Optional newline separated text file of BAM filenames.
#' @param annotation_file GFF3 genome annotation file, or a GFF cache.
#' @param annot_file_dir The directory containing the annotation file.
#' @param original_sRNA_annotation Biotype of pre-annotated ncRNA, or "unknown".
#' @param paired_end_data A boolean indicating if the reads are paired-end.
#'   (Default: FALSE, matching \code{feature_file_editor()}.)
#' @param strandedness "stranded" or "reversely_stranded" (default
#'   "stranded"). It describes the library: a dUTP protocol, which most
#'   bacterial RNA-seq kits use, is "reversely_stranded". A wrong value neither
#'   stops the scan nor warns; it assigns every fragment to the opposite
#'   strand, and the percentiles stop describing the intergenic background.
#' @param scanbamparam Optional \code{ScanBamParam}. When \code{NULL}, the
#'   \code{peak_union_calc()} default filter is built from \code{mapqFilter}.
#' @param mapqFilter Integer. Minimum mapping quality, default 10. See
#'   \code{\link{peak_union_calc}} for the aligner-specific guidance.
#' @param coverage_model "fragment" (default) or "footprint". Applies only to
#'   paired-end reads.
#' @param probs Numeric probabilities for the percentile table.
#' @param low_prob,high_prob The percentile probabilities used by the
#'   suggested-cutoff rule. Both must appear in \code{probs}.
#' @param plot Logical. With an \code{out_dir}, also draw the figure into it.
#'   (Default: TRUE)
#' @param label Short dataset name used in the figure filenames. Defaults to
#'   the name of the BAM directory. Two datasets written to one directory
#'   under one label collide.
#' @param out_dir Destination directory, created if absent. Default
#'   \code{"scout"}, holding the two tables, the suggested pairs, the run log
#'   and, unless \code{plot = FALSE}, the PNG and PDF. \code{NULL} returns
#'   everything while writing nothing.
#' @return Invisibly, a list with \code{percentiles}, \code{distribution},
#'   \code{suggested}, \code{pairs}, \code{scanbamparam} (the filter this run
#'   used, to pass to \code{feature_file_editor()}), \code{parameters} and
#'   \code{files} (the absolute paths written).
#' @seealso \code{\link{suggest_cutoffs}}, \code{\link{peak_union_calc}}
#'
#' @import IRanges
#' @import GenomicAlignments
#' @import Rsamtools
#' @import GenomicRanges
#' @importFrom S4Vectors Rle runValue runLength
#' @importFrom Seqinfo seqnames seqlengths seqinfo
#' @export
parameter_scout <- function(bam_location = ".", bam_txt_list = "",
                            annotation_file, annot_file_dir = ".",
                            original_sRNA_annotation,
                            paired_end_data = FALSE,
                            strandedness = "stranded",
                            scanbamparam = NULL, mapqFilter = 10,
                            coverage_model = c("fragment", "footprint"),
                            probs = c(0.25, 0.5, 0.75, 0.8, 0.85, 0.9),
                            low_prob = 0.5, high_prob = 0.75,
                            plot = TRUE,
                            label = basename(normalizePath(bam_location,
                                                           mustWork = FALSE)),
                            out_dir = "scout") {
  coverage_model <- match.arg(coverage_model)
  ## Same guard as peak_union_calc: prediction cutoffs need stranded coverage.
  valid_strandedness <- c("stranded", "reversely_stranded")
  if (!strandedness %in% valid_strandedness) {
    stop("Invalid 'strandedness' value: '", strandedness,
         "'. Must be one of: ", paste(valid_strandedness, collapse = ", "), ".",
         call. = FALSE)
  }
  stopifnot("probs must lie in [0, 1]" = all(probs >= 0 & probs <= 1),
            "low_prob and high_prob must be in probs" =
              all(c(low_prob, high_prob) %in% probs))

  ## An empty label, from a drive root, would give scout_figure_.png.
  if (!nzchar(label)) label <- "scout"

  ## Before the scan: a mistyped out_dir must fail in seconds.
  if (!is.null(out_dir)) {
    stopifnot("`out_dir` must be a single directory path" =
                is.character(out_dir) && length(out_dir) == 1L && nzchar(out_dir))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(out_dir)) {
      stop("parameter_scout: could not create out_dir: ", out_dir, call. = FALSE)
    }
  }
  .t0 <- Sys.time()

  bam_files <- unlist(.list_bam_files(bam_location, bam_txt_list))
  stopifnot("No BAM files found" = length(bam_files) > 0)
  ## Before the loop, for the same reason.
  .check_bam_limit(length(bam_files), "scan")
  gff_cache <- .resolve_gff_cache(annotation_file, annot_file_dir)
  if (is.null(scanbamparam)) {
    scanbamparam <- .default_scanbamparam(paired_end_data, mapqFilter)
  }

  igr <- NULL
  igr_seqinfo <- NULL
  pct_rows  <- list()
  dist_rows <- list()

  for (f in bam_files) {
    sr <- .strand_split_reads(f, paired_end_data, strandedness,
                              scanbamparam, coverage_model)
    ## Built once from the first BAM's header; every later BAM must agree.
    if (is.null(igr)) {
      igr_seqinfo <- sr$seqinfo
      igr <- .igr_regions(gff_cache, original_sRNA_annotation, igr_seqinfo)
    } else if (!identical(seqlengths(sr$seqinfo), seqlengths(igr_seqinfo))) {
      stop("BAM files disagree on the reference sequence (", basename(f),
           " differs from the first BAM).", call. = FALSE)
    }
    if (sr$n_reads == 0L) {
      warning(basename(f), " yielded no reads after filtering; its percentiles ",
              "will be NA.", call. = FALSE, immediate. = TRUE)
    }

    igr_plus  <- .igr_coverage_rle(sr$plus_reads,  igr$plus)
    igr_minus <- .igr_coverage_rle(sr$minus_reads, igr$minus)

    ## Percentiles are computed on both strands pooled.
    pooled <- c(igr_plus, igr_minus)
    q <- .rle_weighted_quantile(pooled, probs)
    pct_rows[[f]] <- data.frame(
      File = basename(f),
      as.list(stats::setNames(as.numeric(q), .pct_colnames(probs))),
      n_positions = attr(q, "n_positions"),
      stringsAsFactors = FALSE, check.names = FALSE
    )
    dist_rows[[f]] <- rbind(
      .rle_summary(igr_plus,  basename(f), "+"),
      .rle_summary(igr_minus, basename(f), "-")
    )
  }

  percentiles <- do.call(rbind, c(pct_rows, list(make.row.names = FALSE)))
  out <- list(
    percentiles  = percentiles,
    distribution = do.call(rbind, c(dist_rows, list(make.row.names = FALSE))),
    suggested    = suggest_cutoffs(percentiles,
                                   low  = .pct_colnames(low_prob),
                                   high = .pct_colnames(high_prob)),
    ## The offered pairs, from the definition the figure also reads.
    pairs        = .cutoff_pairs(percentiles),
    ## The filter this run used, to hand to feature_file_editor(scanbamparam =).
    scanbamparam = scanbamparam,
    parameters   = list(paired_end_data = paired_end_data,
                        strandedness    = strandedness,
                        mapqFilter      = bamMapqFilter(scanbamparam),
                        coverage_model  = coverage_model,
                        probs           = probs,
                        annotation      = gff_cache$path,
                        n_bam           = length(bam_files),
                        ## Basenames rather than a count: a table found later
                        ## with n_bam = 6 cannot be checked against a directory
                        ## that now holds seven. They match the File column of
                        ## `percentiles`, so the two can be joined.
                        bam_files       = basename(bam_files))
  )

  ## An out_dir writes the four text artefacts and, unless plot = FALSE, the
  ## figure. NULL computes and returns while writing nothing.
  out$files <- character(0)
  if (!is.null(out_dir)) {
    out$files <- .write_scout_outputs(
      out, out_dir,
      bam_location    = bam_location,
      annotation_file = if (is.character(annotation_file)) annotation_file
                        else NA_character_,
      elapsed_min     = round(as.numeric(difftime(Sys.time(), .t0,
                                                  units = "mins")), 2),
      label           = label)
    if (isTRUE(plot)) {
      ## After the tables are on disk, so a failed figure cannot cost the scan.
      fig <- plot_scout_distribution(scout = out, out_dir = out_dir,
                                     label = label)
      out$files <- c(out$files, fig$files)
    }
    ## Already absolute: both writers build them that way.
    for (f in out$files) cat("wrote:", f, "\n")
  }
  ## Invisible: the value carries the whole distribution table.
  invisible(out)
}

## --- section: suggest_cutoffs(), exported ----------------------------------

#' Suggest coverage cutoffs from an intergenic percentile table
#'
#' Applies the percentile rule to the table returned by
#' \code{parameter_scout()}. By default the low cutoff is the median across
#' files of the per-file median intergenic coverage and the high cutoff the
#' median across files of the per-file upper quartile; \code{low} and
#' \code{high} name any other columns. A fractional median is taken up to the
#' next integer.
#'
#' Returns numbers to inspect and record; it does not feed them into the
#' pipeline.
#'
#' @param x The list returned by \code{parameter_scout()}, or its
#'   \code{percentiles} dataframe.
#' @param low The percentile column used for \code{low_coverage_cutoff}.
#' @param high The percentile column used for \code{high_coverage_cutoff}.
#' @return A named integer vector, \code{low_coverage_cutoff} and
#'   \code{high_coverage_cutoff}.
#' @seealso \code{\link{parameter_scout}}
#' @export
suggest_cutoffs <- function(x, low = "Median", high = "P75") {
  pct <- if (is.data.frame(x)) x else x$percentiles
  stopifnot("Column not found in percentile table" =
              all(c(low, high) %in% names(pct)))
  ## The unrounded medians, so a pair differing below the integer is caught.
  lo <- stats::median(pct[[low]],  na.rm = TRUE)
  hi <- stats::median(pct[[high]], na.rm = TRUE)
  if (!is.na(lo) && !is.na(hi) && hi <= lo) {
    warning("The suggested high cutoff (", hi, ") is not above the low cutoff (",
            lo, "). Check the percentile table before using these values.",
            call. = FALSE)
  }
  ## One aggregation rule, held in .agg_percentile().
  c(low_coverage_cutoff  = .agg_percentile(pct[[low]]),
    high_coverage_cutoff = .agg_percentile(pct[[high]]))
}

## --- section: plot_scout_distribution(), exported --------------------------

#' Figure of the intergenic coverage distributions, with the offered cutoffs
#'
#' Draws the percentile table with its \code{Median across BAMs} aggregate row
#' above, one violin-and-box row per BAM below, and the two offered cutoff
#' pairs annotated. Returns the two ggplot objects.
#'
#' Supply exactly one of \code{scout} and \code{in_dir}; \code{in_dir} reads
#' the saved TSVs, which reaches the figure without scanning.
#' \code{out_dir = NULL}, the default, returns the objects without writing.
#' Every other setting is a constant at the head of this file.
#'
#' @param scout The list returned by \code{parameter_scout()}.
#' @param in_dir A directory holding \code{scout_percentiles.tsv} and
#'   \code{scout_distribution.tsv}.
#' @param out_dir Where to write the PNG and PDF, or \code{NULL} to write none.
#' @param label Used in the output filenames.
#' @param plot_title Title above the distributions panel, or \code{NULL}.
#' @return Invisibly, a list with \code{table}, \code{plot} and \code{files}.
#'
#' @importFrom ggplot2 aes annotate coord_cartesian element_blank element_text
#' @importFrom ggplot2 expansion geom_label geom_polygon geom_rect
#' @importFrom ggplot2 geom_segment geom_text geom_vline ggplot labs
#' @importFrom ggplot2 scale_fill_manual scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 theme theme_minimal theme_void
#' @export
plot_scout_distribution <- function(scout = NULL,
                                    in_dir = NULL,
                                    out_dir = NULL,
                                    label = "scout",
                                    plot_title = "Intergenic coverage distribution, expressed positions") {
  ## Exactly one of scout, built in memory, and in_dir, a directory of TSVs.
  if (is.null(scout) == is.null(in_dir)) {
    stop("plot_scout_distribution: supply exactly one of `scout` and `in_dir`. ",
         if (is.null(scout)) "Neither was given." else "Both were given.",
         call. = FALSE)
  }

  if (!is.null(scout)) {
    stopifnot("`scout` must be the list returned by the scan" =
                is.list(scout) && all(c("percentiles", "distribution") %in% names(scout)))
    pct  <- scout$percentiles
    dist <- scout$distribution
  } else {
    dist_path <- file.path(in_dir, "scout_distribution.tsv")
    pct_path  <- file.path(in_dir, "scout_percentiles.tsv")
    stopifnot("scout_distribution.tsv not found in in_dir" = file.exists(dist_path))
    stopifnot("scout_percentiles.tsv not found in in_dir"  = file.exists(pct_path))
    dist <- utils::read.delim(dist_path, stringsAsFactors = FALSE,
                              check.names = FALSE)
    pct  <- utils::read.delim(pct_path,  stringsAsFactors = FALSE,
                              check.names = FALSE)
  }

  ## out_dir = NULL builds the figure and writes nothing, so writing is opt-in.
  if (!is.null(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  stopifnot("distribution table lacks the expected columns" =
              all(c("File", "Coverage", "Weight") %in% names(dist)))
  stopifnot("percentile table lacks the expected columns" =
              all(c("File", "P25", "Median", "P75") %in% names(pct)))

  ## When every BAM is emptied by the read filter the distribution table has no
  ## rows and stats::aggregate() below would stop. Return rather than stop: the
  ## percentile table, the suggested pairs and the run log are already written
  ## and correctly report NA. A single dead BAM draws as an empty row.
  if (nrow(dist) == 0L) {
    warning("No coverage to plot: every BAM in this run was left empty by the ",
            "read filter, so no figure was drawn. The percentile table and the ",
            "run log were still written and report NA. Check mapqFilter against ",
            "your aligner, and check that this module matches your library, ",
            "since the paired-end module keeps only properly paired reads and ",
            "discards single-end data entirely.",
            call. = FALSE, immediate. = TRUE)
    return(invisible(list(table = NULL, plot = NULL, files = character(0))))
  }

  ## The same ceiling the scan applies; a caller may arrive with a saved table.
  .check_bam_limit(nrow(pct), "figure")

  cat("=== plot_scout_distribution ===\n")
  cat("Source    :", if (!is.null(scout)) "scout object, in memory" else in_dir, "\n")
  cat("Label     :", label, "\n")
  cat("Files     :", nrow(pct), "\n\n")


  ## Pool the strands, as the cutoff rule does. Per-strand rows stay in the TSV.

  if ("Strand" %in% names(dist)) {
    pooled <- stats::aggregate(list(Weight = dist$Weight),
                               by = list(File = dist$File, Coverage = dist$Coverage),
                               FUN = sum)
  } else {
    pooled <- dist[, c("File", "Coverage", "Weight")]
  }
  pooled <- pooled[order(pooled$File, pooled$Coverage), ]

  ## File order follows the percentile table. The drawn row is flipped below and
  ## the axis labels reversed to match; the two must stay in step.
  file_levels <- pct$File
  pooled <- pooled[pooled$File %in% file_levels, ]
  stopifnot("no rows left after matching the two tables on File" = nrow(pooled) > 0)


  ## Nested rather than hoisted, because nesting is what preserves .gap()'s
  ## <<- against .y_cur: hoisting it would send that superassignment to the
  ## global environment and the annotation block would misposition.

  ## Weighted standard deviation, for the bandwidth rule.
  .wsd <- function(x, w) {
    mu <- sum(x * w) / sum(w)
    sqrt(sum(w * (x - mu)^2) / sum(w))
  }

  ## Silverman's rule on the weighted sample; n is the total weight.
  .bandwidth <- function(x, w, iqr) {
    n <- sum(w)
    s <- .wsd(x, w)
    a <- min(s, iqr / 1.349)
    if (!is.finite(a) || a <= 0) a <- s
    if (!is.finite(a) || a <= 0) a <- 1
    0.9 * a * n^(-0.2)
  }

  ## Tukey whiskers: extreme values inside the 1.5 IQR fences, hinges from pct.
  .whiskers <- function(cov, q1, q3) {
    iqr <- q3 - q1
    lo_fence <- q1 - 1.5 * iqr
    hi_fence <- q3 + 1.5 * iqr
    inside <- cov[cov >= lo_fence & cov <= hi_fence]
    if (!length(inside)) return(c(min(cov), max(cov)))
    c(min(inside), max(inside))
  }


  ## --- per-file geometry ----------------------------------------------------

  x_max <- max(log2(pooled$Coverage))
  grid_n <- 512

  violin_df <- NULL
  box_df    <- NULL
  med_df    <- NULL
  whisk_df  <- NULL

  ## Row 1 is the bottom of the panel, so the first file takes the highest row.
  row_of_file <- function(i) length(file_levels) - i + 1L

  for (i in seq_along(file_levels)) {
    f <- file_levels[i]
    fr <- row_of_file(i)
    d <- pooled[pooled$File == f, ]
    x <- log2(d$Coverage)
    w <- d$Weight

    prow <- pct[pct$File == f, ]
    q1 <- as.numeric(prow$P25)
    q2 <- as.numeric(prow$Median)
    q3 <- as.numeric(prow$P75)

    ## Density, weighted, on the log2 scale.
    bw <- if (is.null(BW_OVERRIDE)) {
      .bandwidth(x, w, iqr = log2(max(q3, 1)) - log2(max(q1, 1)))
    } else {
      BW_OVERRIDE
    }
    dens <- stats::density(x, weights = w / sum(w), bw = bw,
                           from = 0, to = x_max, n = grid_n)

    violin_df <- rbind(violin_df, data.frame(
      File = f, row = fr, x = dens$x, dens = dens$y, stringsAsFactors = FALSE))

    wk <- .whiskers(d$Coverage, q1, q3)

    box_df <- rbind(box_df, data.frame(
      File = f, row = fr,
      xmin = log2(max(q1, 1)), xmax = log2(max(q3, 1)),
      stringsAsFactors = FALSE))
    med_df <- rbind(med_df, data.frame(
      File = f, row = fr, x = log2(max(q2, 1)), stringsAsFactors = FALSE))
    whisk_df <- rbind(whisk_df, data.frame(
      File = f, row = fr,
      xmin = log2(max(wk[1], 1)), xmax = log2(max(wk[2], 1)),
      stringsAsFactors = FALSE))
  }

  ## Scale the densities into the row height.
  if (identical(VIOLIN_SCALE, "max")) {
    peak <- stats::aggregate(list(peak = violin_df$dens),
                             by = list(File = violin_df$File), FUN = max)
    violin_df$peak <- peak$peak[match(violin_df$File, peak$File)]
    violin_df$h <- VIOLIN_H * violin_df$dens / violin_df$peak
  } else {
    violin_df$h <- VIOLIN_H * violin_df$dens / max(violin_df$dens)
  }
  violin_df$ytop <- violin_df$row + violin_df$h

  ## Close each polygon along its baseline so geom_polygon fills correctly.
  poly_df <- do.call(rbind, lapply(file_levels, function(f) {
    v <- violin_df[violin_df$File == f, ]
    v <- v[order(v$x), ]
    data.frame(File = f,
               x = c(v$x, rev(v$x)),
               y = c(v$ytop, rep(v$row[1], nrow(v))),
               stringsAsFactors = FALSE)
  }))
  poly_df$File <- factor(poly_df$File, levels = file_levels)

  ## Okabe-Ito as hex, reordered for fills. RColorBrewer is not a dependency.
  .OKABE_ITO_FILL <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
                       "#56B4E9", "#E69F00", "#999999", "#F0E442")

  .resolve_palette <- function(name, n) {
    if (identical(name, "grey")) {
      return(rep("grey60", n))
    }
    if (identical(name, "Okabe-Ito")) {
      if (n <= length(.OKABE_ITO_FILL)) return(.OKABE_ITO_FILL[seq_len(n)])
      ## Beyond eight the palette recycles rather than interpolating.
      warning("Okabe-Ito provides 8 distinct colours but ", n,
              " files were supplied; colours repeat from the ninth row. ",
              "Every row carries its own axis label, so this is cosmetic, but ",
              "PALETTE = \"grey\" avoids it if the repeat is distracting.",
              call. = FALSE, immediate. = TRUE)
      return(rep_len(.OKABE_ITO_FILL, n))
    }
    if (identical(name, "viridis")) {
      return(grDevices::hcl.colors(n, palette = "viridis"))
    }
    max_n <- c(Dark2 = 8L, Set2 = 8L, Paired = 12L)
    if (requireNamespace("RColorBrewer", quietly = TRUE) && name %in% names(max_n)) {
      cap <- max_n[[name]]
      base_cols <- RColorBrewer::brewer.pal(min(max(n, 3L), cap), name)
      if (n <= cap) return(base_cols[seq_len(n)])
      ## Recycles rather than interpolating, for the reason given above.
      warning("PALETTE '", name, "' provides only ", cap, " distinct colours but ",
              n, " files were supplied; colours repeat from row ", cap + 1, ". ",
              "Every row carries its own axis label, so this is cosmetic, but ",
              "PALETTE = \"grey\" avoids it if the repeat is distracting.",
              call. = FALSE, immediate. = TRUE)
      return(rep_len(base_cols, n))
    }
    fallback <- switch(name, Dark2 = "Dark 3", Set2 = "Set 2", Paired = "Set 3", "Dark 3")
    grDevices::hcl.colors(max(n, 3L), palette = fallback)[seq_len(n)]
  }

  ## ColorBrewer flags its qualitative palettes safe only at three classes.
  if (PALETTE %in% c("Dark2", "Set2", "Paired") && length(file_levels) > 3) {
    message("Note: ColorBrewer flags '", PALETTE, "' colourblind-safe only at ",
            "n = 3; this figure has ", length(file_levels),
            " files. PALETTE = \"Okabe-Ito\" is safe at 8.")
  }

  pal <- .resolve_palette(PALETTE, length(file_levels))
  names(pal) <- file_levels


  ## --- the plot -------------------------------------------------------------

  box_y <- function(row) row - BOX_GAP

  p_plot <- ggplot() +
    geom_polygon(data = poly_df,
                 aes(x = x, y = y, group = File, fill = File),
                 colour = NA, alpha = FILL_ALPHA) +
    geom_segment(data = whisk_df,
                 aes(x = xmin, xend = xmax,
                     y = box_y(row), yend = box_y(row)),
                 linewidth = 0.4, colour = "grey25") +
    geom_rect(data = box_df,
              aes(xmin = xmin, xmax = xmax,
                  ymin = box_y(row) - BOX_H / 2, ymax = box_y(row) + BOX_H / 2,
                  fill = File),
              colour = "grey25", linewidth = 0.4,
              alpha = min(FILL_ALPHA + 0.15, 1)) +
    geom_segment(data = med_df,
                 aes(x = x, xend = x,
                     y = box_y(row) - BOX_H / 2, yend = box_y(row) + BOX_H / 2),
                 linewidth = 0.7, colour = "grey15") +
    scale_fill_manual(values = pal, guide = "none") +
    ## Expansion is zero because coord_cartesian() pins the y view below. The
    ## labels are reversed to match row_of_file(); the two must stay in step.
    scale_y_continuous(breaks = seq_along(file_levels), labels = rev(file_levels),
                       expand = expansion(mult = c(0, 0))) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.04))) +
    labs(x = expression(log[2] * "(coverage)"), y = NULL, title = plot_title) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor   = element_blank(),
          axis.text.y        = element_text(hjust = 1),
          plot.title         = element_text(size = 11, face = "plain"))

  ## --- the offered values, from the percentile table ------------------------

  ## Percentile columns present, for the console echo below.
  pct_cols_preview <- intersect(c("P25", "Median", "P75", "P80", "P85", "P90"),
                                names(pct))

  lad_p25 <- .agg_percentile(pct$P25)
  lad_p50 <- .agg_percentile(pct$Median)
  lad_p75 <- .agg_percentile(pct$P75)
  lad_p80 <- .agg_percentile(pct$P80)

  ## Aggregates by percentile name, so a pair is looked up not hardcoded.
  agg_by_name <- c(P25 = lad_p25, Median = lad_p50, P75 = lad_p75, P80 = lad_p80)

  ## From .cutoff_pairs(), the definition the writers also read.
  pairs_df <- .cutoff_pairs(pct)

  ## Warn if the stringent pair is not at or above the relaxed one.
  if (!any(is.na(c(pairs_df$low, pairs_df$high)))) {
    if (!(pairs_df$low[2] >= pairs_df$low[1] && pairs_df$high[2] >= pairs_df$high[1])) {
      warning("The pair labelled '", pairs_df$pair[2], "' does not have both cutoffs ",
              "at or above the pair labelled '", pairs_df$pair[1], "'. Check ",
              "PAIR_LOW and PAIR_HIGH: the labels assume the second is the ",
              "stricter of the two.", call. = FALSE, immediate. = TRUE)
    }
  }
  ## Each pair must also be internally valid, high strictly above low.
  for (i in seq_len(nrow(pairs_df))) {
    if (!is.na(pairs_df$low[i]) && !is.na(pairs_df$high[i]) &&
        pairs_df$high[i] <= pairs_df$low[i]) {
      warning("The '", pairs_df$pair[i], "' pair has a high cutoff (",
              pairs_df$high[i], ") that is not above its low cutoff (",
              pairs_df$low[i], "). suggest_cutoffs() warns on this too.",
              call. = FALSE, immediate. = TRUE)
    }
  }

  cat("--- Median across BAMs (the table's aggregate row) ---\n")
  print(data.frame(as.list(stats::setNames(
          vapply(pct_cols_preview, function(nm) .agg_percentile(pct[[nm]]), numeric(1)),
          pct_cols_preview)), check.names = FALSE),
        row.names = FALSE)
  cat("\n--- Possible parameter cutoffs ---\n")
  print(pairs_df, row.names = FALSE)

  ## An aggregate of 0 or 1 imposes no floor; the message below reports it.
  cat("\n")

  ## One marker per percentile in MARKER_PCTS. A value of 1 is not drawn:
  ## log2(1) sits on the panel edge. Tags carry a white backing.
  y_top <- length(file_levels) + VIOLIN_H
  ## The values the pairs use, deduplicated: P75 appears in both.
  offered_vals <- agg_by_name[unique(c(PAIR_LOW, PAIR_HIGH))]

  ## Markers come from all computed percentiles, not only the offered ones.
  stopifnot("MARKER_PCTS must name computed percentiles" =
              all(MARKER_PCTS %in% names(agg_by_name)))

  line_df <- data.frame(
    value = unname(agg_by_name[MARKER_PCTS]),
    tag   = MARKER_PCTS,
    stringsAsFactors = FALSE)
  ## NA-safe: an NA index would give a phantom marker that ggplot2 drops.
  line_df <- line_df[!is.na(line_df$value) & line_df$value > 1, , drop = FALSE]
  line_df$x  <- log2(line_df$value)
  line_df$dy <- (seq_len(nrow(line_df)) - 1) * 0.17
  lab_base   <- y_top + 0.38

  ## Report every offered value with no marker beside it, and why.
  .why_unmarked <- function(nm, v) {
    if (is.na(v))            return("percentile is NA (a BAM yielded no reads)")
    if (v == 0)              return("log2(0) is undefined; use the Median for the low cutoff")
    if (v == 1)              return("log2(1) is the panel edge; a cutoff of 1 keeps everything shown")
    if (!nm %in% MARKER_PCTS) return("not in MARKER_PCTS")
    "no line drawn"
  }
  .unmarked <- setdiff(names(offered_vals), line_df$tag)
  if (length(.unmarked) > 0) {
    cat("NOTE: values offered without a marker on the plot:\n")
    for (nm in .unmarked) {
      cat(sprintf("        %-7s = %-5s %s\n", nm, offered_vals[[nm]],
                  .why_unmarked(nm, offered_vals[[nm]])))
    }
  }
  cat("\n")

  ## --- the annotation block -------------------------------------------------
  ##
  ## Every line is its own layer, and plotmath when BOLD_AGG is TRUE: ggplot2
  ## applies one fontface per text layer and plotmath has no line break.

  ## plotmath helpers: quote a literal run, bold() a bold one, paste() them.
  .q  <- function(s) sprintf('"%s"', s)
  .bd <- function(s) sprintf('bold("%s")', s)
  .pm <- function(...) paste0("paste(", paste(c(...), collapse = ", "), ")")

  ## One wrapped sentence as plotmath, the leading label bold on line one only.
  .block_lines <- function(sentence, label, width) {
    ## strwrap() wraps with the indent; the leading spaces are then replaced by
    ## an x offset, plotmath dropping them. width = NULL gives one line.
    raw <- if (is.null(width)) sentence else
             sub("^ +", "", strwrap(sentence, width = width, exdent = ANN_EXDENT))
    if (isTRUE(BOLD_AGG)) {
      rest <- sub(label, "", raw[1], fixed = TRUE)
      out  <- c(.pm(.bd(label), .q(rest)),
                if (length(raw) > 1) vapply(raw[-1], function(l) .pm(.q(l)),
                                            character(1)))
    } else {
      out <- raw
    }
    list(text = unname(out), plain = raw)
  }

  ## One pair line: "Relaxed: Low <low>, High <high> (Median / P75)".
  .pair_line <- function(label, low, high, from) {
    plain <- sprintf("%s: Low %s, High %s (%s)", label, low, high, from)
    if (!isTRUE(BOLD_AGG)) return(list(text = plain, plain = plain))
    list(text = .pm(.bd(paste0(label, ":")),
                    .q(" Low "), .bd(low),
                    .q(", High "), .bd(high),
                    .q(sprintf(" (%s)", from))),
         plain = plain)
  }

  ## Spacing is a fraction of the plotted y range, so constant on any dataset.
  y_span <- (length(file_levels) + VIOLIN_H) - (1 - BOX_GAP - BOX_H / 2)
  line_h <- ANN_LINE_FRAC * y_span
  gap_h  <- ANN_GAP_FRAC  * y_span

  ## Assemble the block, top down, accumulating one row per drawn line.
  ann <- data.frame(text = character(0), plain = character(0), size = numeric(0),
                    indent = numeric(0), y = numeric(0), stringsAsFactors = FALSE)
  ## lab_base is untouched by the offset, so the marker tags do not move.
  .y_cur <- lab_base + ANN_TOP_OFFSET * line_h
  .push <- function(blk, size, exdent_after_first = TRUE) {
    n <- length(blk$text)
    ys <- .y_cur - (seq_len(n) - 1) * line_h
    ann <<- rbind(ann, data.frame(
      text   = unname(blk$text),
      plain  = unname(blk$plain),
      size   = size,
      indent = c(0, rep(if (exdent_after_first) ANN_EXDENT else 0,
                        max(0, n - 1))),
      y      = ys,
      stringsAsFactors = FALSE))
    .y_cur <<- .y_cur - n * line_h
  }
  .gap <- function(n = 1) .y_cur <<- .y_cur - n * line_h

  lead <- .block_lines(ANN_LEAD, "Possible Cutoff Parameters:", WRAP_CHARS)
  .push(lead, ANN_SIZE)

  .y_cur <- .y_cur - gap_h
  .push(.pair_line(pairs_df$pair[1], pairs_df$low[1], pairs_df$high[1],
                   pairs_df$from[1]), ANN_SIZE)
  .y_cur <- .y_cur - gap_h
  .push(.pair_line(pairs_df$pair[2], pairs_df$low[2], pairs_df$high[2],
                   pairs_df$from[2]), ANN_SIZE)

  ## The note, below the second pair, when SHOW_ANN_NOTE is TRUE.
  if (isTRUE(SHOW_ANN_NOTE) && nzchar(ANN_NOTE)) {
    .gap(1)
    .y_cur <- .y_cur - gap_h
    if (is.null(NOTE_WRAP) && nchar(ANN_NOTE) > NOTE_MAX_CHARS) {
      warning("The note is ", nchar(ANN_NOTE), " characters and wrapping is off ",
              "(NOTE_WRAP is NULL); it may run past the panel. Shorten it or set ",
              "NOTE_WRAP to a number.", call. = FALSE, immediate. = TRUE)
    }
    .push(.block_lines(ANN_NOTE, "Note:", NOTE_WRAP), ANN_NOTE_SIZE)
  }

  ## Anchor: the widest line ends at ANN_RIGHT_FRAC unless ANN_X_FRAC overrides.
  .line_frac <- (nchar(ann$plain) + ann$indent) * ANN_CHAR_FRAC *
                  (ann$size / ANN_SIZE)
  ann_x_frac <- if (is.null(ANN_X_FRAC)) {
    max(0, ANN_RIGHT_FRAC - max(.line_frac))
  } else {
    ANN_X_FRAC
  }
  ann_x  <- ann_x_frac * x_max
  ann$x  <- ann_x + ann$indent * ANN_CHAR_FRAC * x_max

  cat("Annotation: ", nrow(ann), " lines, anchored at ",
      sprintf("%.2f", ann_x_frac), " of the x range",
      if (is.null(ANN_X_FRAC)) " (computed)" else " (ANN_X_FRAC)",
      ", lifted ", sprintf("%.2f", ANN_TOP_OFFSET), " line heights (",
      sprintf("%.2f", ANN_TOP_OFFSET * line_h), " data units)\n", sep = "")
  ## Printed even when overridden, so a large gap between the two is visible.
  if (!is.null(ANN_X_FRAC)) {
    cat("            computed anchor would have been ",
        sprintf("%.2f", max(0, ANN_RIGHT_FRAC - max(.line_frac))),
        ", widest line is ", nchar(ann$plain[which.max(.line_frac)]),
        " characters\n", sep = "")
  }

  p_plot <- p_plot +
    geom_vline(data = line_df, aes(xintercept = x),
               linetype = CUT_LINETYPE, linewidth = 0.5, colour = "grey20") +
    geom_label(data = line_df,
               aes(x = x, y = lab_base - dy, label = tag),
               hjust = -0.12, vjust = 1, size = 2.9, colour = "grey35",
               ## linewidth = 0, not the deprecated label.size = 0.
               fill = "white", linewidth = 0,
               label.padding = grid::unit(0.08, "lines")) +
    annotate("text", x = ann$x, y = ann$y, label = ann$text,
             parse = isTRUE(BOLD_AGG),
             hjust = 0, vjust = 1, size = ann$size, colour = "grey20")


  ## Pin the y view so the annotation cannot resize the distributions. Limits
  ## come from the data; coord_cartesian() limits the view, dropping no rows.
  y_floor  <- 1 - BOX_GAP - BOX_H / 2
  y_ceil   <- lab_base
  y_extent <- y_ceil - y_floor
  y_view   <- c(y_floor - 0.06 * y_extent,
                y_ceil  + ANN_HEADROOM_MULT * y_extent)

  ## Does the block fit in the headroom? The first line is drawn with vjust = 1.
  ann_top <- max(ann$y)
  if (ann_top > y_view[2]) {
    over_lines <- (ann_top - y_view[2]) / line_h
    warning("The annotation block overflows the panel by ",
            sprintf("%.2f", over_lines), " line heights and will be clipped. ",
            "Either lower ANN_TOP_OFFSET to ",
            sprintf("%.2f", ANN_TOP_OFFSET - over_lines),
            " or raise ANN_HEADROOM_MULT to ",
            sprintf("%.3f", (ann_top - y_ceil) / y_extent),
            ". Raising the headroom makes the distributions shorter in ",
            "proportion, so the first is the cheaper remedy.",
            call. = FALSE, immediate. = TRUE)
  }

  p_plot <- p_plot + coord_cartesian(ylim = y_view)

  cat("Panel      : y view ", sprintf("%.2f", y_view[1]), " to ",
      sprintf("%.2f", y_view[2]), ", headroom ",
      sprintf("%.2f", (y_view[2] - y_ceil) / line_h),
      " line heights above the marker tags, block uses ",
      sprintf("%.2f", (ann_top - y_ceil) / line_h), "\n", sep = "")


  ## --- the table, drawn as a ggplot so no extra package is needed ------------

  pct_cols <- intersect(c("P25", "Median", "P75", "P80", "P85", "P90"), names(pct))
  tab <- pct[, c("File", pct_cols), drop = FALSE]
  n_bam <- nrow(tab)

  ## Appended as an ordinary row so fmt() formats it with its column.
  agg_row <- data.frame(File = AGG_ROW_LABEL, stringsAsFactors = FALSE)
  for (nm in pct_cols) agg_row[[nm]] <- .agg_percentile(tab[[nm]])
  tab <- rbind(tab, agg_row)

  tab_cols <- names(tab)

  fmt <- function(v, nm) {
    if (nm == "File") return(as.character(v))
    format(as.numeric(v), trim = TRUE)
  }

  ## The aggregate sits AGG_ROW_GAP lower so its rule reads as a rule.
  row_of <- c(seq_len(n_bam) + 1L, n_bam + 2L + AGG_ROW_GAP)
  is_agg <- c(rep(FALSE, n_bam), TRUE)

  cells <- NULL
  for (j in seq_along(tab_cols)) {
    nm <- tab_cols[j]
    cells <- rbind(cells, data.frame(
      col = j, row = row_of,
      label = fmt(tab[[nm]], nm),
      header = FALSE, agg = is_agg, stringsAsFactors = FALSE))
    cells <- rbind(cells, data.frame(
      col = j, row = 1L, label = nm,
      header = TRUE, agg = FALSE, stringsAsFactors = FALSE))
  }
  ## Left-align the File column, right-align everything else.
  cells$hjust <- ifelse(cells$col == 1, 0, 1)
  cells$x <- ifelse(cells$col == 1, 0,
                    (cells$col - 1) / (length(tab_cols) - 1) * 0.72 + 0.28)

  ## The rule between the last BAM and the aggregate, midway across the gap.
  agg_rule_y <- -(n_bam + 1.5 + AGG_ROW_GAP / 2)

  p_table <- ggplot(cells, aes(x = x, y = -row)) +
    geom_text(aes(label = label, hjust = hjust,
                  fontface = ifelse(header | agg, "bold", "plain")),
              size = 3.2) +
    annotate("segment", x = -0.02, xend = 1.02,
             y = -1.5, yend = -1.5, linewidth = 0.4, colour = "grey40") +
    annotate("segment", x = -0.02, xend = 1.02,
             y = agg_rule_y, yend = agg_rule_y,
             linewidth = AGG_RULE_WIDTH, colour = AGG_RULE_COL) +
    scale_x_continuous(limits = c(-0.03, 1.06), expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0.12, 0.12))) +
    theme_void()


  ## --- combine and write ----------------------------------------------------

  draw_figure <- function() {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(
      2, 1, heights = grid::unit(c(TABLE_FRAC, 1 - TABLE_FRAC), "npc"))))
    print(p_table, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(p_plot,  vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
    grid::popViewport()
  }

  files <- character(0)
  if (is.null(out_dir)) {
    ## Say so rather than falling silent: silence here reads as a failure.
    cat("No out_dir given, so nothing was written.",
        "The ggplot objects are returned.\n")
  } else {
    fig_stem <- label

    png_path <- file.path(out_dir, sprintf("scout_figure_%s.png", fig_stem))
    grDevices::png(png_path, width = FIG_WIDTH, height = FIG_HEIGHT,
                   units = "in", res = FIG_RES)
    draw_figure()
    grDevices::dev.off()
    cat("wrote:", normalizePath(png_path, winslash = "/"), "\n")

    pdf_path <- file.path(out_dir, sprintf("scout_figure_%s.pdf", fig_stem))
    grDevices::pdf(pdf_path, width = FIG_WIDTH, height = FIG_HEIGHT)
    draw_figure()
    grDevices::dev.off()
    cat("wrote:", normalizePath(pdf_path, winslash = "/"), "\n")

    ## Absolute, matching the two lines printed above.
    files <- c(png = png_path, pdf = pdf_path)
    files[] <- normalizePath(files, winslash = "/")
  }

  cat("\nBandwidth rule :", if (is.null(BW_OVERRIDE)) "Silverman, weighted" else BW_OVERRIDE, "\n")
  cat("Violin scaling :", VIOLIN_SCALE, "\n")
  cat("Palette        :", PALETTE, "at alpha", FILL_ALPHA, "\n")
  cat("Pairs offered  :",
      paste(sprintf("%s %s/%s", pairs_df$pair, pairs_df$low, pairs_df$high),
            collapse = " | "), "\n")
  cat("Bold numbers   :", if (isTRUE(BOLD_AGG)) "yes, via plotmath" else "no", "\n")
  cat("Markers drawn  :", paste(line_df$tag, collapse = ", "), "\n")
  cat("Note           :", if (isTRUE(SHOW_ANN_NOTE) && nzchar(ANN_NOTE))
                          "shown" else "none", "\n")
  cat("=== done ===\n")

  ## Returned so a caller can inspect or re-render without scanning again.
  invisible(list(table = p_table, plot = p_plot, files = files))
}
