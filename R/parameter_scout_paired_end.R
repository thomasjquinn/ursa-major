# =============================================================================
# parameter_scout_paired_end.R
# Version 1.0
#
# Optional parameter-scouting module for paired-end data. Advisory only.
#
# Reads BAMs and an annotation, reports the coverage percentiles of the
# intergenic regions, and draws one figure of those distributions, so that
# low_coverage_cutoff and high_coverage_cutoff can be chosen from the data
# rather than guessed. It writes no annotation and sets no parameter. Nothing
# else in baerhunter calls it, and deleting it leaves the package unchanged.
#
# Three exported functions:
#
#   parameter_scout_paired_end()  the scan, and the one most callers want.
#                                 Minutes, scaling with depth and BAM count.
#   plot_scout_distribution()     the figure alone, from the scan's object or
#                                 from saved tables. Seconds.
#   suggest_cutoffs()             the cutoff rule alone, on a percentile table.
#
# The sibling parameter_scout_single_end.R differs only in its defaults and in
# this header, not in its logic: every branch on library type is inside
# .strand_split_reads(), which takes paired_end_data as an argument. The two
# export DIFFERENT names because a package collates all of R/ into one
# namespace, and two files defining one name leave only the later-collated one
# standing, silently.
#
# This file carries its own copies of .list_bam_files(),
# .default_scanbamparam() and .strand_split_reads(), taken from
# peak_union_calc() in feature_file_editor.R. Their bodies must stay identical
# to the code there, or the percentiles stop being on the scale the pipeline
# later applies, which is the whole point of the module. Method after
# Sivasankaran (2024).
#
# What this file needs. It ATTACHES the packages listed below and installs
# none, so a missing one stops the source() line and is named there rather
# than surfacing mid-scan as a missing function. Install the Bioconductor
# ones with BiocManager::install(), ggplot2 with install.packages(). It also
# needs two baerhunter files sourced: count_features.R for
# .resolve_gff_cache() and feature_file_editor.R for major_features(). Being
# files rather than packages these are not attached, so omitting one stops
# the scan at could not find function.
#
# Rle and runValue come from S4Vectors; IRanges exports runLength but neither
# of those, which is why the importFrom is needed. gaps is qualified
# GenomicRanges::gaps, being defined in two packages.
#
# The percentile rule is a convention, not a derived noise floor: the intergenic
# coverage distribution has no antimode separating a noise population from a
# signal one, so the median and P75 were chosen because they give sensible
# output. Anyone replacing it with a data-driven estimate should evaluate noisyR
# first (Moutsopoulos et al. 2021, NAR 49(14), e83), which works on per-base
# coverage and so applies directly here.
# =============================================================================


## --- section: imports, the whole stack -------------------------------------

## ATTACHED RATHER THAN IMPORTED. There is no NAMESPACE while this file is
## sourced, so every Bioconductor generic called unqualified below resolves only
## against the search path, and the figure body makes many unqualified ggplot2
## calls. All EIGHT of these become @importFrom directives, with matching
## DESCRIPTION Imports entries, once the package is built: a package must not
## attach packages at top level.
##
## WHAT THIS DOES AND DOES NOT DO. It does not help a caller who has not
## INSTALLED these: library() fails either way. It helps one who has them and has
## not attached them, and it moves that failure from the middle of a scan, where
## it surfaces as `could not find function "scanBamFlag"` and reads as a defect
## here, to this line, where it names the package that is missing.
##
## THE ORDER IS DELIBERATE. Attaching in a different order changes which package
## masks which, and every percentile this project has recorded was produced under
## this one.
suppressPackageStartupMessages({
  library(S4Vectors)
  library(IRanges)
  library(GenomicRanges)
  library(GenomicAlignments)
  library(Rsamtools)
  ## The seqinfo generics moved to their own package. Older Bioconductor carries
  ## them in GenomeInfoDb, and choosing by what is installed rather than by a
  ## version number keeps this working in both directions.
  if (requireNamespace("Seqinfo", quietly = TRUE)) {
    library(Seqinfo)
  } else {
    library(GenomeInfoDb)
  }
  library(ggplot2)
  library(grid)
})


## --- section: internal constants -------------------------------------------

#' The maximum number of BAM files one scan can carry
#'
#' Internal. Nine. The reasoning is with \code{.check_bam_limit()} below rather
#' than repeated here.
#'
#' Not an argument, deliberately: the remedy in the message, scanning in subsets,
#' gives a caller the numbers without an escape hatch. Relaxing this later is
#' easy; withdrawing an override would not be.
#'
#' @keywords internal
.BH_MAX_BAMS <- 9L


## ---------------------------------------------------------------------------
## Figure configuration. plot_scout_distribution() takes out_dir, label and
## plot_title, plus exactly one of scout and in_dir. Everything else is a plain
## constant below, so the figure is changed by editing a named value.
## ---------------------------------------------------------------------------

## The two cutoffs do not act on one scale. Raising the LOW cutoff moves the
## slice boundaries, so features start and stop in different places and one long
## slice can break into several, each of which may pass: the count can rise while
## every feature shortens. Raising the HIGH cutoff leaves the boundaries alone
## and only changes how much evidence a feature needs, so its output is a strict
## subset of the lower setting's. Only the second is monotone. A reader who takes
## the stringent pair to mean fewer predictions will be wrong.

## Fill palette. Colour is decorative: every row carries its own axis label, so
## the fill encodes nothing. The files are unordered categories, so the palette
## must be qualitative, distinct hues at similar lightness.
##
##   "Okabe-Ito" 8 hues, colour-vision safe at all 8, base R. Default.
##   "Dark2"     ColorBrewer, saturated. Flagged safe only at n = 3.
##   "Set2"      ColorBrewer, pastel. Same n = 3 caveat; washes out at low alpha.
##   "Paired"    ColorBrewer, 12 colours, encodes pairs. Not wanted here.
##   "viridis"   sequential ramp sampled discretely. Varying lightness makes
##               some rows read as heavier.
##   "grey"      one neutral fill, safest in greyscale.
PALETTE <- "Okabe-Ito"

## Fill transparency.
FILL_ALPHA <- 0.75

## Violin scaling. "area" gives every violin the same area and preserves the
## relative peak heights, which is what ggplot2's geom_violin does by default.
## "max" scales each violin to the same height instead.
VIOLIN_SCALE <- "area"

## Kernel bandwidth on the log2 scale. NULL uses Silverman's rule computed from
## the weighted distribution, which is narrow at these position counts and
## reproduces the lumpy low end. Set a number (try 0.15 to 0.3) to smooth.
BW_OVERRIDE <- NULL

## Device settings.
FIG_WIDTH   <- 9      # inches
FIG_HEIGHT  <- 7      # inches
FIG_RES     <- 300    # dpi, PNG only
TABLE_FRAC  <- 0.32   # fraction of figure height given to the table

## The aggregate row at the foot of the percentile table: one value per column,
## from .agg_percentile(). It is a median per COLUMN, so under Median it is the
## median of the per-BAM medians and under P75 the median of the per-BAM upper
## quartiles.
AGG_ROW_LABEL  <- "Median across BAMs"
AGG_ROW_GAP    <- 0.35     # row-heights of clear space above the aggregate row
AGG_RULE_COL   <- "grey20" # darker than the header rule's grey40
AGG_RULE_WIDTH <- 0.7

## --- the annotation block on the right of the plot -------------------------

## Bold the pair labels and their values. TRUE renders the WHOLE annotation
## block with plotmath, which is the only way to bold part of a text layer in
## ggplot2; FALSE draws every line as plain text.
BOLD_AGG <- TRUE

## The lead sentence and the width it wraps to. Wrapping is computed, so the
## sentence can be edited without re-breaking it by hand. If the label is
## lengthened, check the block still breaks into two lines: its height decides
## whether it clears the distribution tails below it.
ANN_LEAD <- paste("Possible Cutoff Parameters: Use the shape of the",
                  "distributions to inform your choice.")
WRAP_CHARS <- 52

## Text size for the block, and one size down for the notes so they read as a
## footnote to the numbers rather than as a third instruction.
ANN_SIZE      <- 3.1
ANN_NOTE_SIZE <- 2.7

## Hanging indent, in characters, for wrapped continuation lines. Implemented as
## an x offset rather than leading spaces, because plotmath drops leading
## whitespace on some devices. 0 because no line in the block needs one; set 4
## if a future line must hang under a bold label.
ANN_EXDENT <- 0

## Note wrapping. NULL does not wrap, which suits single-line notes; set a
## number for a longer one. A check below warns if a note may run off the panel.
NOTE_WRAP <- NULL
NOTE_MAX_CHARS <- 70

## Horizontal anchor of the annotation block, as a fraction of the x range. All
## elements are left-aligned from it.
##
## NULL computes it, placing the widest line to end at ANN_RIGHT_FRAC. The
## estimator assumes a uniform character width and over-estimates a mostly
## lower-case sentence, which pulls the block left into the distributions, so a
## number is given here instead.
ANN_X_FRAC    <- 0.62
ANN_RIGHT_FRAC <- 0.98

## Vertical start of the annotation block, in LINE HEIGHTS above its default.
## 0 is unlifted. Lifting moves the block beside a shorter tail, where there is
## more clear space to its right.
##
## The unit is a line height and not a data unit, because the panel is a fixed
## physical height whatever the row count, so the same data offset is a
## different distance on a three-row figure than on a six-row one.
##
## It does not move the marker tags, which would push them into the title.
ANN_TOP_OFFSET <- 5.4

## Headroom above the marker tags, as a fraction of the plotted range.
##
## The annotation is an ordinary layer, so ggplot2 would otherwise fit the y
## scale to it and draw every distribution shorter to make room. The y view is
## pinned by coord_cartesian() to limits computed from the data alone, and the
## block is drawn into this headroom instead.
##
## Raising this buys room for a taller block at a proportional cost in
## distribution height. If the block still does not fit, a warning names both
## remedies rather than clipping the text.
ANN_HEADROOM_MULT <- 0.26

## Width of one character at ANN_SIZE, as a fraction of panel width. Used only
## to place the block, never to lay out the text, so an error here shifts the
## block slightly and breaks nothing.
ANN_CHAR_FRAC <- 0.0091

## Vertical spacing, as a fraction of the plotted y range rather than in data
## units, because the panel is a fixed physical height whatever the number of
## rows. A fraction is therefore a constant physical distance on any dataset.
ANN_LINE_FRAC <- 0.050    # between consecutive lines
ANN_GAP_FRAC  <- 0.027    # extra space above each pair line and the note

## Which percentiles get a dashed marker. This is NOT the offered set: P25 is
## marked and not offered, P80 is offered and not marked.
##
## A marker shows where a landmark sits in the distribution; the aggregate row
## above the plot is what identifies the numbers, so an offered value without a
## marker is still findable. P80 is unmarked because its tag collides with P75,
## which sits about 0.47 log2 units away. Add "P80" to restore it.
MARKER_PCTS <- c("P25", "Median", "P75")

## The note below the two pairs. Its wording does not vary with the data.
## "Try these first" rather than "these are the best", since no
## percentile here is validated against prediction output. "Cutoffs" rather than
## "parameters", to avoid colliding with the bold label in the lead sentence.
ANN_NOTE      <- "Note: Try these first. Other cutoffs may suit your data better."
SHOW_ANN_NOTE <- TRUE

## THE TWO PAIRS, following Sivasankaran (2024):
##
##   Relaxed     low = Median   high = P75
##   Stringent   low = P75      high = P80
##
## P75 does both jobs, as the high of the relaxed pair and the low of the
## stringent one, so choosing P75 for both is not a valid pair.
##
## Values are looked up by percentile name from the aggregate row, so the labels
## and the numbers cannot drift apart.
PAIR_LABELS <- c(relaxed = "Relaxed", stringent = "Stringent")
PAIR_LOW    <- c(relaxed = "Median",  stringent = "P75")
PAIR_HIGH   <- c(relaxed = "P75",     stringent = "P80")

## Line style for the percentile markers. All markers share one style, so no
## marker looks more authoritative than another. "dotted" is the other sensible
## choice; "longdash" and "twodash" read poorly at print size.
CUT_LINETYPE <- "dashed"

## Geometry of one file's row: violin above the baseline, box below.
VIOLIN_H <- 0.72
BOX_H    <- 0.16
BOX_GAP  <- 0.14


## --- section: internal helpers, BAM reading --------------------------------

## Taken from peak_union_calc() in feature_file_editor.R. Their bodies must stay
## identical to the code there: the module measures coverage on the scale
## peak_union_calc() later slices, and that holds only while the two agree.

#' List the BAM files for a run
#'
#' Internal helper. Resolves the BAM set from an explicit newline-separated list
#' file, or by scanning a directory.
#'
#' @param bam_location The directory containing BAM files.
#' @param bam_txt_list Optional newline separated text file of BAM filenames.
#'
#' @return A character vector of BAM file paths, or a list of them when
#'   \code{bam_txt_list} is given.
#'
#' @keywords internal
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
#' Internal helper. Drops unmapped, QC-failing, secondary and supplementary
#' alignments; paired-end reads must additionally be properly paired with a
#' mapped mate.
#'
#' @param paired_end_data A boolean indicating if the reads are paired-end.
#' @param mapqFilter Integer. Minimum mapping quality. See
#'   \code{\link{peak_union_calc}} for the aligner-specific guidance.
#'
#' @return An object of class \code{ScanBamParam}.
#'
#' @keywords internal
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
#' Internal helper. Reads one BAM under the supplied filter, applies the
#' paired-end coverage model, and splits the reads by strand, inverting for
#' single-end reversely-stranded libraries.
#'
#' Returns reads, not coverage: the caller runs \code{coverage()} itself, which
#' keeps \code{compute_strand_peaks()} inside \code{peak_union_calc()} unchanged.
#'
#' strandMode maps the library type to the pair's strand, matching featureCounts
#' strandSpecific: 1 takes the pair strand from the FIRST mate ("stranded"),
#' 2 from the LAST mate ("reversely_stranded": dUTP, NSR/NNSR, Illumina stranded
#' TruSeq).
#'
#' @param f Path to a BAM file.
#' @param paired_end_data A boolean indicating if the reads are paired-end.
#' @param strandedness "stranded" or "reversely_stranded".
#' @param scanbamparam A \code{ScanBamParam} object.
#' @param coverage_model "fragment" or "footprint".
#'
#' @return A list with \code{plus_reads}, \code{minus_reads}, \code{n_reads}
#'   (alignments retained by the filter, counted before the split) and
#'   \code{seqinfo} (from the BAM header).
#'
#' @keywords internal
.strand_split_reads <- function(f, paired_end_data, strandedness,
                                scanbamparam, coverage_model) {
  if (paired_end_data) {
    strand_mode <- if (strandedness == "reversely_stranded") 2 else 1
    read_pairs <- readGAlignmentPairs(f, strandMode = strand_mode, param = scanbamparam)
    ## coverage_model, paired-end only: "fragment" counts the whole pair from
    ## leftmost to rightmost base, "footprint" only the aligned blocks of each
    ## mate. The GRanges is built here and coverage() run on it downstream,
    ## because coverage() on a GAlignmentPairs uses footprint unconditionally.
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
#' Internal. The single definition of what this module offers. Both the figure
#' and the writers read it, so a text artefact cannot state a different offer
#' from the figure beside it in the same directory.
#'
#' @param percentiles The per-BAM percentile table.
#'
#' @return A data frame of one row per pair, with \code{pair}, \code{low},
#'   \code{high} and \code{from}.
#'
#' @keywords internal
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
#' Internal helper. Builds the per-strand intergenic regions as the complement
#' of the major annotated features. tRNA and rRNA are retained as features by
#' \code{major_features()} and so are masked OUT of the intergenic space, which
#' keeps very high rRNA coverage out of the background. No other RNA class is
#' retained: tmRNA, RNase P RNA and SRP RNA fall INTO the background, because the
#' \code{[^tr]RNA} pattern excludes only a t or an r immediately before "RNA".
#'
#' The reference name and length come from the BAM header, so the regions key
#' against the BAM coverage on any genome, and \code{seqlengths} is set so
#' \code{gaps()} runs to the chromosome end and the terminal region is kept.
#'
#' @param gff_cache A GFF cache (see \code{load_gff_cache}).
#' @param original_sRNA_annotation Biotype of pre-annotated ncRNA, or "unknown".
#' @param bam_seqinfo A \code{Seqinfo} from the BAM header (single sequence).
#'
#' @return A list with slots \code{plus} and \code{minus}, each a GRanges of the
#'   intergenic regions on that strand.
#'
#' @keywords internal
.igr_regions <- function(gff_cache, original_sRNA_annotation, bam_seqinfo) {
  stopifnot("parameter scout expects a single reference sequence" =
              length(seqnames(bam_seqinfo)) == 1)
  seqname   <- as.character(seqnames(bam_seqinfo))
  seqlength <- unname(seqlengths(bam_seqinfo))
  ## A missing seqlength makes gaps() stop at the last feature and silently
  ## drop the terminal intergenic region, so require it from the BAM header.
  stopifnot("BAM header carries no sequence length; cannot define intergenic regions safely" =
              !is.na(seqlength))

  ## target_strand "." returns the full major-feature set, both strands.
  major_f <- major_features(gff_cache, annot_file_directory = ".", ".",
                            original_sRNA_annotation)
  ## GFF3 permits "." in the strand column; keep only stranded features.
  major_f <- major_f[major_f[, 7] %in% c("+", "-"), ]
  stopifnot("No stranded major features found in the annotation" = nrow(major_f) > 0)

  ## Feature POSITIONS come from the GFF and the seqname from the BAM, so the
  ## gaps key against the BAM coverage whatever the GFF seqid string says.
  feat <- GenomicRanges::GRanges(
    seqnames   = seqname,
    ranges     = IRanges(start = as.integer(major_f[, 4]),
                         end   = as.integer(major_f[, 5])),
    strand     = major_f[, 7],
    seqlengths = stats::setNames(seqlength, seqname)
  )
  ## gaps() on a stranded GRanges with known seqlengths returns the per-strand
  ## complement to the chromosome ends, plus a "*"-strand complement spanning
  ## the whole sequence. Taking only "+" and "-" discards the "*" entries.
  g <- GenomicRanges::gaps(feat)
  list(plus  = g[strand(g) == "+"],
       minus = g[strand(g) == "-"])
}


#' Intergenic coverage values as a single Rle
#'
#' Internal helper. Restricts a strand's coverage to the intergenic regions and
#' returns the values as one Rle, without expanding to a per-base vector.
#'
#' @param reads A GAlignments (or GRanges) for one strand.
#' @param igr_granges A GRanges of intergenic regions on that strand.
#'
#' @return An Rle of the coverage values inside \code{igr_granges}.
#'
#' @keywords internal
.igr_coverage_rle <- function(reads, igr_granges) {
  cvg <- coverage(reads)                 # RleList keyed by seqname; BAM extends to seqlength
  if (length(igr_granges) == 0L) return(Rle(integer(0)))
  restricted <- cvg[igr_granges]         # RleList, one element per intergenic range
  unlist(restricted, use.names = FALSE)  # one Rle over all intergenic positions
}


#' Weighted quantiles of a run-length encoded vector
#'
#' Internal helper. Computes quantiles directly from run values and run lengths,
#' reproducing \code{stats::quantile(type = 7)} on the expanded vector without
#' expanding it. Zero-coverage positions are excluded, so a percentile describes
#' the distribution of expressed intergenic positions.
#'
#' @param x An Rle of coverage values.
#' @param probs Numeric vector of probabilities.
#'
#' @return A named numeric vector of quantiles, with attribute
#'   \code{n_positions} (the number of non-zero positions).
#'
#' @keywords internal
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
#' @keywords internal
.pct_colnames <- function(probs) {
  vapply(probs, function(p) {
    if (isTRUE(all.equal(p, 0.5))) "Median" else paste0("P", round(p * 100))
  }, character(1))
}


#' Run-length summary of a coverage Rle, for plotting
#'
#' Internal helper. Collapses a coverage Rle to (value, weight) rows, dropping
#' zeros, which is what a weighted violin or density plot needs without
#' materialising every position.
#'
#' @keywords internal
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
#' Internal helper. The single implementation of the cross-BAM aggregation rule:
#' the median across files of one percentile column, taken up to the next whole
#' number. Both \code{suggest_cutoffs()} and the figure call it, so the pair the
#' module reports and the values the figure draws cannot disagree.
#'
#' \code{ceiling()} rather than \code{round()}, for two reasons. These cutoffs
#' are a noise and expression floor, and rounding a floor up is the conservative
#' direction, trading sensitivity for specificity. And base \code{round()} is
#' round-half-to-even, so a median of 16.5 would round DOWN to 16 while 15.5
#' rounds up; the cutoff would then depend on the parity of the integer below it.
#'
#' \code{na.rm = TRUE} so that one BAM yielding no reads leaves the aggregate
#' usable rather than propagating NA into every offered value.
#'
#' @param v A numeric vector: one percentile column, one value per BAM.
#'
#' @return A single integer.
#'
#' @keywords internal
.agg_percentile <- function(v) {
  as.integer(ceiling(stats::median(as.numeric(v), na.rm = TRUE)))
}

## --- section: .check_bam_limit(), the single BAM-count ceiling -------------

#' Enforce the BAM-count ceiling
#'
#' Internal. The single implementation of the limit, called once before any BAM
#' is opened and once before the figure is drawn, since a caller can reach the
#' figure from a saved percentile table without scanning.
#'
#' Nine, because past it the annotation block overflows the panel and the
#' violins are too compressed to read the shape the figure asks the reader to
#' use. The percentile table and the cutoffs would be sound at any count, and
#' the message says so. A stop rather than a warning: a reader handed a clipped
#' figure has already been misled.
#'
#' @param n Integer. Number of BAM files, or rows in the percentile table.
#' @param context "scan" or "figure", naming where the count came from.
#'
#' @return Invisibly TRUE, or stops.
#'
#' @keywords internal
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
#' Internal. \code{write.table(eol = "\\n")} is NOT sufficient on Windows: a
#' text-mode connection translates the newline regardless, so the file arrives with
#' CRLF. Opening the connection with \code{open = "wb"} is what actually holds.
#'
#' \code{scipen = 999} because a coverage value in the millions would otherwise be
#' written as \code{1e+06} into a table meant to be read by other tools, and the
#' option is restored on exit rather than left set for the session.
#'
#' @keywords internal
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
#' @keywords internal
.write_txt_lf <- function(lines, path) {
  con <- file(path, open = "wb")
  on.exit(close(con), add = TRUE)
  writeLines(lines, con, sep = "\n")
  invisible(path)
}

#' Write a scout result to a directory
#'
#' Internal. Four artefacts: the two tables, the suggested pair in a form meant
#' to be copied by hand, and a run log recording what produced them.
#'
#' Paths are recorded exactly as the caller gave them, since an anchor computed
#' here would be an anchor to somebody else's machine.
#'
#' @param scout The list returned by \code{parameter_scout_paired_end()}.
#' @param out_dir Destination directory, created if absent.
#' @param bam_location,annotation_file As passed to the scan, recorded in the log.
#' @param elapsed_min Wall-clock minutes, or NA.
#' @param label Short dataset name, recorded in the log.
#'
#' @return Invisibly, a named character vector of the ABSOLUTE paths written,
#'   with forward slashes on every platform.
#'
#' @keywords internal
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
    ## Both pairs, from scout$pairs, the same object the figure draws from. The
    ## percentile names sit in the headings, so the shared value is visible.
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
  ## ABSOLUTE, and built that way rather than converted on the way out. A caller
  ## who is given a relative path holds something that stops resolving the moment
  ## they change directory, and out_dir is now a default rather than a value they
  ## typed, so this vector is where they look to find out where the files went.
  ## files[] rather than files preserves the names.
  files[] <- normalizePath(files, winslash = "/")
  invisible(files)
}

## --- section: parameter_scout_paired_end(), exported -----------------------

#' Scout coverage cutoffs from intergenic coverage
#'
#' Reports, for each BAM, the distribution of per-base coverage restricted to the
#' intergenic regions, so that \code{low_coverage_cutoff} and
#' \code{high_coverage_cutoff} can be chosen from the data rather than guessed,
#' and returns a suggested pair. Optional and advisory: it sets no parameter.
#'
#' Coverage is built by the same machinery \code{peak_union_calc()} uses, so the
#' percentiles are on the scale the pipeline later applies. The arguments
#' controlling it (\code{paired_end_data}, \code{strandedness},
#' \code{scanbamparam}, \code{mapqFilter}, \code{coverage_model}) must be given
#' the same values here as in \code{feature_file_editor()}, or the derived
#' cutoffs will not transfer.
#'
#' Percentiles are pooled across strands and zero-coverage positions are
#' excluded, so each describes expressed intergenic positions. tRNA and rRNA are
#' masked out of the intergenic space; tmRNA, RNase P RNA and SRP RNA are not, so
#' their loci remain in the background distribution.
#'
#' Reads every BAM once, at roughly the cost of one \code{feature_file_editor()}
#' pass.
#'
#' @param bam_location The directory containing BAM files.
#' @param bam_txt_list Optional newline separated text file of BAM filenames.
#' @param annotation_file GFF3 genome annotation file, or a GFF cache.
#' @param annot_file_dir The directory containing the annotation file.
#' @param original_sRNA_annotation Biotype of pre-annotated ncRNA, or "unknown".
#' @param paired_end_data A boolean indicating if the reads are paired-end.
#'   (Default: TRUE, this being the paired-end module.)
#' @param strandedness "stranded" or "reversely_stranded". (Default:
#'   "stranded", matching \code{peak_union_calc()}.) It describes the
#'   library rather than a preference: a dUTP protocol, which most
#'   bacterial RNA-seq kits use, is "reversely_stranded". A wrong value
#'   neither stops the scan nor warns. It assigns every fragment to the
#'   opposite strand, so an intergenic region collects coverage from gene
#'   bodies on the other strand and the percentiles stop describing the
#'   intergenic background.
#' @param scanbamparam Optional \code{ScanBamParam}. When \code{NULL}, the same
#'   default filter used by \code{peak_union_calc()} is built from
#'   \code{mapqFilter}.
#' @param mapqFilter Integer. Minimum mapping quality. Default 10. See
#'   \code{\link{peak_union_calc}} for the aligner-specific guidance.
#' @param coverage_model "fragment" (default) or "footprint". It applies only
#'   to paired-end reads, which is what this module reads.
#' @param probs Numeric probabilities for the percentile table.
#' @param low_prob,high_prob The percentile probabilities used by the
#'   suggested-cutoff rule. Both must appear in \code{probs}.
#' @param plot Logical. With an \code{out_dir}, also draw the figure into it.
#'   Ignored when \code{out_dir} is \code{NULL}. (Default: TRUE)
#' @param label Short dataset name used in the figure filenames,
#'   \code{scout_figure_<label>.png} and \code{.pdf}. Defaults to the NAME OF
#'   THE BAM DIRECTORY, so a scan run inside a folder called \code{bovis}
#'   produces \code{scout_figure_bovis.png} with nothing passed. Name the
#'   dataset rather than the tool if you set it: two datasets written to one
#'   directory under one label collide.
#' @param out_dir Destination directory, created if absent, and resolved
#'   against the working directory like any other relative path. The default is
#'   \code{"scout"}, so a scan run inside the folder holding the BAM files and
#'   the annotation leaves a \code{scout} subfolder beside them holding six
#'   files: the two tables, the suggested pair, the run log, and (unless
#'   \code{plot = FALSE}) the PNG and PDF of the figure. Pass \code{NULL} to
#'   compute and return everything while writing nothing.
#'
#' @return Invisibly, a list with \code{percentiles} (one row per BAM), \code{distribution}
#'   (run-length summary for weighted plotting), \code{suggested} (named integer
#'   pair from \code{suggest_cutoffs()}), \code{pairs} (both offered cutoff
#'   pairs, the same table the figure draws), \code{scanbamparam} (the filter
#'   this run used, to pass to \code{feature_file_editor()} so both stages filter
#'   identically), \code{parameters} (the coverage-construction arguments) and
#'   \code{files} (the ABSOLUTE paths written, empty when \code{out_dir} is
#'   \code{NULL}).
#'
#' @seealso \code{\link{suggest_cutoffs}}, \code{\link{peak_union_calc}}
#'
#' @import IRanges
#' @import GenomicAlignments
#' @import Rsamtools
#' @import GenomicRanges
#' @importFrom S4Vectors Rle runValue runLength
#' @importFrom Seqinfo seqnames seqlengths
#' @export
parameter_scout_paired_end <- function(bam_location = ".", bam_txt_list = "",
                            annotation_file, annot_file_dir = ".",
                            original_sRNA_annotation,
                            paired_end_data = TRUE,
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
  ## Same guard as peak_union_calc: a cutoff for prediction must come from
  ## stranded coverage.
  valid_strandedness <- c("stranded", "reversely_stranded")
  if (!strandedness %in% valid_strandedness) {
    stop("Invalid 'strandedness' value: '", strandedness,
         "'. Must be one of: ", paste(valid_strandedness, collapse = ", "), ".",
         call. = FALSE)
  }
  stopifnot("probs must lie in [0, 1]" = all(probs >= 0 & probs <= 1),
            "low_prob and high_prob must be in probs" =
              all(c(low_prob, high_prob) %in% probs))

  ## The default label is the BAM directory's name, which is empty in the one
  ## case where that directory is a drive or filesystem root. An empty label
  ## would produce a figure called scout_figure_.png, so fall back to the tool's
  ## own name. Guarded here rather than in the default because a helper would be
  ## a new top-level definition.
  if (!nzchar(label)) label <- "scout"

  ## Before the scan, not after. A mistyped out_dir must fail in seconds rather
  ## than after minutes of BAM reading with the result then thrown away. Same
  ## principle as the BAM ceiling below.
  if (!is.null(out_dir)) {
    stopifnot("`out_dir` must be a single directory path" =
                is.character(out_dir) && length(out_dir) == 1L && nzchar(out_dir))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(out_dir)) {
      stop("parameter_scout_paired_end: could not create out_dir: ", out_dir, call. = FALSE)
    }
  }
  .t0 <- Sys.time()

  bam_files <- unlist(.list_bam_files(bam_location, bam_txt_list))
  stopifnot("No BAM files found" = length(bam_files) > 0)
  ## Before the loop, not after: a scan of ten BAMs is tens of minutes of reading
  ## that would end in the same refusal.
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
    ## Build the intergenic regions once, from the first BAM's header, then
    ## require every subsequent BAM to agree.
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
    ## The offered pairs, from the single definition the figure also reads, so the
    ## artefacts written below cannot state a different offer from the figure.
    pairs        = .cutoff_pairs(percentiles),
    ## The filter this run used, returned so the caller can hand the SAME object
    ## to feature_file_editor(scanbamparam = ...). Both stages then filter
    ## identically by value rather than by two constructions agreeing.
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

  ## An out_dir, the default being "scout", writes the four text artefacts and,
  ## unless plot = FALSE, draws the figure into the same directory, so one call
  ## yields six files. out_dir = NULL computes and returns while writing
  ## nothing, as plot_scout_distribution() also does, and is the way to ask for
  ## that rather than the way to get it by omission.
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
      ## Rendered AFTER the tables are on disk, so a figure that fails to draw
      ## cannot cost a completed scan. Its own writer reports the two files.
      fig <- plot_scout_distribution(scout = out, out_dir = out_dir,
                                     label = label)
      out$files <- c(out$files, fig$files)
    }
    ## The paths are already absolute: both writers build them that way, so the
    ## console and out$files are the same strings rather than two conversions of
    ## one path. A caller working by hand has the console and nothing else, and a
    ## relative path is ambiguous exactly when the working directory is not the
    ## one they assumed.
    for (f in out$files) cat("wrote:", f, "\n")
  }
  ## Invisibly, because this is called for its artefacts and for the console
  ## summary above. The value carries the whole distribution table, tens of
  ## thousands of rows, so a bare call at the prompt would print it.
  invisible(out)
}

## --- section: suggest_cutoffs(), exported ----------------------------------

#' Suggest coverage cutoffs from an intergenic percentile table
#'
#' Applies the pooled percentile rule to the table returned by
#' \code{parameter_scout_paired_end()} and returns a pair of cutoffs. By default
#' the low cutoff is the median across files of the per-file median intergenic
#' coverage and the high cutoff the median across files of the per-file upper
#' quartile; \code{low} and \code{high} name any other columns.
#' A fractional median is taken up to the next integer, so the cutoff is rounded
#' in the conservative direction (see \code{.agg_percentile}).
#'
#' Returns numbers for the caller to inspect and record. It does not feed them
#' into the pipeline: cutoffs kept as explicit, written-down values are what
#' preserves reproducibility.
#'
#' @param x The list returned by \code{parameter_scout_paired_end()}, or its
#'   \code{percentiles} dataframe.
#' @param low The percentile column used for \code{low_coverage_cutoff}.
#' @param high The percentile column used for \code{high_coverage_cutoff}.
#'
#' @return A named integer vector, \code{low_coverage_cutoff} and
#'   \code{high_coverage_cutoff}.
#'
#' @seealso \code{\link{parameter_scout_paired_end}}
#' @export
suggest_cutoffs <- function(x, low = "Median", high = "P75") {
  pct <- if (is.data.frame(x)) x else x$percentiles
  stopifnot("Column not found in percentile table" =
              all(c(low, high) %in% names(pct)))
  ## The warning compares the UNROUNDED medians, as it always has, so that a
  ## high and low which differ only below the integer are still caught.
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
#' Draws a percentile table with a \code{Median across BAMs} aggregate row above,
#' and one violin-and-box row per BAM below, with the two offered cutoff pairs
#' annotated. Returns the two ggplot objects.
#'
#' Supply exactly one of \code{scout} and \code{in_dir}. \code{in_dir} reads the
#' saved TSVs, which is how a caller reaches the figure without scanning.
#'
#' \code{out_dir = NULL}, the default, builds the figure and returns the objects
#' without writing, so writing is opt-in. Every other setting is an internal
#' constant near the head of this file: editing the source is how the figure is
#' changed.
#'
#' @param scout The list returned by \code{parameter_scout_paired_end()}.
#' @param in_dir A directory holding \code{scout_percentiles.tsv} and
#'   \code{scout_distribution.tsv}.
#' @param out_dir Where to write the PNG and PDF, or \code{NULL} to write none.
#' @param label Used in the output filenames.
#' @param plot_title Title above the distributions panel, or \code{NULL}.
#'
#' @return Invisibly, a list with \code{table}, \code{plot} and \code{files}.
#' @export
plot_scout_distribution <- function(scout = NULL,
                                    in_dir = NULL,
                                    out_dir = NULL,
                                    label = "scout",
                                    plot_title = "Intergenic coverage distribution, expressed positions") {
  ## The tables arrive one of two ways and exactly one must be given. `scout` is
  ## what the scan has just built in memory; `in_dir` is a directory of saved
  ## TSVs, which is how a caller reaches the figure without scanning at all.
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
    dist <- read.delim(dist_path, stringsAsFactors = FALSE, check.names = FALSE)
    pct  <- read.delim(pct_path,  stringsAsFactors = FALSE, check.names = FALSE)
  }

  ## out_dir = NULL means build the figure and return the objects, writing
  ## nothing. Writing is therefore opt-in.
  if (!is.null(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  stopifnot("distribution table lacks the expected columns" =
              all(c("File", "Coverage", "Weight") %in% names(dist)))
  stopifnot("percentile table lacks the expected columns" =
              all(c("File", "P25", "Median", "P75") %in% names(pct)))

  ## NOTHING TO PLOT. When every BAM yields no reads after filtering, every
  ## .rle_summary() returns a zero-row frame and the distribution table is
  ## empty. There is no distribution to draw, and stats::aggregate() below stops
  ## with "no rows to aggregate", which surfaces as a stack trace from inside
  ## stats rather than as anything a caller can act on.
  ##
  ## Return instead of stopping. The scan that called this has already written
  ## the percentile table, the suggested pair and the run log, all of which
  ## correctly say NA, and an error here would lose them and take the whole run
  ## with it. The warning names the two settings that actually cause it.
  ##
  ## A SINGLE dead BAM among several is not this case and is left alone: that
  ## file's row draws flat and empty, which is the honest picture, and the
  ## others draw normally.
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

  ## The same ceiling the scan applies, reached independently because a caller
  ## can arrive with a saved percentile table and never open a BAM.
  .check_bam_limit(nrow(pct), "figure")

  cat("=== plot_scout_distribution ===\n")
  cat("Source    :", if (!is.null(scout)) "scout object, in memory" else in_dir, "\n")
  cat("Label     :", label, "\n")
  cat("Files     :", nrow(pct), "\n\n")


  ## ---------------------------------------------------------------------------
  ## Pool the strands. The cutoff is computed on both strands pooled and the
  ## figure follows it. The per-strand rows stay in scout_distribution.tsv.
  ## ---------------------------------------------------------------------------

  if ("Strand" %in% names(dist)) {
    pooled <- stats::aggregate(list(Weight = dist$Weight),
                               by = list(File = dist$File, Coverage = dist$Coverage),
                               FUN = sum)
  } else {
    pooled <- dist[, c("File", "Coverage", "Weight")]
  }
  pooled <- pooled[order(pooled$File, pooled$Coverage), ]

  ## File order: as given in the percentile table, so the figure and the table
  ## above it list the samples identically. That takes two steps and not one. The
  ## row a file is drawn on is FLIPPED below, because a continuous y axis puts row
  ## 1 at the bottom, and the axis labels are reversed to match. Either change
  ## alone is worse than neither: flipping the rows without reversing the labels
  ## puts each violin under another file's name.
  file_levels <- pct$File
  pooled <- pooled[pooled$File %in% file_levels, ]
  stopifnot("no rows left after matching the two tables on File" = nrow(pooled) > 0)


  ## ---------------------------------------------------------------------------
  ## Helpers. These are nested rather than hoisted, because nesting is what
  ## preserves .gap()'s <<- against .y_cur: hoisting .gap would send that
  ## superassignment to the global environment and the annotation block would
  ## silently misposition.
  ## ---------------------------------------------------------------------------

  ## Weighted standard deviation, for the bandwidth rule.
  .wsd <- function(x, w) {
    mu <- sum(x * w) / sum(w)
    sqrt(sum(w * (x - mu)^2) / sum(w))
  }

  ## Silverman's rule of thumb on the weighted sample. n is the total weight, that
  ## is the number of positions, not the number of distinct coverage values.
  .bandwidth <- function(x, w, iqr) {
    n <- sum(w)
    s <- .wsd(x, w)
    a <- min(s, iqr / 1.349)
    if (!is.finite(a) || a <= 0) a <- s
    if (!is.finite(a) || a <= 0) a <- 1
    0.9 * a * n^(-0.2)
  }

  ## Tukey whiskers on a weighted set of values: the extreme observed values still
  ## inside the 1.5 IQR fences. Hinges come from the percentile table.
  .whiskers <- function(cov, q1, q3) {
    iqr <- q3 - q1
    lo_fence <- q1 - 1.5 * iqr
    hi_fence <- q3 + 1.5 * iqr
    inside <- cov[cov >= lo_fence & cov <= hi_fence]
    if (!length(inside)) return(c(min(cov), max(cov)))
    c(min(inside), max(inside))
  }


  ## ---------------------------------------------------------------------------
  ## Per-file geometry
  ## ---------------------------------------------------------------------------

  x_max <- max(log2(pooled$Coverage))
  grid_n <- 512

  violin_df <- NULL
  box_df    <- NULL
  med_df    <- NULL
  whisk_df  <- NULL

  ## Row 1 is the bottom of the panel, so the first file in the table takes the
  ## HIGHEST row. Read with the reversed axis labels below, which name row 1 as
  ## the last file in the table.
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

  ## Okabe-Ito as a plain hex vector, so nothing needs installing. Reordered for
  ## FILLS: the canonical order starts with grey and puts a pale yellow fifth,
  ## both of which read poorly as a filled area below full opacity.
  ##
  ## RColorBrewer is not a declared dependency; it arrives with ggplot2. If it is
  ## absent, the base grDevices HCL analogues are used.
  .OKABE_ITO_FILL <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
                       "#56B4E9", "#E69F00", "#999999", "#F0E442")

  .resolve_palette <- function(name, n) {
    if (identical(name, "grey")) {
      return(rep("grey60", n))
    }
    if (identical(name, "Okabe-Ito")) {
      if (n <= length(.OKABE_ITO_FILL)) return(.OKABE_ITO_FILL[seq_len(n)])
      ## Beyond eight files the palette recycles rather than interpolating: the
      ## ninth row takes the first colour again. Interpolated qualitative colours
      ## look distinct without being reliably distinguishable, which is a worse
      ## failure than an honest repeat, and the fill encodes nothing anyway.
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

  ## Colourblind-safety note, emitted rather than buried: ColorBrewer flags its
  ## qualitative palettes safe only at three classes.
  if (PALETTE %in% c("Dark2", "Set2", "Paired") && length(file_levels) > 3) {
    message("Note: ColorBrewer flags '", PALETTE, "' colourblind-safe only at ",
            "n = 3; this figure has ", length(file_levels),
            " files. PALETTE = \"Okabe-Ito\" is safe at 8.")
  }

  pal <- .resolve_palette(PALETTE, length(file_levels))
  names(pal) <- file_levels


  ## ---------------------------------------------------------------------------
  ## The plot
  ## ---------------------------------------------------------------------------

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
    ## Expansion is zero here because the y view is pinned by coord_cartesian()
    ## below, from limits computed off the data alone. Leaving the expansion in
    ## would apply it a second time on top of those limits.
    ## REVERSED, and it must stay in step with row_of_file() above. Row 1 is the
    ## bottom of the panel and carries the LAST file in the percentile table, so
    ## the labels run in the opposite direction to the table's rows while the
    ## drawn order matches it.
    scale_y_continuous(breaks = seq_along(file_levels), labels = rev(file_levels),
                       expand = expansion(mult = c(0, 0))) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.04))) +
    labs(x = expression(log[2] * "(coverage)"), y = NULL, title = plot_title) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor   = element_blank(),
          axis.text.y        = element_text(hjust = 1),
          plot.title         = element_text(size = 11, face = "plain"))

  ## ---------------------------------------------------------------------------
  ## The values offered, computed from the percentile table
  ## ---------------------------------------------------------------------------

  ## Percentile columns actually present, for the console echo below. The drawn
  ## table computes the same intersection separately at pct_cols.
  pct_cols_preview <- intersect(c("P25", "Median", "P75", "P80", "P85", "P90"),
                                names(pct))

  lad_p25 <- .agg_percentile(pct$P25)
  lad_p50 <- .agg_percentile(pct$Median)
  lad_p75 <- .agg_percentile(pct$P75)
  lad_p80 <- .agg_percentile(pct$P80)

  ## The aggregates, by percentile name, so a pair can be looked up rather than
  ## hardcoded. P25 is computed because the table's aggregate row shows it, and
  ## is not among the offered pairs.
  agg_by_name <- c(P25 = lad_p25, Median = lad_p50, P75 = lad_p75, P80 = lad_p80)

  ## From .cutoff_pairs(), the single definition the writers also read, so the
  ## figure and the text artefact beside it cannot state different pairs.
  pairs_df <- .cutoff_pairs(pct)

  ## The stringent pair raises both cutoffs, so it retains less intergenic
  ## sequence, but it may still report MORE features because raising the low
  ## cutoff moves slice boundaries. Warn if the arithmetic stops holding,
  ## which would mean the labels are the wrong way round.
  if (!any(is.na(c(pairs_df$low, pairs_df$high)))) {
    if (!(pairs_df$low[2] >= pairs_df$low[1] && pairs_df$high[2] >= pairs_df$high[1])) {
      warning("The pair labelled '", pairs_df$pair[2], "' does not have both cutoffs ",
              "at or above the pair labelled '", pairs_df$pair[1], "'. Check ",
              "PAIR_LOW and PAIR_HIGH: the labels assume the second is the ",
              "stricter of the two.", call. = FALSE, immediate. = TRUE)
    }
  }
  ## Each pair must also be internally valid, high strictly above low. This is the
  ## condition suggest_cutoffs() warns on, checked here so the figure does not
  ## print a pair the module would object to.
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

  ## An aggregate of 0 or 1 imposes no floor, because every position in the
  ## distribution carries at least one read by construction. The unmarked-values
  ## message below reports it if one ever appears at an offered percentile.
  cat("\n")

  ## One vertical marker per percentile in MARKER_PCTS. Every marked value is a
  ## cell of the aggregate row, so a line can be traced to a visible number. The
  ## tag reads "Median" rather than "P50", matching the table heading.
  ##
  ## A value of 1 is not drawn: log2(1) is zero, so the marker would sit on the
  ## panel edge and read as the axis. The console reports any offered value left
  ## without a line.
  ##
  ## Tags carry a white backing because the markers are full-height lines, so a
  ## tag beside one marker is crossed by the next when two percentiles are close.
  y_top <- length(file_levels) + VIOLIN_H
  ## The values the two pairs actually use, for the reporting below. P75 appears in
  ## both pairs, as the high of the relaxed one and the low of the stringent one,
  ## so it is deduplicated or it would be counted twice.
  offered_vals <- agg_by_name[unique(c(PAIR_LOW, PAIR_HIGH))]

  ## Markers are drawn from all computed percentiles, not only the offered
  ## ones, so a landmark can be seen whether or not anyone would choose it.
  stopifnot("MARKER_PCTS must name computed percentiles" =
              all(MARKER_PCTS %in% names(agg_by_name)))

  line_df <- data.frame(
    value = unname(agg_by_name[MARKER_PCTS]),
    tag   = MARKER_PCTS,
    stringsAsFactors = FALSE)
  ## NA-safe: a BAM that yielded no reads gives NA percentiles by design, and
  ## subsetting with an NA index would produce a row of NAs, i.e. a phantom marker
  ## that ggplot2 silently drops with a warning about missing values.
  line_df <- line_df[!is.na(line_df$value) & line_df$value > 1, , drop = FALSE]
  line_df$x  <- log2(line_df$value)
  line_df$dy <- (seq_len(nrow(line_df)) - 1) * 0.17
  lab_base   <- y_top + 0.38

  ## Report every offered value with no line beside it, and why. Such a value
  ## cannot be located on the distribution by eye, so it should be a deliberate
  ## choice rather than something noticed later from the figure.
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

  ## ---------------------------------------------------------------------------
  ## The annotation block.
  ##
  ## Every line is its own layer, and every line is plotmath when BOLD_AGG is TRUE.
  ## That is forced rather than chosen: ggplot2 applies one fontface per text
  ## layer, so bolding a label needs plotmath, and plotmath has no line break, so a
  ## bold label followed by wrapping text cannot be one layer. Sending the
  ## unbolded lines through plotmath too keeps every left edge on the same
  ## renderer, which is what makes every line's first character align.
  ## ---------------------------------------------------------------------------

  ## plotmath helpers. A literal run is quoted; a bold run is wrapped in bold();
  ## paste() concatenates them. Spaces inside the quoted runs are preserved, which
  ## is why the padding sits inside the quotes.
  .q  <- function(s) sprintf('"%s"', s)
  .bd <- function(s) sprintf('bold("%s")', s)
  .pm <- function(...) paste0("paste(", paste(c(...), collapse = ", "), ")")

  ## One wrapped sentence, as plotmath lines, with a leading label emboldened on
  ## the first line only. fixed = TRUE on the label sub(), so a label containing
  ## a regex metacharacter cannot misbehave.
  .block_lines <- function(sentence, label, width) {
    ## strwrap() wraps WITH the indent, so continuation lines are measured
    ## against the space they will occupy. The leading spaces are then stripped
    ## and replaced by an x offset, because plotmath drops leading whitespace on
    ## some devices. width = NULL draws the sentence as a single line.
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

  ## One pair line: bold label, bold values, plain connective text.
  ##
  ##     Relaxed: Low <low>, High <high> (Median / P75)
  ##
  ## The percentile names are in parentheses so the pair traces to two cells of
  ## the aggregate row.
  .pair_line <- function(label, low, high, from) {
    plain <- sprintf("%s: Low %s, High %s (%s)", label, low, high, from)
    if (!isTRUE(BOLD_AGG)) return(list(text = plain, plain = plain))
    list(text = .pm(.bd(paste0(label, ":")),
                    .q(" Low "), .bd(low),
                    .q(", High "), .bd(high),
                    .q(sprintf(" (%s)", from))),
         plain = plain)
  }

  ## Layout. Spacing is a fraction of the plotted y range, so it is a constant
  ## physical distance whatever the number of BAMs (see ANN_LINE_FRAC).
  y_span <- (length(file_levels) + VIOLIN_H) - (1 - BOX_GAP - BOX_H / 2)
  line_h <- ANN_LINE_FRAC * y_span
  gap_h  <- ANN_GAP_FRAC  * y_span

  ## Assemble the block, top down, accumulating one row per drawn line.
  ann <- data.frame(text = character(0), plain = character(0), size = numeric(0),
                    indent = numeric(0), y = numeric(0), stringsAsFactors = FALSE)
  ## The block starts at lab_base plus the offset, in line heights, so the lift
  ## is the same physical distance whatever the row count. lab_base itself is
  ## untouched, which keeps the marker tags where they are.
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

  ## The note, below the second pair with a blank line between. Drawn when
  ## SHOW_ANN_NOTE is TRUE and the text is non-empty.
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

  ## Horizontal anchor: the widest line ends at ANN_RIGHT_FRAC unless
  ## ANN_X_FRAC overrides. Width is estimated from the visible character count
  ## scaled by each line's size, which is why every row carries its plain text.
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
  ## What the computed anchor WOULD have been, printed even when overridden, so the
  ## gap between the estimate and the override is visible rather than hidden. A
  ## large gap means ANN_CHAR_FRAC is mis-estimating the widest line.
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
               ## linewidth = 0, not label.size = 0, which is deprecated from
               ## ggplot2 3.5.0. Render-neutral: ggplot2 translates the old
               ## name before the layer data is built.
               fill = "white", linewidth = 0,
               label.padding = grid::unit(0.08, "lines")) +
    annotate("text", x = ann$x, y = ann$y, label = ann$text,
             parse = isTRUE(BOLD_AGG),
             hjust = 0, vjust = 1, size = ann$size, colour = "grey20")


  ## Pin the y view so the annotation cannot resize the distributions. Limits
  ## come from the data alone: boxes and whiskers set the floor, marker tags the
  ## ceiling. coord_cartesian() limits the view rather than the data, so nothing
  ## is dropped and no rows are removed with a warning.
  y_floor  <- 1 - BOX_GAP - BOX_H / 2
  y_ceil   <- lab_base
  y_extent <- y_ceil - y_floor
  y_view   <- c(y_floor - 0.06 * y_extent,
                y_ceil  + ANN_HEADROOM_MULT * y_extent)

  ## Does the block fit in the headroom? The first line is drawn with vjust = 1, so
  ## its top edge is its y position; the last line's baseline is the block's foot.
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


  ## ---------------------------------------------------------------------------
  ## The table, drawn as a ggplot so no extra package is needed
  ## ---------------------------------------------------------------------------

  pct_cols <- intersect(c("P25", "Median", "P75", "P80", "P85", "P90"), names(pct))
  tab <- pct[, c("File", pct_cols), drop = FALSE]
  n_bam <- nrow(tab)

  ## Appended as an ordinary row so it is formatted by the same fmt() call as
  ## the column it sits under and cannot drift from it.
  agg_row <- data.frame(File = AGG_ROW_LABEL, stringsAsFactors = FALSE)
  for (nm in pct_cols) agg_row[[nm]] <- .agg_percentile(tab[[nm]])
  tab <- rbind(tab, agg_row)

  tab_cols <- names(tab)

  fmt <- function(v, nm) {
    if (nm == "File") return(as.character(v))
    format(as.numeric(v), trim = TRUE)
  }

  ## Row 1 is the header and rows 2 to n_bam + 1 the BAMs. The aggregate sits
  ## AGG_ROW_GAP further down again, so the separating rule has room to read as a
  ## rule rather than as an underline on the last BAM.
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

  ## The rule between the last BAM and the aggregate, placed midway across the
  ## gap the row offset opened.
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


  ## ---------------------------------------------------------------------------
  ## Combine and write
  ## ---------------------------------------------------------------------------

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

    ## Absolute, matching the two lines printed above and the four paths the
    ## table writer returns, so every path this module hands back resolves from
    ## anywhere.
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

  ## The two ggplot objects are returned so a caller can inspect or re-render the
  ## figure without scanning again.
  invisible(list(table = p_table, plot = p_plot, files = files))
}
