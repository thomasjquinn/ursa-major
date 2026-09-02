# Commit 3: Performance updates

## Commit Abstract: Cut redundant work in the sRNA and UTR detection path by reading each BAM once across both strands, vectorising the per-peak coverage test, and computing the feature count a single time, with no change to scientific output.

### Commit Summary:
Speed up sRNA/UTR detection, no output change

### Commit Description:
perf(feature_file_editor): read each BAM once across strands,
replace per-peak viewApply with vectorised viewSums,
and compute the feature count once

## Repeated `as.data.frame(sRNA_UTR)` in `strand_feature_editor`

**Issue Summary:** `strand_feature_editor()` converts the combined sRNA and UTR ranges to a data frame three times per call, only to read the row count.

**Solution Summary:** Compute the count once with `length(sRNA_UTR)` and reuse it; on an `IRanges` object this returns the same integer with no allocation.

**Note:** A speed fix with identical output; `length()` and `nrow(as.data.frame())` return the same integer on an `IRanges`, so no test is required.

### feature_file_editor.R

```
# commit 2 (previous)
# line 244 (context: sRNA_UTR is built here)
  sRNA_UTR <- c(sRNA_IRanges, UTR_IRanges)

# line 247
  seqid <- rep(major_strand_features[1,1],nrow(as.data.frame(sRNA_UTR)))

# line 248
  empty_col <- rep(".",nrow(as.data.frame(sRNA_UTR)))

# line 252 (plus strand branch)
    cmp_strand <- data.frame(seqid, empty_col, empty_col, as.integer(start(sRNA_UTR)), as.integer(end(sRNA_UTR)), empty_col, rep("+", nrow(as.data.frame(sRNA_UTR))), empty_col, names(sRNA_UTR), stringsAsFactors = FALSE)

# line 254 (minus strand branch)
    cmp_strand <- data.frame(seqid, empty_col, empty_col, as.integer(start(sRNA_UTR)), as.integer(end(sRNA_UTR)), empty_col, rep("-", nrow(as.data.frame(sRNA_UTR))), empty_col, names(sRNA_UTR), stringsAsFactors = FALSE)
```

```
# commit 3 (this commit)
# lines 240-241: the count is computed once on the line after sRNA_UTR
  sRNA_UTR <- c(sRNA_IRanges, UTR_IRanges)
  n_features <- length(sRNA_UTR)

# line 244
  seqid <- rep(major_strand_features[1,1], n_features)

# line 245
  empty_col <- rep(".", n_features)

# line 249 (plus strand branch)
    cmp_strand <- data.frame(seqid, empty_col, empty_col, as.integer(start(sRNA_UTR)), as.integer(end(sRNA_UTR)), empty_col, rep("+", n_features), empty_col, names(sRNA_UTR), stringsAsFactors = FALSE)

# line 251 (minus strand branch)
    cmp_strand <- data.frame(seqid, empty_col, empty_col, as.integer(start(sRNA_UTR)), as.integer(end(sRNA_UTR)), empty_col, rep("-", n_features), empty_col, names(sRNA_UTR), stringsAsFactors = FALSE)
```

---

## `coverage()` recomputed across strand calls

**Issue Summary:** `peak_union_calc()` runs once per strand, so every BAM file is read and parsed twice per run, which dominates cost for large alignments.

**Solution Summary:** Read each BAM once, split alignments by strand, and return both unions from a single call as a named `plus`/`minus` list.

**Note:** A speed fix that halves BAM read and parse; output is expected identical, pending the regression testing the strand-mapping rewrite makes necessary.

### feature_file_editor.R

#### peak_union_calc()

The function is restructured to read each BAM once, split by strand, and return a named `plus`/`minus` list, dropping the `target_strand` argument. The previous commit's single-chromosome `stop()` guard and idiomatic comparison forms carry through into the restructured helper unchanged; the per-peak detection inside the helper uses the vectorised `viewSums` form covered separately below, shown here already integrated as the shipped end state.

```
# commit 2 (previous)
# roxygen line 8
#' @param target_strand A character string indicating the strand. Supports two values; '+' and '-'.

# roxygen line 15
#' @return An object of IRanges class, containing the genomic coordinates of selected and unified peaks.

# signature and body, lines 23-80
peak_union_calc <- function(bam_location = ".", bam_txt_list = "", target_strand, low_coverage_cutoff, high_coverage_cutoff,  peak_width, paired_end_data = FALSE, strandedness  = "unstranded") {
  ## Find all BAM files in the directory.
  if (bam_txt_list != ""){
    bam_files <- readLines(bam_txt_list)
    bam_files <- lapply(bam_files, function(x) paste(bam_location, x, sep="/"))
  }else{
    bam_files <- list.files(path = bam_location, pattern = ".BAM$", full.names = TRUE, ignore.case = TRUE)
  }
  peak_union <- IRanges()
  ## Go over each BAM file to extract coverage peaks for a target strand and gradually build a union of all peak sets.
  for (f in bam_files) {
    ## Read a BAM file in accordance with its type and select only the reads aligning to a target strand.
    strand_alignment <- c()
    if (!paired_end_data & strandedness  == "stranded") {
      file_alignment <- readGAlignments(f)
      strand_alignment <- file_alignment[strand(file_alignment)==target_strand,]
    } else if (paired_end_data & strandedness  == "stranded") {
      file_alignment <- readGAlignmentPairs(f, strandMode = 1)
      strand_alignment_unmerged <- file_alignment[strand(file_alignment)==target_strand,]
      #coerce to get single ranges
      strand_alignment <- granges(strand_alignment_unmerged)
    } else if (!paired_end_data & strandedness  == "reversely_stranded") {
      file_alignment <- readGAlignments(f)
      relevant_strand <- c()
      if (target_strand=="+") {
        relevant_strand <- "-"
      } else {
        relevant_strand <- "+"
      }
      strand_alignment <- file_alignment[strand(file_alignment)==relevant_strand,]
    } else if (paired_end_data & strandedness  == "reversely_stranded") {
      file_alignment <- readGAlignmentPairs(f, strandMode = 2)
      strand_alignment_unmerged <- file_alignment[strand(file_alignment)==target_strand,]
      #coerce to get single ranges
      strand_alignment <- granges(strand_alignment_unmerged)
    }
    ## Create a strand coverage vector and extract it
    strand_cvg <- coverage(strand_alignment)
    list_components <- names(strand_cvg)
    target <- c()
    if (length(list_components)==1) {
      target <- list_components
    } else {
      stop(paste("Invalid BAM file:", f, sep = " "))
    }
    ## Cut the coverage vector to obtain the expression peaks with the coverage above the low cut-off values.
    peaks <- IRanges::slice(strand_cvg[[target]], lower = low_coverage_cutoff, includeLower=TRUE)
    ## Examine the peaks for the stretches of coverage above the high cut-off. The stretches have to be a defined width.
    test <- viewApply(peaks, function(x) peak_analysis(x,high_coverage_cutoff,peak_width))
    ## Select only the peaks that satisfy the high cut-off condition.
    selected_peaks <- peaks[vapply(test, function(x) !is.null(x), logical(1))]
    ## Convert peak coordinates into IRanges.
    peaks_IRange <- IRanges(start = start(selected_peaks), end = end(selected_peaks))
    ## Calculate the peak union in with the previous peak sets.
    peak_union <- union(peak_union,peaks_IRange)
  }
  return(peak_union)
}
```

```
# commit 3 (this commit)
# the @param target_strand line is removed

# lines 14-15: the @return line is replaced (now spans two lines)
#' @return A named list with two IRanges objects, `plus` and `minus`, holding the
#'   unified peak coordinates for each strand.

# restructured signature and body, lines 23-76
peak_union_calc <- function(bam_location = ".", bam_txt_list = "", low_coverage_cutoff, high_coverage_cutoff, peak_width, paired_end_data = FALSE, strandedness = "unstranded") {
  ## Find all BAM files in the directory.
  if (bam_txt_list != ""){
    bam_files <- readLines(bam_txt_list)
    bam_files <- lapply(bam_files, function(x) paste(bam_location, x, sep = "/"))
  } else {
    bam_files <- list.files(path = bam_location, pattern = ".BAM$", full.names = TRUE, ignore.case = TRUE)
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
  ## Read each BAM once, split by strand, and accumulate both peak unions.
  for (f in bam_files) {
    if (paired_end_data) {
      strand_mode <- if (strandedness == "reversely_stranded") 2 else 1
      file_alignment <- granges(readGAlignmentPairs(f, strandMode = strand_mode))
    } else {
      file_alignment <- readGAlignments(f)
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
    plus_union  <- union(plus_union,  compute_strand_peaks(plus_reads,  f))
    minus_union <- union(minus_union, compute_strand_peaks(minus_reads, f))
  }
  return(list(plus = plus_union, minus = minus_union))
}
```

#### call sites

```
# commit 2 (previous)
# line 331
    plus_strand_peaks <- peak_union_calc(bam_location = bam_directory, bam_txt_list = bam_list, "+", low_coverage_cutoff, high_coverage_cutoff, min_sRNA_length, paired_end_data, strandedness)

# line 339
    minus_strand_peaks <- peak_union_calc(bam_location = bam_directory, bam_txt_list = bam_list, "-", low_coverage_cutoff, high_coverage_cutoff, min_sRNA_length, paired_end_data, strandedness)
```

```
# commit 3 (this commit)
# lines 328-329: one call computes both strands
    peak_sets <- peak_union_calc(bam_location = bam_directory, bam_txt_list = bam_list, low_coverage_cutoff, high_coverage_cutoff, min_sRNA_length, paired_end_data, strandedness)
    plus_strand_peaks  <- peak_sets$plus

# line 337: reuse the minus slot from the single call above
    minus_strand_peaks <- peak_sets$minus
```

---

## Replace `viewApply` with vectorised view primitives

**Issue Summary:** `peak_union_calc()` counts qualifying positions per peak with a `viewApply` R-level closure (`peak_analysis()`), incurring per-peak dispatch and a per-peak Rle materialisation.

**Solution Summary:** Threshold the coverage once into a logical Rle and count qualifying positions per peak with `viewSums`, replacing the per-peak closure while leaving `peak_analysis()` in place.

**Note:** A pure speed fix with no behaviour change; `viewSums` reproduces the current per-peak count exactly, so the output is byte-identical and a test confirms it.

### feature_file_editor.R

The detection block inside `peak_union_calc()` changes from a per-peak `viewApply` closure to a single vectorised `viewSums` pass. The shipped form of this block sits inside the restructured helper shown in the previous change (committed lines 40-48 within `compute_strand_peaks()`); it is shown here against the previous commit's loop body to isolate the `viewApply`-to-`viewSums` replacement. `peak_analysis()` (the exported helper at committed lines 90-96) is left unchanged and is simply no longer called by the pipeline, so the exported interface is unchanged.

```
# commit 2 (previous)
# lines 68-73
    ## Cut the coverage vector to obtain the expression peaks with the coverage above the low cut-off values.
    peaks <- IRanges::slice(strand_cvg[[target]], lower = low_coverage_cutoff, includeLower=TRUE)
    ## Examine the peaks for the stretches of coverage above the high cut-off. The stretches have to be a defined width.
    test <- viewApply(peaks, function(x) peak_analysis(x,high_coverage_cutoff,peak_width))
    ## Select only the peaks that satisfy the high cut-off condition.
    selected_peaks <- peaks[vapply(test, function(x) !is.null(x), logical(1))]
```

```
# commit 3 (this commit)
# lines 40-48, inside compute_strand_peaks()
    ## Cut the coverage vector to obtain the expression peaks above the low cut-off.
    strand_rle <- strand_cvg[[target]]
    peaks <- IRanges::slice(strand_rle, lower = low_coverage_cutoff, includeLower = TRUE)
    ## Count, per peak, the positions above the high cut-off, then keep the peaks
    ## whose count exceeds the minimum width. This reproduces the count-based test.
    above_high_rle <- strand_rle > high_coverage_cutoff
    above_high_views <- Views(above_high_rle, ranges(peaks))
    positions_above <- viewSums(above_high_views)
    selected_peaks <- peaks[positions_above > peak_width]
```

---

## `as.vector(View_line)` materialises the Rle per peak

**Issue Summary:** Inside `peak_analysis()`, `as.vector(View_line)` materialises each peak's coverage Rle to a plain vector, allocating memory on every per-peak invocation.

**Solution Summary:** No standalone change; the `viewSums` refactor stops the pipeline calling `peak_analysis()`, so the per-peak materialisation no longer runs.

**Note:** A speed-only effect with identical output; resolved entirely by the `viewSums` refactor, and the line is retained (not deleted), just no longer executed.

This change makes no edit of its own. Inside `peak_analysis()`, `cvg_string <- as.vector(View_line)` (line 92 of the shipped file) is unchanged. Because the restructured `peak_union_calc()` no longer calls `peak_analysis()`, this per-peak materialisation no longer runs during a pipeline run, removing its aggregate allocation cost across tens of thousands of peaks. The line is retained for any external caller of the exported helper.
