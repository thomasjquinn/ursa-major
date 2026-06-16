# Commit 8: Paired-End Methodology

## Commit Abstract: Filter low-quality, multi-mapping, secondary and supplementary reads out of single-end and paired-end peak detection through a configurable ScanBamParam (default mapqFilter 10 plus flag exclusions), add a fragment-versus-footprint paired-end coverage model, make the per-strand union accumulators explicit as IRanges::union, and warn when an over-strict filter empties a BAM.

### Commit Summary:
Add read-quality filter and paired-end coverage

### Commit Description:
feat(peak-detection): add a configurable read-quality filter to
single-end and paired-end peak detection with a fragment/footprint
coverage model and explicit IRanges::union calls

## Read-quality filter for paired-end peak detection

**Issue Summary:** `peak_union_calc()` reads every BAM alignment with no quality or flag filtering, so low-mapping-quality, multi-mapping, secondary and supplementary reads inflate the coverage that peaks are called from.

**Solution Summary:** Add a user-facing `scanbamparam` argument defaulting to an internally built filter (`mapqFilter = 10` plus flag exclusions), with the threshold also settable through a dedicated `mapqFilter` argument, applied to the paired-end read call; the single-end branch is extended by the single-end read-quality filter. The same patch makes the two bare `union` accumulator calls explicit as `IRanges::union`.

**Note:** Output is not byte-identical under the new default on paired-end data, so this is a methodological read-quality fix and carries the largest behaviour-change footprint of any approved change. The `IRanges::union` part is identity-preserving on its own. On the single-end Cortes benchmark this change alone is a no-op; the observable single-end effect is carried by the single-end read-quality filter.

### feature_file_editor.R

#### `peak_union_calc()` signature

The signature gains the `scanbamparam` filter argument; the trailing `mapqFilter` and `coverage_model` arguments shown here are added by the convenience-argument and paired-end coverage changes.

```
# commit 7 (previous)
# line 23
peak_union_calc <- function(bam_location = ".", bam_txt_list = "", low_coverage_cutoff, high_coverage_cutoff, peak_width, paired_end_data = FALSE, strandedness = "unstranded") {
```

```
# commit 8 (this commit)
# line 52
peak_union_calc <- function(bam_location = ".", bam_txt_list = "", low_coverage_cutoff, high_coverage_cutoff, peak_width, paired_end_data = FALSE, strandedness = "unstranded", scanbamparam = NULL, mapqFilter = 10, coverage_model = c("fragment", "footprint")) {
```

#### Default read-quality filter

A default `ScanBamParam` is built before the BAM loop when the caller passes none. The single-end `else` branch is added by the single-end read-quality filter; the combined `if`/`else` form is shown here.

```
# commit 7 (previous)
# no commit 7 equivalent; new block, inserted after the minus_union initialiser (commit 7 line 53)
```

```
# commit 8 (this commit)
# lines 85-113 (filter-builder comment shown abbreviated to its first line; comment 85-96, if-block 97-113)
  ## Build the default read-quality filter when none is supplied.
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
```

#### Convenience `mapqFilter` argument

The threshold is exposed as a dedicated `mapqFilter` argument (default 10) on both `peak_union_calc()` and `feature_file_editor()`, so a user can change the value without rebuilding the whole `ScanBamParam` and risking the loss of the flag exclusions. It feeds the default-builder only: when `scanbamparam = NULL`, the builder uses `mapqFilter = mapqFilter`; when the caller supplies their own `scanbamparam`, that object takes precedence and `mapqFilter` is ignored. So a call with neither argument runs at 10, `mapqFilter = N` keeps the default flags with a different threshold, and `scanbamparam = ScanBamParam(...)` is full control. `mapqFilter = NA` disables mapping-quality filtering while keeping the flag exclusions.

```
# commit 8 (this commit)
# the argument on both signatures (default 10); the builder consumes it (line 112)
    scanbamparam <- ScanBamParam(flag = scanbamflag, mapqFilter = mapqFilter)
```

#### Zero-read safety net

A too-high `mapqFilter` is made to fail loudly rather than silently: an accumulator before the loop, a per-BAM check inside it, and one summary warning after it, gated so it never fires under a permissive `ScanBamParam()`.

```
# commit 7 (previous)
# no commit 7 equivalent; new code at the three locations below
```

```
# commit 8 (this commit)
# line 84 (accumulator, alongside the union accumulators)
  empty_bams  <- character(0)

# line 130 (inside the loop, after the per-BAM read)
    if (length(file_alignment) == 0L) empty_bams <- c(empty_bams, f)

# lines 144-151 (after the loop, before the return)
  ## Warn once if an active mapqFilter left any BAM with no reads.
  if (length(empty_bams) > 0L && !is.null(scanbamparam) && !is.na(bamMapqFilter(scanbamparam))) {
    warning(paste0(length(empty_bams), " BAM file(s) yielded no reads after filtering with ",
                   "mapqFilter = ", bamMapqFilter(scanbamparam), " (e.g. ", basename(empty_bams[1]),
                   "). This usually means the threshold exceeds the aligner's maximum MAPQ ",
                   "(bwa aln ~37, bowtie2 ~42); lower mapqFilter or set it to NA."),
            call. = FALSE, immediate. = TRUE)
  }
```

#### `union` accumulators

The two bare per-strand accumulator calls are made explicit to force resolution to the IRanges method.

```
# commit 7 (previous)
# lines 72-73
    plus_union  <- union(plus_union,  compute_strand_peaks(plus_reads,  f))
    minus_union <- union(minus_union, compute_strand_peaks(minus_reads, f))
```

```
# commit 8 (this commit)
# lines 141-142
    plus_union  <- IRanges::union(plus_union,  compute_strand_peaks(plus_reads,  f))
    minus_union <- IRanges::union(minus_union, compute_strand_peaks(minus_reads, f))
```

#### `feature_file_editor()` signature

The wrapper signature gains the same `scanbamparam` argument; the trailing `mapqFilter` and `coverage_model` arguments shown here are added by the convenience-argument and paired-end coverage changes.

```
# commit 7 (previous)
# line 338
feature_file_editor <- function(bam_directory = ".", bam_list = "", original_annotation_file, annot_file_dir = ".", output_file, original_sRNA_annotation, low_coverage_cutoff, high_coverage_cutoff, min_sRNA_length, min_UTR_length, paired_end_data = FALSE, strandedness  = "stranded") {
```

```
# commit 8 (this commit)
# line 444
feature_file_editor <- function(bam_directory = ".", bam_list = "", original_annotation_file, annot_file_dir = ".", output_file, original_sRNA_annotation, low_coverage_cutoff, high_coverage_cutoff, min_sRNA_length, min_UTR_length, paired_end_data = FALSE, strandedness  = "stranded", scanbamparam = NULL, mapqFilter = 10, coverage_model = c("fragment", "footprint")) {
```

#### `peak_union_calc()` call in `feature_file_editor()`

The single call threads the new arguments through as named arguments; `coverage_model` is added by the paired-end coverage change.

```
# commit 7 (previous)
# line 342
    peak_sets <- peak_union_calc(bam_location = bam_directory, bam_txt_list = bam_list, low_coverage_cutoff, high_coverage_cutoff, min_sRNA_length, paired_end_data, strandedness)
```

```
# commit 8 (this commit)
# line 448
    peak_sets <- peak_union_calc(bam_location = bam_directory, bam_txt_list = bam_list, low_coverage_cutoff, high_coverage_cutoff, min_sRNA_length, paired_end_data, strandedness, scanbamparam = scanbamparam, mapqFilter = mapqFilter, coverage_model = coverage_model)
```

---

## Read-quality filter for single-end peak detection

**Issue Summary:** Dr Nobeli's drafted read-quality filter is wired only to the paired-end read call, so the single-end branch (`readGAlignments(f)`, line 60) still reads every alignment unfiltered.

**Solution Summary:** Extend the same `scanbamparam` filter to the single-end read call, and generalise the paired-end default-builder to an `if`/`else` with a single-end default flag set that drops the three paired-only flags which would otherwise reject every single-end read.

**Note:** Not identity-preserving on single-end data under the new default. This is the change that affects the single-end Cortes benchmark, where the paired-end filter itself has no effect.

### feature_file_editor.R

#### Single-end read call

The single-end read gains the filter parameter. The single-end `else` branch of the default-builder it relies on is shown in the paired-end filter's default-builder block, where the combined `if`/`else` form appears.

```
# commit 7 (previous)
# line 60
      file_alignment <- readGAlignments(f)
```

```
# commit 8 (this commit)
# line 128
      file_alignment <- readGAlignments(f, param = scanbamparam)
```

---

## Paired-end coverage includes the inter-mate gap

**Issue Summary:** For paired-end data, `peak_union_calc()` collapses each read pair to a single range with `granges()`, so coverage is incremented across the unsequenced inter-mate gap as well as the two read footprints. This is the documented fragment-span behaviour, equivalent to deepTools `bamCoverage --extendReads`.

**Solution Summary:** Expose the coverage model as a `coverage_model` argument. The default `"fragment"` preserves the current `granges()` behaviour exactly. The `"footprint"` option computes coverage from the mate footprints with `unlist(grglist(...))`, excluding the inter-mate gap, which is the model standard RNA-seq practice favours.

**Note:** Identity-preserving by default. The footprint option is not identity-preserving on paired-end data, but it is reached only by explicitly passing `coverage_model = "footprint"`. Single-end data is unaffected and the parameter is a no-op on the single-end path.

### feature_file_editor.R

#### `coverage_model` argument check

The argument is resolved at the top of the function body. The `coverage_model` argument itself is added to both signatures, shown in the paired-end read-quality filter's signature blocks.

```
# commit 7 (previous)
# no commit 7 equivalent; new line at the top of the function body
```

```
# commit 8 (this commit)
# line 53
  coverage_model <- match.arg(coverage_model)
```

#### Paired-end coverage branch

The single read of each pair is bound once and converted by the chosen coverage model. `param = scanbamparam` on the read is added by the paired-end read-quality filter; this change adds the `read_pairs` binding and the fragment/footprint branch.

```
# commit 7 (previous)
# lines 56-58
    if (paired_end_data) {
      strand_mode <- if (strandedness == "reversely_stranded") 2 else 1
      file_alignment <- granges(readGAlignmentPairs(f, strandMode = strand_mode))
```

```
# commit 8 (this commit)
# lines 116-126 (coverage_model comment at 119-121 elided)
    if (paired_end_data) {
      strand_mode <- if (strandedness == "reversely_stranded") 2 else 1
      read_pairs <- readGAlignmentPairs(f, strandMode = strand_mode, param = scanbamparam)
      if (coverage_model == "footprint") {
        file_alignment <- unlist(grglist(read_pairs))
      } else {
        file_alignment <- granges(read_pairs)
      }
```

---

## The default `mapqFilter` value and per-aligner guidance

The package default is `mapqFilter = 10`, decided by Dr Nobeli (June 2026), superseding the interim 1 (7 June 2026). A threshold of 10 keeps uniquely-mapped reads on the common bacterial aligners while discarding the ambiguous and multi-mapping reads that cluster at MAPQ 0 to 3, and it matches the long-standing HTSeq default. MAPQ is aligner-specific, so the value must not exceed the MAPQ an aligner gives a uniquely-mapped read, or coverage empties. Maximum recommended values for the bacterial aligners:

| Aligner | Unique-read MAPQ | Max `mapqFilter` |
| --- | --- | --- |
| BWA-MEM | 60 | 60 |
| minimap2 | 60 | 60 |
| Rsubread align/subjunc | 40 | 40 |
| BWA aln (backtrack) | ~37 | 37 |
| Bowtie2 | 42 (non-monotonic) | use 1 |

Bowtie2 is the documented exception: its non-monotonic MAPQ scores some uniquely-mapped reads carrying mismatches 3 or 8, below 10, so Bowtie2 users should set `mapqFilter = 1`. STAR and HISAT2 are excluded as non-bacterial (splice-aware) aligners.

Testing is pinned at `mapqFilter = 1`, separate from the package default of 10, because the NH-tagged Test 3 BAMs are Bowtie2 (where 10 is unsafe) and 1 holds the already-characterised baselines. The shipped default (10) and the test value (1) are different numbers serving different roles, so raising the default does not change the testing methodology. The zero-read safety-net warning catches a too-high value that empties a BAM; it does not catch the silent Bowtie2 unique-with-mismatch case, which is why the per-aligner guidance is documented rather than left to the warning.

## Inline documentation (Commit 8)

The roxygen and inline comments for the filter and coverage parameters were brought forward from the final documentation commit to Commit 8, because users meet these parameters for the first time here. Both `peak_union_calc()` and `feature_file_editor()` carry `@param` blocks for `scanbamparam`, `mapqFilter`, and `coverage_model`.

```
# roxygen on both functions (abbreviated)
#' @param scanbamparam An optional Rsamtools::ScanBamParam object giving full
#'   control of the BAM read filter. When NULL (the default) an internal filter
#'   is built from mapqFilter plus alignment-flag exclusions (unmapped,
#'   QC-failing, secondary and supplementary dropped; paired-end also requires
#'   a properly paired read with a mapped mate). A supplied scanbamparam takes
#'   precedence and mapqFilter is then ignored.
#' @param mapqFilter Integer. Minimum MAPQ a read must have to be retained, used
#'   only when scanbamparam = NULL. Default 10, which keeps uniquely-mapped reads
#'   on the common bacterial aligners while dropping the MAPQ 0 to 3 tier. Do NOT
#'   set above the MAPQ your aligner gives a unique read, or coverage is empty.
#'   Ceilings: BWA-MEM 60, minimap2 60, Rsubread align/subjunc 40, BWA aln 37,
#'   Bowtie2 42 (but see note). Note: Bowtie2's non-monotonic MAPQ scores some
#'   unique-with-mismatch reads 3 or 8, below 10, so Bowtie2 users should set
#'   mapqFilter = 1. Set mapqFilter = NA to disable filtering while keeping the
#'   flag exclusions. A filter that retains no reads triggers a warning.
#' @param coverage_model Character, "fragment" or "footprint"; paired-end only,
#'   ignored for single-end. "fragment" (default) counts each pair from leftmost
#'   to rightmost base (the inter-mate gap is covered); "footprint" counts only
#'   the mates' aligned blocks (gap uncovered).
```

The default-builder and the paired-end coverage branch also carry inline comments restating the per-aligner ceilings, the Bowtie2 caveat, and the fragment-versus-footprint distinction at the point of use.
