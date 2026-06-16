# Commit 7: Behavioural Fixes

## Commit Abstract: Correct the plus-strand UTR length filter to apply the UTR cut-off, and pin featureCounts' multi-mapping and read-pair counting so quantification is deterministic across Rsubread versions.

### Commit Summary:
Correct UTR filter, make counts deterministic

### Commit Description:
fix: set countMultiMappingReads and countReadPairs explicitly for
deterministic counts across Rsubread versions, and correct the
plus-strand UTR length filter

## Plus-strand UTR length filter mismatch

**Issue Summary:** The plus-strand `UTR_calc()` call passes `min_sRNA_length` where the signature expects `min_UTR_length`, so plus-strand UTRs filter at the wrong cut-off.

**Solution Summary:** Change the fourth positional argument on line 347 from `min_sRNA_length` to `min_UTR_length`, matching the signature and the correct minus-strand call.

**Note:** Output is not identical; this is a behavioural fix that removes plus-strand putative UTRs whose length falls between the two parameters.

### feature_file_editor.R

#### Line 347 — plus-strand `UTR_calc()` call in `feature_file_editor`

The plus-strand call passes the sRNA-length cut-off where the signature expects the UTR-length cut-off; the minus-strand call already passes the correct argument and is unchanged.

```
# commit 6 (previous)
# line 347
    plus_UTR <- UTR_calc(maj_plus_features, "+", plus_strand_peaks, min_sRNA_length)
```

```
# commit 7 (this commit)
# line 347
    plus_UTR <- UTR_calc(maj_plus_features, "+", plus_strand_peaks, min_UTR_length)
```

---

## Set `countMultiMappingReads` and `countReadPairs` explicitly

**Issue Summary:** `count_features()` sets neither `countMultiMappingReads` nor `countReadPairs`, so the same baerhunter code produces version-dependent counts after Rsubread's R-side defaults changed between 2019 and 2026.

**Solution Summary:** Expose both as snake_case `count_features()` arguments forwarded to `featureCounts()`, defaulting to `FALSE` and `TRUE` so output is deterministic across Rsubread versions and reproduces 2019 behaviour.

**Note:** Not identity-preserving relative to current unpatched 2.x code; multi-mapper counts at repetitive loci drop, and a `Rsubread (>= 2.4.3)` requirement is documented (an in-code annotation in this commit, plus a `DESCRIPTION` `Imports` floor at the documentation commit) because `2.4.3` is the lowest Rsubread version confirmed to carry `countReadPairs`.

### count_features.R

Relative to commit 6, this change adds the `count_multi_mapping_reads` and `count_read_pairs` arguments, their two forwarding lines and the preceding version comment, the two matching roxygen entries, and the correction to the `...` roxygen entry. Every other argument, forwarding line and roxygen entry shown is unchanged from commit 6 and is reproduced only for context.

#### Signature

```
# commit 6 (previous)
# lines 105-117
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
                           ...){
```

```
# commit 7 (this commit)
# lines 107-121
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
```

#### featureCounts call

```
# commit 6 (previous)
# lines 146-156
  fc <- featureCounts(bam_files,
                      annot.ext = nsaf_df,
                      chrAliases = chromosome_alias_file,
                      strandSpecific = strand_specific,
                      isPairedEnd = paired_end,
                      allowMultiOverlap = TRUE,
                      fraction = TRUE,
                      largestOverlap = largest_overlap,
                      fracOverlapFeature = frac_overlap_feature,
                      read2pos = read_to_pos,
                      ...)
```

```
# commit 7 (this commit)
# lines 150-164
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
```

#### Roxygen @param entries

```
# commit 6 (previous)
# lines 93-97
#' @param excl_rna A boolean indicating if misc RNA features (rRNA, tRNA) are excluded from quantification. (Defaults=TRUE)
#' @param largest_overlap A boolean; if TRUE, assigns each read to the feature with the largest number of overlapping bases. Maps to featureCounts largestOverlap. Combined with the package's `fraction = TRUE`, Rsubread >= 2.14.0 is recommended, since earlier versions silently miscount that combination. (Default: FALSE)
#' @param frac_overlap_feature Minimum fraction of a feature that must be overlapped before a read is assigned to it. Maps to featureCounts fracOverlapFeature. (Default: 0)
#' @param read_to_pos Reduce each read to a single base before counting: 5 for the 5' end, 3 for the 3' end, or NULL to count the whole read. Maps to featureCounts read2pos. (Default: NULL)
#' @param ... Optional parameters passed on to featureCounts()  Default: allowMultiOverlap = TRUE, fraction = TRUE
```

```
# commit 7 (this commit)
# lines 93-99
#' @param excl_rna A boolean indicating if misc RNA features (rRNA, tRNA) are excluded from quantification. (Defaults=TRUE)
#' @param largest_overlap A boolean; if TRUE, assigns each read to the feature with the largest number of overlapping bases. Maps to featureCounts largestOverlap. Combined with the package's `fraction = TRUE`, Rsubread >= 2.14.0 is recommended, since earlier versions silently miscount that combination. (Default: FALSE)
#' @param frac_overlap_feature Minimum fraction of a feature that must be overlapped before a read is assigned to it. Maps to featureCounts fracOverlapFeature. (Default: 0)
#' @param read_to_pos Reduce each read to a single base before counting: 5 for the 5' end, 3 for the 3' end, or NULL to count the whole read. Maps to featureCounts read2pos. (Default: NULL)
#' @param count_multi_mapping_reads A boolean; if FALSE, reads mapping to multiple locations are excluded from counts. Maps to featureCounts countMultiMappingReads. (Default: FALSE, reproducing 2019 behaviour)
#' @param count_read_pairs A boolean; for paired-end data, if TRUE each fragment is counted once rather than each mate separately. Ignored for single-end data. Maps to featureCounts countReadPairs. Requires Rsubread >= 2.4.3. (Default: TRUE)
#' @param ... Optional parameters passed on to featureCounts(). Note that allowMultiOverlap and fraction are set internally to TRUE and cannot be overridden via ...
```
