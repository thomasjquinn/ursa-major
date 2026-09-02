# Commit 6: Additional Syntax and Robustness

## Commit Abstract: Add three opt-in featureCounts overlap parameters to count_features, inline IRanges::match in the sRNA and UTR predictors, remove the dead feature_type parameter from tpm_normalisation, and apply small robustness fixes (single-branch read_annotation_file, unquoted table output, and an invisible path return from tpm_flagging).

### Commit Summary:
Additional syntax and robustness fixes

### Commit Description:
feat(count_features): expose largestOverlap, fracOverlapFeature and
read2pos, inline IRanges::match in sRNA_calc and UTR_calc, remove the
dead feature_type parameter, collapse read_annotation_file to one
branch, write count and TPM tables without quoting, and return the
output path from tpm_flagging

## Local `match` shadowing in `sRNA_calc` and `UTR_calc`

**Issue Summary:** `sRNA_calc` and `UTR_calc` each define a local `match` that shadows base R's `match`, wrapping `IRanges::match` with `nomatch = 0`, which is easy to misread at the use site.

**Solution Summary:** Delete the two local definitions and call `IRanges::match(..., nomatch = 0)` directly at each use site, making the dispatch explicit and removing the shadow.

**Note:** Output is byte-identical, because the inlined call is exactly the local wrapper's body; this is a clarity and robustness fix, not a behaviour change.

### feature_file_editor.R

#### sRNA_calc

```
# commit 5 (previous)
# sRNA_calc: comment line 163, definition line 164, use site line 169
  ## define function to make sure match is IRanges::match and not base
  match <- function(x, table) IRanges::match(x, table, nomatch = 0)
  ...
  IGR_sRNAs <- union_peak_ranges[match(union_peak_ranges, subsetByOverlaps(union_peak_ranges, strand_IRange, maxgap = 1L)) == 0,]
```

```
# commit 6 (this commit)
# sRNA_calc: comment and definition deleted; IRanges::match inlined at the use site, line 167
  IGR_sRNAs <- union_peak_ranges[IRanges::match(union_peak_ranges, subsetByOverlaps(union_peak_ranges, strand_IRange, maxgap = 1L), nomatch = 0) == 0,]
```

#### UTR_calc

```
# commit 5 (previous)
# UTR_calc: comment line 197, definition line 198, use site line 208
  ## Define function to make sure match is IRanges::match and not base
  match <- function(x, table) IRanges::match(x, table, nomatch = 0)
  ...
  UTRs <- split_features[match(split_features, subsetByOverlaps(split_features, strand_IRange)) == 0]
```

```
# commit 6 (this commit)
# UTR_calc: comment and definition deleted; IRanges::match inlined at the use site, line 204
  UTRs <- split_features[IRanges::match(split_features, subsetByOverlaps(split_features, strand_IRange), nomatch = 0) == 0]
```

---

## Optionally expose `largestOverlap`, `fracOverlapFeature`, `read2pos`

**Issue Summary:** `count_features()` does not expose featureCounts' `largestOverlap`, `fracOverlapFeature`, and `read2pos` overlap arguments as documented parameters, so users cannot tune read-to-feature overlap behaviour.

**Solution Summary:** Add three snake_case formal arguments to `count_features()`, mapped to the camelCase featureCounts arguments, each defaulted to featureCounts' own default so existing behaviour is preserved exactly.

**Note:** Output is byte-identical with the defaults in place, because each forwarded default equals featureCounts' own default; this is an additive usability and robustness fix, not a behaviour change.

### count_features.R

#### Signature

```
# commit 5 (previous)
# lines 109-118
count_features <- function(bam_dir=".",
                           annotation_dir=".",
                           annotation_file,
                           output_dir=".",
                           output_filename="dataset",
                           chromosome_alias_file,
                           strandedness,
                           is_paired_end,
                           excl_rna = TRUE,
                           ...){
```

```
# commit 6 (this commit)
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

#### featureCounts call

```
# commit 5 (previous)
# lines 147-154
  fc <- featureCounts(bam_files,
                      annot.ext = nsaf_df,
                      chrAliases = chromosome_alias_file,
                      strandSpecific = strand_specific,
                      isPairedEnd = paired_end,
                      allowMultiOverlap = TRUE,
                      fraction = TRUE,
                      ...)
```

```
# commit 6 (this commit)
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

#### Roxygen @param entries

```
# commit 5 (previous)
# lines 100-101
#' @param excl_rna A boolean indicating if misc RNA features (rRNA, tRNA) are excluded from quantification. (Defaults=TRUE)
#' @param ... Optional parameters passed on to featureCounts()  Default: allowMultiOverlap = TRUE, fraction = TRUE
```

```
# commit 6 (this commit)
# lines 93-97
#' @param excl_rna A boolean indicating if misc RNA features (rRNA, tRNA) are excluded from quantification. (Defaults=TRUE)
#' @param largest_overlap A boolean; if TRUE, assigns each read to the feature with the largest number of overlapping bases. Maps to featureCounts largestOverlap. Combined with the package's `fraction = TRUE`, Rsubread >= 2.14.0 is recommended, since earlier versions silently miscount that combination. (Default: FALSE)
#' @param frac_overlap_feature Minimum fraction of a feature that must be overlapped before a read is assigned to it. Maps to featureCounts fracOverlapFeature. (Default: 0)
#' @param read_to_pos Reduce each read to a single base before counting: 5 for the 5' end, 3 for the 3' end, or NULL to count the whole read. Maps to featureCounts read2pos. (Default: NULL)
#' @param ... Optional parameters passed on to featureCounts()  Default: allowMultiOverlap = TRUE, fraction = TRUE
```

The three new named arguments supersede the `...` route for these settings: a caller who previously passed `largestOverlap`, `fracOverlapFeature` or `read2pos` through `...` in camelCase will now hit a loud duplicate-argument error and must switch to the snake_case argument. This is the one non-default consequence of an otherwise additive change.

---

## `feature_type` parameter is dead code

**Issue Summary:** `tpm_normalisation()` declares and documents a `feature_type` parameter that the function body never references, so a user who passes it expecting feature-type filtering is silently ignored.

**Solution Summary:** The decided resolution (Thomas, 5 June 2026) is to remove the parameter from the signature and delete its roxygen line, since it does nothing and the documentation is misleading. As a separate adjunct on the remove path, tidy the adjacent stale comment above the `read.delim(count_table)` line, which promises feature-type filtering the body has never performed.

**Note:** Under removal the output is byte-identical, because the body never reads the parameter; this is a dead-code and documentation fix. The not-taken alternative, implementing the filter, would have been a methodology change.

### tpm_norm_flagging.R

#### Roxygen @param line

```
# commit 5 (previous)
# line 7
#' @param feature_type A string indicating desired feature type(s) from annotation.
```

```
# commit 6 (this commit)
# delete the @param feature_type line; no replacement
```

#### Signature

```
# commit 5 (previous)
# line 18
tpm_normalisation <- function(count_table, complete_ann, feature_type = c("putative_sRNA", "putative_UTR"), is_gff = TRUE, output_file = NA, excl_rna = TRUE) {
```

```
# commit 6 (this commit)
# line 17
tpm_normalisation <- function(count_table, complete_ann, is_gff = TRUE, output_file = NA, excl_rna = TRUE) {
```

#### Stale comment above `read.delim(count_table)`

The comment promises feature-type filtering the body has never performed; only the comment changes, and the `read.delim` line is shown for context.

```
# commit 5 (previous)
# comment line 32 (read.delim line 33 unchanged)
  ## Load in the count table and filter for feature types
  count_df <- read.delim(count_table)
```

```
# commit 6 (this commit)
# comment line 31 (read.delim line 32)
  ## Load in the count table
  count_df <- read.delim(count_table)
```

---

## Other small style and cosmetic fixes

**Issue Summary:** Minor readability defects in `count_features.R` and `tpm_norm_flagging.R`: duplicated branch logic in `read_annotation_file`, default-quoting `write.table` calls, and a return value that contradicts the roxygen.

**Solution Summary:** Collapse `read_annotation_file` to one branch, set `quote = FALSE` on three `write.table` calls, and return the output path from `tpm_flagging`.

**Note:** Not byte-identical; the `quote = FALSE` change removes quotation marks from three output files, while the other two edits affect only an unused return value and a path string.

### count_features.R

#### `read_annotation_file`

This issue collapses the two identical branches into one and replaces `file.access(..., mode = 0) == 0` with `file.exists(...)`.

```
# commit 5 (previous)
# lines 11-22
read_annotation_file <- function(annot_dir, annot_file){
  #read in annotation file
  annot_file_loc <- c()
  if (annot_dir==".") {
    annot_file_loc <- annot_file
    stopifnot("Annotation file not found" = file.access(annot_file_loc, mode=0) == 0)
  } else {
    annot_file_loc <- paste(annot_dir, annot_file, sep = "/")
    stopifnot("Annotation file not found" = file.access(annot_file_loc, mode=0) == 0)
  }
  return(annot_file_loc)
}
```

```
# commit 6 (this commit)
# lines 11-15
read_annotation_file <- function(annot_dir, annot_file){
  annot_file_loc <- file.path(annot_dir, annot_file)
  stopifnot("Annotation file not found" = file.exists(annot_file_loc))
  return(annot_file_loc)
}
```

#### `write.table` quoting

```
# commit 5 (previous)
# line 156
  write.table(fc$counts, count_file_name, sep = "\t")
# line 158
  write.table(fc$stat, summary_file_name, sep = "\t")
```

```
# commit 6 (this commit)
# line 158
  write.table(fc$counts, count_file_name, sep = "\t", quote = FALSE)
# line 160
  write.table(fc$stat, summary_file_name, sep = "\t", quote = FALSE)
```

### tpm_norm_flagging.R

#### `write.table` quoting

```
# commit 5 (previous)
# line 54
    write.table(tpm_df, output_file, sep = "\t")
```

```
# commit 6 (this commit)
# line 53
    write.table(tpm_df, output_file, sep = "\t", quote = FALSE)
```

#### `tpm_flagging` return value

The preceding `write.table(new_annot, ...)` line is unchanged and shown for context; the edit adds the invisible return.

```
# commit 5 (previous)
# write.table line 107, closing brace line 109
  write.table(new_annot, output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

}
```

```
# commit 6 (this commit)
# lines 106-108
  write.table(new_annot, output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  invisible(output_file)
}
```

The matching `@return` roxygen is aligned in the same edit so it describes the returned path.

```
# commit 5 (previous)
# line 68
#' @return A Gff3 file, where a target feature type has an expression flag added to its attribute column.
```

```
# commit 6 (this commit)
# line 67
#' @return The path to the output GFF3 file, returned invisibly. The written file is the input annotation with an expression flag added to the attribute column of each flagged feature.
```
