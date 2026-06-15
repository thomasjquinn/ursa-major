# Commit 1: Error Handling Fixes

## Commit Abstract: Modernise the error handling by replacing the assertthat dependency with base R checks and making failures halt at the point of detection.

### Commit Summary:
Modernise error handling and drop assertthat

### Commit Description:
refactor(errors): swap assertthat for base stopifnot checks
and replace error-string returns with stop() on failure

## Drop `assertthat`, use base `stopifnot`

**Issue Summary:** The package depends on `assertthat`, a frozen third-party dependency last updated in 2019, solely for input checks that base R now provides natively.

**Solution Summary:** Replace each `assert_that()` call with base `stopifnot()` named-message syntax and remove the `assertthat` import tags and its DESCRIPTION entry.

**Note:** This is a dependency-reduction fix rather than a speed optimisation; outputs and user-facing error messages are identical, and only the internal traceback context differs.

### count_features.R

#### Roxygen import tags (lines 10, 34, 62, 112)

```
# v0-baseline (original)
# lines 10, 34, 62, 112 (identical at each location)
#' @importFrom assertthat assert_that
```

```
# commit 1 (this commit)
# delete the line at each of the four locations; no replacement
```

#### Assertion calls

```
# v0-baseline (original)
# line 17
assert_that(file.access(annot_file_loc, mode=0) == 0, msg="Annotation file not found")

# line 20
assert_that(file.access(annot_file_loc, mode=0) == 0, msg="Annotation file not found")

# lines 39-40
assert_that(strand_param %in% strand_strings,
            msg="Invalid strandedness parameter: must be either 'unstranded', 'stranded' or 'reversely_stranded'")

# line 68
assert_that(ncol(gff)==9, msg="annotation file format is invalid")

# line 131
assert_that(length(bam_files) > 0 , msg="Empty bam directory")

# line 137
assert_that(dir.exists(output_dir), msg="Output directory doesn't exist")
```

```
# commit 1 (this commit)
# line 16
stopifnot("Annotation file not found" = file.access(annot_file_loc, mode=0) == 0)

# line 19
stopifnot("Annotation file not found" = file.access(annot_file_loc, mode=0) == 0)

# line 37
stopifnot("Invalid strandedness parameter: must be either 'unstranded', 'stranded' or 'reversely_stranded'" = strand_param %in% strand_strings)

# line 64
stopifnot("annotation file format is invalid" = ncol(gff)==9)

# line 126
stopifnot("Empty bam directory" = length(bam_files) > 0)

# line 132
stopifnot("Output directory doesn't exist" = dir.exists(output_dir))
```

### tpm_norm_flagging.R

#### Roxygen import tag (line 16)

```
# v0-baseline (original)
# line 16
#' @importFrom assertthat assert_that
```

```
# commit 1 (this commit)
# delete the line; no replacement
```

#### Assertion calls

```
# v0-baseline (original)
# line 23
assert_that(dir.exists(out_dir), msg="Output directory doesn't exist.")

# line 37
assert_that(all(rownames(count_df) == feature_names), msg="Wrong feature order. Is excl_rna param set the same as for count_features?")
```

```
# commit 1 (this commit)
# line 22
stopifnot("Output directory doesn't exist." = dir.exists(out_dir))

# line 36
stopifnot("Wrong feature order. Is excl_rna param set the same as for count_features?" = all(rownames(count_df) == feature_names))
```

### DESCRIPTION

```
# v0-baseline (original)
# Imports field, line 23 within the block
Imports: 
    assertthat,
    stringr,
    tools,
    IRanges,
    GenomicAlignments,
    Rsamtools,
    Rsubread,
    DESeq2
```

```
# commit 1 (this commit)
Imports: 
    stringr,
    tools,
    IRanges,
    GenomicAlignments,
    Rsamtools,
    Rsubread,
    DESeq2
```

---

## Replace error-string returns with `stop()`

**Issue Summary:** Several functions signal failure by returning an error string instead of calling `stop()`, so failures propagate silently downstream rather than halting execution.

**Solution Summary:** Replace each error-string return with `stop()`, which raises a proper condition at the point of detection and preserves the original message text.

**Note:** This is a robustness fix, not a speed optimisation; every error message string is unchanged, and only the control flow on failure differs.

### feature_file_editor.R

#### Line 66 — single-chromosome rejection in `peak_union_calc`

The check assumes a single-chromosome genome. When the coverage vector carries more than one sequence name, the function returns an error string instead of stopping.

```
# v0-baseline (original)
# lines 63-67 (only line 66 changes; the rest is shown for context)
if (length(list_components)==1) {
  target <- list_components
} else {
  return(paste("Invalid BAM file:",f, sep = " "))
}
```

```
# commit 1 (this commit)
# lines 63-67 (only line 66 changes; the rest is shown for context)
if (length(list_components)==1) {
  target <- list_components
} else {
  stop(paste("Invalid BAM file:", f, sep = " "))
}
```

#### Lines 180, 221, 256, 290 — invalid-strand fallthrough

Each of these is the `else` branch of a `target_strand` check that fires when the strand is neither `"+"` nor `"-"`. Line 180 is in `sRNA_calc`, line 221 in `UTR_calc`, and lines 256 and 290 in `strand_feature_editor`.

```
# v0-baseline (original)
# line 180 (sRNA_calc)
return("Select strand")

# line 221 (UTR_calc)
return("Select strand")

# line 256 (strand_feature_editor)
return("Select strand")

# line 290 (strand_feature_editor)
return("Select strand")
```

```
# commit 1 (this commit)
# line 180
stop("Select strand")

# line 221
stop("Select strand")

# line 256
stop("Select strand")

# line 290
stop("Select strand")
```

#### Line 383 — empty BAM directory guard in `feature_file_editor`

This is the `else` branch reached when no BAM files are found in the directory.

```
# v0-baseline (original)
# line 383
return("No BAMs in bam directory!")
```

```
# commit 1 (this commit)
# line 383
stop("No BAMs in bam directory!")
```
