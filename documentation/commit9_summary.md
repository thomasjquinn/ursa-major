# Commit 9: Repeated GFF parsing across the pipeline

## Commit Abstract: Replace the repeated GFF reads scattered across counting, prediction and flagging with a single parse per file, held in a small reusable cache and threaded through the consumers by a path-or-cache dispatch that leaves every public signature and all six pipeline outputs unchanged.

### Commit Summary:
Parse each annotation GFF once and cache it

### Commit Description:
perf(gff): parse each annotation GFF once and reuse a cache across the
pipeline

## Repeated GFF parsing across the pipeline

**Issue Summary:** A single end-to-end pipeline run parses the same GFF file at least five times: once in each of the two `major_features()` strand calls, again where `feature_file_editor()` combines the results, in each `make_saf()` call from counting and TPM normalisation, and through the `readLines` reads in `tpm_flagging()` and `tpm_flag_filtering()`. On the small Mtb annotation the cost is minor, but on larger augmented bacterial GFFs the repeated parsing becomes a real fraction of runtime.

**Solution Summary:** Parse each distinct GFF once into a small cache holding its raw lines and its parsed dataframe, and thread that cache through the consumers via a path-or-cache dispatch, so a file path is read once and reused while a pre-built cache triggers no disk read. Public signatures are unchanged.

**Note:** All six scientific outputs stay byte-identical to the commit 8 baseline, so the cache is a pure performance change that a current baerhunter user does not need to adapt to.

### count_features.R

#### Cache helpers `load_gff_cache` and `.resolve_gff_cache`

Two new functions provide the cache and the path-or-cache dispatch. `load_gff_cache` builds the three-slot cache and is exported; `.resolve_gff_cache` returns a cache unchanged or builds one from a path, and is internal.

```
# commit 8 (previous)
# neither helper exists at commit 8
```

```
# commit 9 (this commit)
# lines 41-54 (load_gff_cache, exported)
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

# lines 72-78 (.resolve_gff_cache, internal)
.resolve_gff_cache <- function(x, annot_file_directory = ".") {
  if (is.list(x) &&
      all(c("path", "raw_lines", "parsed") %in% names(x))) {
    return(x)
  }
  load_gff_cache(x, annot_file_directory)
}
```

#### `make_saf`

The function takes the parsed dataframe from the cache instead of reading the file itself, so a pre-built cache is reused and a path is parsed once on entry.

```
# commit 8 (previous)
# lines 54-55
#function to create SAF file from gff from feature file editor
gff <- read.delim(ann_file, header = FALSE, comment.char = "#")
```

```
# commit 9 (this commit)
# lines 117-118
gff_cache <- .resolve_gff_cache(ann_file)
gff <- gff_cache$parsed
```

#### `count_features`

The function builds the cache once and passes it to `make_saf`, replacing the explicit annotation read.

```
# commit 8 (previous)
# lines 124-125 (read in annotation)
##read in annotation file (in gff3 form)
annot_file_loc <- read_annotation_file(annotation_dir, annotation_file)

# lines 131-132 (convert to SAF)
## Convert gff to SAF dataframe
nsaf_df <- make_saf(ann_file=annot_file_loc, exclude=excl_rna)
```

```
# commit 9 (this commit)
# lines 187-188 (load the cache once)
## Load the annotation once (path resolved, existence checked, parsed and raw lines cached).
gff_cache <- load_gff_cache(annotation_file, annotation_dir)

# lines 194-195 (convert to SAF from the cache)
## Convert to SAF from the cache
nsaf_df <- make_saf(ann_file = gff_cache, exclude = excl_rna)
```

### feature_file_editor.R

#### `major_features`

The function resolves its annotation argument through the cache dispatch instead of building a path and parsing the file itself.

```
# commit 8 (previous)
# lines 191-198
annot_file_loc <- c()
if (annot_file_directory==".") {
  annot_file_loc <- annotation_file
} else {
  annot_file_loc <- paste(annot_file_directory, annotation_file, sep = "/")
}

gff <- read.delim(annot_file_loc, header = FALSE, comment.char = "#")
```

```
# commit 9 (this commit)
# lines 191-192
gff_cache <- .resolve_gff_cache(annotation_file, annot_file_directory)
gff <- gff_cache$parsed
```

#### `feature_file_editor`

The wrapper loads the original GFF once at the top of the BAM-present branch and feeds both `major_features` calls, the combine step, and the header restoration from that one cache. The inline path-build-and-parse block collapses into the single cache load.

```
# commit 8 (previous)
# line 451 (plus-strand call)
maj_plus_features <- major_features(original_annotation_file, annot_file_directory = annot_file_dir, "+", original_sRNA_annotation)

# line 459 (minus-strand call)
maj_minus_features <- major_features(original_annotation_file, annot_file_directory = annot_file_dir, "-", original_sRNA_annotation)

# lines 466-474 (combine results)
annot_file_loc <- c()
if (annot_file_dir==".") {
  annot_file_loc <- original_annotation_file
} else {
  annot_file_loc <- paste(annot_file_dir, original_annotation_file, sep = "/")
}

gff <- read.delim(annot_file_loc, header = FALSE, comment.char = "#")
annotation_dataframe <- rbind(gff, plus_annot_dataframe, minus_annot_dataframe)

# line 482 (header restoration)
f <- readLines(annot_file_loc)
```

```
# commit 9 (this commit)
# line 442 (new: load the original GFF once)
gff_cache <- load_gff_cache(original_annotation_file, annot_file_dir)

# line 448 (plus-strand call)
maj_plus_features <- major_features(gff_cache, annot_file_directory = annot_file_dir, "+", original_sRNA_annotation)

# line 456 (minus-strand call)
maj_minus_features <- major_features(gff_cache, annot_file_directory = annot_file_dir, "-", original_sRNA_annotation)

# line 463 (combine results, from the cache)
annotation_dataframe <- rbind(gff_cache$parsed, plus_annot_dataframe, minus_annot_dataframe)

# line 471 (header restoration, from the cache)
f <- gff_cache$raw_lines
```

### tpm_norm_flagging.R

#### `tpm_flagging`

The function takes the raw lines from the cache instead of reading the annotation with `readLines`.

```
# commit 8 (previous)
# line 94
ann_file <- readLines(complete_annotation)
```

```
# commit 9 (this commit)
# lines 94-95
gff_cache <- .resolve_gff_cache(complete_annotation)
ann_file  <- gff_cache$raw_lines
```

#### `tpm_flag_filtering`

The function loads the flagged GFF once and takes both the parsed dataframe and the raw lines from that cache, replacing the separate `read.delim` and `readLines` reads of the same file.

```
# commit 8 (previous)
# lines 127-128 (load annotation data)
##load in annotation data.
annot_data <- read.delim(flagged_annotation_file, header = FALSE, comment.char = "#")

# lines 134-135 (restore header)
## Restore the original header.
f <- readLines(flagged_annotation_file)
```

```
# commit 9 (this commit)
# lines 128-130 (load the flagged GFF once)
## Load the flagged GFF once.
gff_cache  <- .resolve_gff_cache(flagged_annotation_file)
annot_data <- gff_cache$parsed

# lines 136-137 (restore header from the cache)
## Restore the original header from the cache.
f <- gff_cache$raw_lines
```

## Coverage model comment text in feature_file_editor.R

**Issue Summary:** The paired-end coverage model parameter carried inline and roxygen comments that were wordier than needed, and the roxygen description characterised the two models in a way that read as a mild lean rather than a neutral statement of what each does.

**Solution Summary:** Tighten all three comment sites to plain definitional wording, remove the evaluative characterisation from the roxygen, and reorder the inline branch comment so the `"fragment"` default is described first; each block keeps its original line count so no surrounding code moves.

**Note:** Comment-only and output-neutral. The executable code is computationally identical to the pre-edit file, confirmed by an identical parse tree and an identical deparsed program under R, so a current baerhunter user sees no behavioural change.

#### Roxygen `@param coverage_model` (identical in `peak_union_calc()` and the `feature_file_editor()` wrapper)

```
# commit 8 (previous)
# lines 34-41 and 425-432
#' @param coverage_model Character, one of "fragment" or "footprint", controlling
#'   how paired-end reads contribute to coverage; applies to paired-end data only
#'   and is ignored for single-end. "fragment" (the default) counts each read
#'   pair as a single fragment spanning its leftmost to rightmost aligned base,
#'   so the unsequenced insert between the mates is treated as covered.
#'   "footprint" counts only the aligned blocks of the two mates, leaving the gap
#'   between them uncovered. Fragment coverage is smoother and fuller; footprint
#'   coverage is more conservative and reflects only sequenced bases.
```

```
# commit 9 (this commit)
# lines 34-41 and 425-432
#' @param coverage_model Character, one of "fragment" or "footprint", setting how
#'   paired-end reads contribute to coverage. It applies only to paired-end data
#'   and is ignored for single-end. "fragment" (the default) treats each read
#'   pair as one span from its leftmost to rightmost aligned base, so the
#'   unsequenced insert between the mates is counted as covered. "footprint"
#'   instead counts only the aligned blocks of the two mates, so the insert
#'   between them is left uncovered and the coverage reflects only the bases
#'   that were actually sequenced.
```

#### Inline comment on the paired-end branch of `peak_union_calc()`

```
# commit 8 (previous)
# lines 119-121
      ## coverage_model (paired-end only): "footprint" counts only the aligned
      ## blocks of each mate, leaving the gap between mates uncovered; "fragment"
      ## (default) counts the whole pair from leftmost to rightmost base.
```

```
# commit 9 (this commit)
# lines 119-121 (reordered so the "fragment" default is described first)
      ## coverage_model (paired-end only): "fragment" (default) counts the whole
      ## pair from leftmost to rightmost base; "footprint" counts only the aligned
      ## blocks of each mate, leaving the gap between mates uncovered.
```
