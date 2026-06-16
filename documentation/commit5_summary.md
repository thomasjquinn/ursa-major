# Commit 5: Regex and Pattern-Matching Robustness

## Commit Abstract: Harden the package's pattern matching by making the TPM flag regex colon-agnostic, escaping the BAM file listing pattern, adding a diagnostic for unparsable GFF IDs, and vectorising the flag lookup in `tpm_flagging`.

### Commit Summary:
Harden regex and pattern matching robustness

### Commit Description:
fix(regex): match colon-less feature IDs in tpm_flagging,
escape the .BAM listing pattern, warn on unparsable GFF IDs,
and vectorise the flag lookup

**Note:** Unlike the earlier commits, this one is not entirely byte-identical on the project data: on the NCBI RefSeq H37Rv annotation used by the Cortes run, the colon-agnostic regex in the first section below newly flags the gene set in `flagged.gff3` and `filtered.gff3` (roughly 4060 features), while every other output stays byte-identical. On a colon-style (Ensembl) annotation the whole commit is byte-identical, and the other two changes here (the `.BAM$` escape and the `strand_feature_editor()` diagnostic) add only warnings.

## TPM-to-GFF matching regex requires colon in feature ID

**Issue Summary:** The ID-extraction regex in `tpm_flagging()` requires a colon, so annotations with bare colon-less IDs are silently passed through unflagged with no warning.

**Solution Summary:** Make the capture colon-agnostic to match how `make_saf()` builds the TPM-table IDs, and add a guard that warns when no feature ID matches the table.

**Note:** Output is byte-identical only on colon-style (Ensembl) annotations. The project's Cortes run uses the NCBI RefSeq H37Rv annotation, whose colon-less IDs (`ID=gene-Rv0001;`) the colon-requiring regex never matched, so on the project data this is a documented correctness change: it flags roughly 4060 gene-level features (4060 of 10848 lines in `flagged.gff3`, 4060 of 8448 in `filtered.gff3`) that the baseline and Commit 4 silently left unflagged, while the four non-flagging outputs stay byte-identical. The only other added behaviour is a no-match warning, which does not fire on Cortes because the gene IDs now match.

### tpm_norm_flagging.R

This change is not committed as a standalone block. In this same commit the per-line loop the colon-agnostic capture would have restructured is dissolved by the flag-lookup vectorisation, so the colon-agnostic `sub()` ships hoisted to a single vectorised call and the no-match guard ships rewritten against the `matched` vector. The committed form of the whole region is shown under "Vectorise the flag lookup in `tpm_flagging`" at lines 97-105 of the shipped file, where the first line, `feature_names <- sub(".*?ID=(.*?);.*", "\\1", ann_file)`, is this issue's colon-agnostic, vectorised capture. Against Commit 4 the contribution here is the pattern, from the colon-requiring `.*?ID=(.*?:.*?);.*` to the colon-agnostic `.*?ID=(.*?);.*`, together with the added no-match warning.

---

## Vectorise the flag lookup in `tpm_flagging`

**Issue Summary:** Commit 4 pre-allocated `new_annot` and removed the O(N²) vector growth, but the per-line `feature_name %in% flag_names` membership test inside the loop was left unchanged. That test is O(lines x features), so `tpm_flagging` stays super-linear. The 1 June single-end scaling run measured it at baseline 121.1 s against Commit-4 78.6 s at 50,000 annotation lines, about 1.5 times faster but not flat.

**Solution Summary:** Dissolve the loop. With `feature_names` already extracted in one vectorised `sub()` call (hoisted by the colon-agnostic regex change in the same commit), compute membership once over the whole vector and assign the flagged lines by logical index, replacing the per-line loop with three whole-vector operations.

**Note:** The vectorisation itself is output-neutral: replacing the per-line loop with whole-vector operations produces an identical flagged GFF in isolation. It sits in the same commit as the colon-agnostic regex, however, so the combined Commit 5 effect on `flagged.gff3` and `filtered.gff3` is not byte-identical on the project's RefSeq annotation, where the colon-agnostic regex newly flags the gene set. This is a pure performance fix that completes the linearisation of `tpm_flagging` begun by the earlier pre-allocation and the colon-agnostic regex; it carries no methodology change of its own.

### tpm_norm_flagging.R

This is the committed form of the region and carries both of this commit's `tpm_flagging` changes: the colon-agnostic capture and no-match guard (the regex correctness change) and the loop dissolution (this performance change). The previous-commit block is Commit 4's colon-requiring, pre-allocated, indexed loop; the this-commit block replaces it with the vectorised extraction, a single `matched` membership test, the no-match guard, an unchanged-line default, and an indexed assignment.

```
# commit 4 (previous)
# lines 97-106
  new_annot <- character(length(ann_file))
  for (i in seq_along(ann_file)) {
    feature_name <- sub(".*?ID=(.*?:.*?);.*", "\\1", ann_file[i])

    if (feature_name %in% flag_names) {
      new_annot[i] <- paste0(ann_file[i], ";expression_flag=", flags[feature_name])
    } else {
      new_annot[i] <- ann_file[i]
    }
  }
```

```
# commit 5 (this commit)
# lines 97-105
  feature_names <- sub(".*?ID=(.*?);.*", "\\1", ann_file)
  matched <- feature_names %in% flag_names
  if (!any(matched)) {
    warning("No annotation feature IDs matched the TPM table; check that the GFF ID format is consistent with the count table.",
            call. = FALSE, immediate. = TRUE)
  }
  new_annot <- ann_file
  new_annot[matched] <- paste0(ann_file[matched], ";expression_flag=",
                               flags[feature_names[matched]])
```

---

## BAM file regex matches more than intended

**Issue Summary:** Three `list.files()` calls use `pattern = ".BAM$"`, where the unescaped dot is a regex wildcard, so files such as `sampleBAM` are matched as BAM input.

**Solution Summary:** Escape the dot to `pattern = "\\.BAM$"` at all three sites so only a literal `.BAM` ending matches, keeping `ignore.case = TRUE`.

**Note:** Output is byte-identical for any dotted BAM filenames, including the Cortes data; this is a robustness fix and no test is required.

### feature_file_editor.R

```
# commit 4 (previous)
# line 29
    bam_files <- list.files(path = bam_location, pattern = ".BAM$", full.names = TRUE, ignore.case = TRUE)
# line 325
  test <- list.files(path = bam_directory, pattern = ".BAM$", full.names = TRUE, ignore.case = TRUE)
```

```
# commit 5 (this commit)
# line 29
    bam_files <- list.files(path = bam_location, pattern = "\\.BAM$", full.names = TRUE, ignore.case = TRUE)
# line 343
  test <- list.files(path = bam_directory, pattern = "\\.BAM$", full.names = TRUE, ignore.case = TRUE)
```

### count_features.R

```
# commit 4 (previous)
# line 125
  bam_files <- list.files(path = bam_dir, pattern = ".BAM$", full.names = TRUE, ignore.case = TRUE)
```

```
# commit 5 (this commit)
# line 125
  bam_files <- list.files(path = bam_dir, pattern = "\\.BAM$", full.names = TRUE, ignore.case = TRUE)
```

---

## ID-parsing regex falls through silently

**Issue Summary:** `strand_feature_editor()` parses GFF IDs with colon-requiring `sub()` patterns that return the input unchanged on no-match, so malformed attributes are carried downstream silently.

**Solution Summary:** Wrap each parse in a helper that runs the same `sub()` and records the offending attribute when the pattern does not match, leaving matches untouched, then emit a single deduplicated summary warning at the end of the function reporting how many feature IDs could not be parsed.

**Note:** Output is byte-identical on the Cortes run; this is a robustness fix that only adds a diagnostic, a single end-of-function summary warning, on the previously silent no-match case.

### feature_file_editor.R

#### `parse_id` helper and accumulator

An accumulator and a local helper are added near the start of the function. The helper runs the same `sub()` and records the attribute on no-match, writing to the accumulator in the enclosing function.

```
# commit 4 (previous)
# no commit-4 equivalent; new code inserted at the head of strand_feature_editor(), after commit-4 line 237
```

```
# commit 5 (this commit)
# lines 238-247
  ## Collect GFF attributes whose ID could not be parsed, to report once after the loop.
  unparsed_attrs <- character(0)
  ## Parse a capture group from a GFF ID attribute, recording the attribute on no-match.
  parse_id <- function(attr, pattern) {
    parsed <- sub(pattern, "\\1", attr)
    if (parsed == attr) {
      unparsed_attrs <<- c(unparsed_attrs, attr)
    }
    parsed
  }
```

#### ID-parsing sites

Each parse keeps its original pattern; only the call is wrapped in `parse_id()`, so the parsed value is unchanged for every matching attribute.

```
# commit 4 (previous)
# line 263
  previous_feature_name <- sub("ID=.*?:(.*?);.*", "\\1", cmp_strand[nrow(cmp_strand),9])
# line 267
    feature_name <- sub("ID=.*?:(.*?);.*", "\\1", cmp_strand[i,9])
# line 269 (type pattern: captures before the colon)
      feature_type <- sub("ID=(.*?):.*?;.*", "\\1", cmp_strand[i,9])
# line 274
        next_feature_name <- sub("ID=.*?:(.*?);.*", "\\1", cmp_strand[i+1,9])
# line 276
        next_feature_name <- sub("ID=.*?:(.*?);.*", "\\1", cmp_strand[1,9])
```

```
# commit 5 (this commit)
# line 273
  previous_feature_name <- parse_id(cmp_strand[nrow(cmp_strand),9], "ID=.*?:(.*?);.*")
# line 277
    feature_name <- parse_id(cmp_strand[i,9], "ID=.*?:(.*?);.*")
# line 279
      feature_type <- parse_id(cmp_strand[i,9], "ID=(.*?):.*?;.*")
# line 284
        next_feature_name <- parse_id(cmp_strand[i+1,9], "ID=.*?:(.*?);.*")
# line 286
        next_feature_name <- parse_id(cmp_strand[1,9], "ID=.*?:(.*?);.*")
```

#### End-of-function summary warning

A single deduplicated summary warning is added after the loop, immediately before `return(cmp_strand)`. It fires only when at least one attribute failed to parse.

```
# commit 4 (previous)
# no commit-4 equivalent; new code inserted before return(cmp_strand), after commit-4 line 294
```

```
# commit 5 (this commit)
# lines 305-312
  ## Report unparsable feature IDs once, as a single deduplicated summary.
  if (length(unparsed_attrs) > 0L) {
    n_failed <- length(unique(unparsed_attrs))
    warning(paste0(n_failed, " of ", nrow(cmp_strand),
                   " feature IDs could not be parsed from the GFF attribute column (e.g. ",
                   unparsed_attrs[1], "); IDs should have the form ID=type:name;."),
            call. = FALSE, immediate. = TRUE)
  }
```
