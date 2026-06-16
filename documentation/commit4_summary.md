# Commit 4: TPM-filtering cluster updates

## Commit Abstract: Vectorise the TPM normalisation and flag-filtering code by replacing the per-column dataframe growth, the RPK summation, the flagging accumulation loop, and the per-row filter predicate with whole-object operations, simplify the filter write path, and switch the flag match to a literal.

### Commit Summary:
Vectorise TPM normalisation and flag filtering

### Commit Description:
perf(tpm): vectorise the normalisation and flag-filtering paths
and match target_flag as a literal

## Column-by-column dataframe growth in `tpm_normalisation`

**Issue Summary:** The RPK and TPM dataframes are built one column at a time across two loops, copying the growing object every iteration, for O(N²) cost in sample count.

**Solution Summary:** Vectorise both loops: `rpk_df <- count_df / feature_lengths` recycles down each column, and `sweep(rpk_df, 2, scaling_fact, "/")` produces TPM in one call.

**Note:** Output is byte-identical; this is a performance fix that scales the function linearly with sample count and is independent of the other Commit 4 edits.

### tpm_norm_flagging.R

The two accumulation loops in `tpm_normalisation()` each collapse to a single vectorised line. The scaling-factor lines between them and the name-restoration lines after them are unchanged.

```
# commit 3 (previous)
# lines 41-47 (RPK)
  rpk_df <- data.frame( count_df[,1] / feature_lengths )
  if (ncol(count_df) > 1){
    for (i in 2:ncol(count_df)) {
      sample_rpk <- count_df[ , i] / feature_lengths
      rpk_df     <- data.frame(rpk_df, sample_rpk)
    }
  }

# lines 53-59 (TPM)
  tpm_df <- data.frame(rpk_df[,1]/scaling_fact[1])
  if (ncol(rpk_df) > 1){
    for (n in 2:ncol(rpk_df)) {
      sample_tpm <- rpk_df[,n]/scaling_fact[n]
      tpm_df <- data.frame(tpm_df, sample_tpm)
    }
  }
```

```
# commit 4 (this commit)
# RPK: divide the dataframe by the length vector, which recycles down each column
  rpk_df <- count_df / feature_lengths

# TPM: divide each column by its per-sample scaling factor
  tpm_df <- sweep(rpk_df, 2, scaling_fact, "/")
```

---

## `apply(rpk_df, 2, sum)` instead of `colSums`

**Issue Summary:** Per-sample RPK totals are computed with `apply(rpk_df, 2, sum)`, which dispatches an R-level call per column instead of using a dedicated compiled primitive.

**Solution Summary:** Replace `apply(rpk_df, 2, sum)` with `colSums(rpk_df)`, a drop-in substitution that returns the same named numeric vector via compiled C.

**Note:** Output is byte-identical; this is a trivial drop-in performance fix with no separate test required (covered incidentally by the Commit 4 regression run).

### tpm_norm_flagging.R

```
# commit 3 (previous)
# line 49
  sample_rpk_sum <- apply(rpk_df, 2, sum)
```

```
# commit 4 (this commit)
  sample_rpk_sum <- colSums(rpk_df)
```

---

## Vector-growing loop in `tpm_flagging`

**Issue Summary:** `new_annot` is initialised empty and grown one element at a time inside the annotation loop, copying the growing vector every iteration for O(N²) cost in annotation-line count.

**Solution Summary:** Pre-allocate `new_annot` to `character(length(ann_file))` and assign into it by index, switching the loop from `for (i in ann_file)` to `for (i in seq_along(ann_file))`.

**Note:** Output is byte-identical; this is a performance fix in `tpm_flagging` that is independent of the other Commit 4 edits and scales linearly with annotation-line count.

### tpm_norm_flagging.R

The empty-then-append accumulation in `tpm_flagging()` is replaced by a pre-allocated vector with indexed assignment. The previous commit already uses `paste0` on the assignment line; this issue's own change is the accumulation strategy, not the paste form.

```
# commit 3 (previous)
# lines 109-119
  new_annot <- c()
  for (i in ann_file) {
    feature_name <- sub(".*?ID=(.*?:.*?);.*", "\\1", i)

    if (feature_name %in% flag_names) {
      new_line <- paste0(i, ";expression_flag=", flags[feature_name])
      new_annot <- c(new_annot, new_line)
    } else {
      new_annot <- c(new_annot, i)
    }
  }
```

```
# commit 4 (this commit): pre-allocate and assign by index
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

---

## Vectorise the row filter in `tpm_flag_filtering`

**Issue Summary:** The filter iterates over every row of `annot_data` to call a predicate built from vectorised primitives, so per-row R-level dispatch dominates as feature count grows.

**Solution Summary:** Replace the loop and the `selection` helper with three whole-column operations (`%in%`, `grepl`, logical combination) and one logical subset.

**Note:** Output is byte-identical; this is the structural performance fix that owns the Commit 4 rewrite of the block and absorbs the local pre-allocation fix recorded separately below.

### tpm_norm_flagging.R

The internal `selection` helper and the per-row loop are replaced by three whole-column operations and a single logical subset. The flag-match line ships with `fixed = TRUE`; that argument change is the one behaviour change in this commit, documented separately below.

```
# commit 3 (previous)
# lines 145-169
  # An internal function to go examine a table row: all target features are checked and filtered by the desired flag; all the other features are kept.
  selection <- function(table_row, target_features, target_flag) {
    if (as.character(table_row[[3]]) %in% target_features) {
      if (grepl(target_flag, as.character(table_row[[9]]), ignore.case = TRUE)){
        return(TRUE)
      }else{
        return(FALSE)
      }
    } else if (!(as.character(table_row[[3]]) %in% target_features)) {
      return(TRUE)
    }
  }
  #apply function across all rows
  # (error using apply with no null rows, this maintains matrix structure)
  filtered_vec <- list(length(annot_data))
  for (i in 1:nrow(annot_data)){
    if (selection(annot_data[i,], target_features, target_flag)){
      filtered_vec[i] <- TRUE
    }else{
      filtered_vec[i] <- FALSE
    }
  }
  filtered_selection <- annot_data[unlist(filtered_vec),]
  selection_not_null <- filtered_selection[,!vapply(filtered_selection, is.null, logical(1))]
  df <- data.frame(matrix(unlist(selection_not_null), nrow=nrow(selection_not_null), byrow=F))
```

```
# commit 4 (this commit): three whole-column operations and one logical subset
  is_target  <- annot_data[, 3] %in% target_features
  flag_match <- grepl(target_flag, annot_data[, 9], fixed = TRUE)
  keep       <- !is_target | flag_match
  filtered_selection <- annot_data[keep, ]
```

---

## Simplify the write path in `tpm_flag_filtering`

**Issue Summary:** A no-op column-null filter is chained to a destructive matrix-unlist round-trip; the dataframe is corrupted in memory and `write.table` then hides the damage.

**Solution Summary:** Delete both lines once the vectorised refactor produces a clean `filtered_selection`, and write that dataframe directly.

**Note:** Output GFF3 is byte-identical; this is a robustness fix that simplifies the output handoff, co-applied after the structural rewrite in Commit 4.

### tpm_norm_flagging.R

With `filtered_selection` already a clean dataframe, the no-op column-null filter and the matrix-unlist round-trip are removed, and the `write.table` consumer is renamed from `df` to `filtered_selection`.

```
# commit 3 (previous)
# lines 168-169
  selection_not_null <- filtered_selection[,!vapply(filtered_selection, is.null, logical(1))]
  df <- data.frame(matrix(unlist(selection_not_null), nrow=nrow(selection_not_null), byrow=F))

# line 180
  write.table(df, output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
```

```
# commit 4 (this commit)
# the two lines above are deleted entirely (no replacement); filtered_selection is already a clean dataframe

# write filtered_selection directly
  write.table(filtered_selection, output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
```

---

## `target_flag` matched as regex via `grepl`

**Issue Summary:** `target_flag` is matched as a case-insensitive regex, so metacharacters in the flag are interpreted as a pattern rather than as literal characters.

**Solution Summary:** Set `fixed = TRUE` and remove `ignore.case = TRUE` so the flag matches literally and case-sensitively, applied to the `grepl` line that the vectorised rewrite introduces.

**Note:** Output is unchanged for typical usage (the standard flag names are lowercase literals); this is the one behaviour-changing item in Commit 4, affecting only case-mismatched or metacharacter-containing flags.

### tpm_norm_flagging.R

In the previous commit the flag was matched by a case-insensitive regex inside the per-row `selection` helper. The vectorised rewrite lifts that match to a whole-column operation, and this issue changes its argument from `ignore.case = TRUE` to `fixed = TRUE`, so the match becomes a literal, case-sensitive substring test. It is the only behaviour change in the commit.

```
# commit 3 (previous)
# line 148 (inside the per-row selection helper)
      if (grepl(target_flag, as.character(table_row[[9]]), ignore.case = TRUE)){
```

```
# commit 4 (this commit): literal, case-sensitive match on the vectorised line
  flag_match <- grepl(target_flag, annot_data[, 9], fixed = TRUE)
```

---

## Broken pre-allocation in `tpm_flag_filtering`

**Issue Summary:** The pre-allocation `filtered_vec <- list(length(annot_data))` is doubly broken: `list(x)` makes a one-element list, and `length()` on a dataframe gives the column count.

**Solution Summary:** Replace the broken pre-allocation and its `for` loop with one `vapply` call that returns a typed logical vector of the right length.

**Note:** Output is byte-identical; this is a performance fix that is superseded inside Commit 4 by the fuller refactor of the same block.

This change makes no independent edit. The broken pre-allocation and the `for` loop that follows are removed entirely by the vectorised rewrite shown above, where the `selection` helper and row loop collapse into the whole-column `keep` vector. The standalone `vapply` fix this issue describes therefore does not land in the shipped code; it is recorded here for audit. The broken line itself is visible in the previous-commit block of the rewrite above (commit 3, line 159).
