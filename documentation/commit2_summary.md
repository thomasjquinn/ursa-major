# Commit 2: Syntax cleanup

## Commit Abstract: Replace a set of working but non-idiomatic R patterns across the three active source files with their standard equivalents, a no-functional-change style pass that leaves every scientific output unchanged.

### Commit Summary:
Tidy non-idiomatic R syntax, no output change

### Commit Description:
refactor(syntax): drop redundant == TRUE/== FALSE comparisons;
use message(), invisible() and isTRUE(); prefer TRUE/FALSE over T/F;
modernise paste0() and fixed = TRUE

---

## Replace `print()` with `message()` for status output

**Issue Summary:** The wrapper reports progress with `print()`, which writes to stdout, cannot be silenced, and mixes status chatter in with genuine result output.

**Solution Summary:** Replace each status `print("...")` with `message("...")`, which writes to stderr and can be turned off with `suppressMessages()`, preserving the text.

**Note:** This is a usability fix rather than a speed optimisation; data output is unchanged, though status lines move to stderr and drop the `[1]` prefix.

### feature_file_editor.R

```
# commit 1 (previous)
# line 332
print("Extracted plus strand data from BAM files")

# line 337
print("Built plus strand annotation dataframe")

# line 340
print("Extracted minus strand data from BAM files")

# line 345
print("Built minus strand annotation dataframe")

# line 361
print("Prepared complete annotation dataframe")

# line 375
print("Building output file now")
```

```
# commit 2 (this commit)
# line 332
message("Extracted plus strand data from BAM files")

# line 337
message("Built plus strand annotation dataframe")

# line 340
message("Extracted minus strand data from BAM files")

# line 345
message("Built minus strand annotation dataframe")

# line 361
message("Prepared complete annotation dataframe")

# line 375
message("Building output file now")
```

---

## Replace `return("Done!")` with `invisible()`

**Issue Summary:** Two functions end by returning the literal string `"Done!"`, a status marker that auto-prints at the console and gives callers nothing usable.

**Solution Summary:** Replace each `return("Done!")` with `invisible()`, returning the output file path where one file is written and `NULL` where two are.

**Note:** This is an idiom fix rather than a speed optimisation; the written files are unchanged, but auto-printing is suppressed and the captured return value changes.

### feature_file_editor.R

This function writes one GFF3 file to `output_file`, which is in scope at the return, so returning that path is the most useful idiom.

```
# commit 1 (previous)
# line 381
return("Done!")
```

```
# commit 2 (this commit)
# line 381
invisible(output_file)
```

### count_features.R

This function writes two files (a counts table and a summary) derived from the `output_file` base path, so no single path is the natural return value; `invisible(NULL)` is the cleaner default.

```
# commit 1 (previous)
# line 163
return("Done!")
```

```
# commit 2 (this commit)
# line 160
invisible(NULL)
```

---

## `paired_end` setup verbosity

**Issue Summary:** `count_features()` sets its `paired_end` flag with a verbose four-line default-and-overwrite block whose `== TRUE` test also errors on malformed input.

**Solution Summary:** Replace the block with a single `paired_end <- isTRUE(is_paired_end)` assignment, which states the intent directly and returns `FALSE` for unexpected inputs.

**Note:** A readability and robustness fix, not a speed optimisation; output is identical for valid scalar inputs and differs only by failing gracefully on malformed ones.

### count_features.R

The preceding comment line is left in place; only the assignment block below it changes.

```
# commit 1 (previous)
# lines 141-144
paired_end <- FALSE
if(is_paired_end==TRUE){
  paired_end <- TRUE
}
```

```
# commit 2 (this commit)
# line 141
paired_end <- isTRUE(is_paired_end)
```

---

## Redundant `== TRUE` and `== FALSE` comparisons

**Issue Summary:** Logical values are compared explicitly to `TRUE` or `FALSE` across the source files, adding visual noise and making conditions harder to read.

**Solution Summary:** Drop `== TRUE` and rewrite `== FALSE` as a leading `!`, and for the two `lapply()` sites build the logical vector with `vapply()`.

**Note:** A clarity and linting fix, not a speed optimisation; all sites are behaviour-preserving, including the two `vapply()` rewrites that replace the list-coercing comparisons.

### feature_file_editor.R

```
# commit 1 (previous)
# line 36
if (paired_end_data == FALSE & strandedness  == "stranded") {

# line 39
} else if (paired_end_data == TRUE & strandedness  == "stranded") {

# line 44
} else if (paired_end_data == FALSE & strandedness  == "reversely_stranded") {

# line 53
} else if (paired_end_data == TRUE & strandedness  == "reversely_stranded") {

# line 73
selected_peaks <- peaks[lapply(test, function(x) !is.null(x))==TRUE]

# line 137 (two comparisons on one line)
major_f <- gff[grepl("Parent", gff[,9], ignore.case = TRUE)==FALSE & gff[,3]!='chromosome' & gff[,3]!='biological_region' & grepl(ori_sRNA_biotype, gff[,9], ignore.case = TRUE)==FALSE & gff[,3]!='region' & gff[,3]!='sequence_feature',]
```

```
# commit 2 (this commit)
# line 36
if (!paired_end_data & strandedness  == "stranded") {

# line 39
} else if (paired_end_data & strandedness  == "stranded") {

# line 44
} else if (!paired_end_data & strandedness  == "reversely_stranded") {

# line 53
} else if (paired_end_data & strandedness  == "reversely_stranded") {

# line 73  -- special case: lapply() returns a list, so == TRUE was coercing
#            it to a logical vector. Build the logical vector directly.
selected_peaks <- peaks[vapply(test, function(x) !is.null(x), logical(1))]

# line 137 (both comparisons fixed)
major_f <- gff[!grepl("Parent", gff[,9], ignore.case = TRUE) & gff[,3]!='chromosome' & gff[,3]!='biological_region' & !grepl(ori_sRNA_biotype, gff[,9], ignore.case = TRUE) & gff[,3]!='region' & gff[,3]!='sequence_feature',]
```

### tpm_norm_flagging.R

```
# commit 1 (previous)
# line 20
if (is.na(output_file)==FALSE) {

# line 65
if (is.na(output_file)==FALSE) {

# line 147
if ((as.character(table_row[[3]]) %in% target_features)==TRUE) {

# line 148
if (grepl(target_flag, as.character(table_row[[9]]), ignore.case = TRUE)==TRUE){

# line 153
} else if ((as.character(table_row[[3]]) %in% target_features)==FALSE) {

# line 161
if (selection(annot_data[i,], target_features, target_flag) == TRUE){

# line 167
filtered_selection <- annot_data[unlist(filtered_vec)==TRUE,]

# line 168
selection_not_null <- filtered_selection[,lapply(filtered_selection, is.null)==FALSE]
```

```
# commit 2 (this commit)
# line 20
if (!is.na(output_file)) {

# line 65
if (!is.na(output_file)) {

# line 147
if (as.character(table_row[[3]]) %in% target_features) {

# line 148
if (grepl(target_flag, as.character(table_row[[9]]), ignore.case = TRUE)){

# line 153
} else if (!(as.character(table_row[[3]]) %in% target_features)) {

# line 161
if (selection(annot_data[i,], target_features, target_flag)){

# line 167  -- unlist() already returns an atomic logical vector, so == TRUE
#             is a pure redundancy and is simply dropped.
filtered_selection <- annot_data[unlist(filtered_vec),]

# line 168  -- special case: lapply() returns a list, so == FALSE was coercing
#             it. "!lapply(...)" would error; build the vector with vapply.
selection_not_null <- filtered_selection[,!vapply(filtered_selection, is.null, logical(1))]
```

### count_features.R

Lines 67 and 71 of the previous commit are the first line of multi-line `gff` subsetting expressions; only the `==FALSE` on each first line changes and the continuation lines are unchanged.

```
# commit 1 (previous)
# line 67
major_f <- gff[grepl("Parent", gff[,9], ignore.case = TRUE)==FALSE &

# line 71
major_f <- gff[grepl("Parent", gff[,9], ignore.case = TRUE)==FALSE &
```

```
# commit 2 (this commit)
# line 67
major_f <- gff[!grepl("Parent", gff[,9], ignore.case = TRUE) &

# line 71
major_f <- gff[!grepl("Parent", gff[,9], ignore.case = TRUE) &
```

---

## `T`/`F` shorthand instead of `TRUE`/`FALSE`

**Issue Summary:** Several sites use `T`/`F` shorthand for `TRUE`/`FALSE`, which are reassignable variable names rather than reserved words and trigger an `R CMD check` NOTE.

**Solution Summary:** Replace each `T` with `TRUE` and each `F` with `FALSE` at the affected sites, a trivial find-and-replace with no logic change.

**Note:** No functional difference, not a speed optimisation; `T` and `F` evaluate to `TRUE` and `FALSE` in a clean session, so output is identical.

### count_features.R

#### Executable sites

Lines 60, 117, 155 and 156 are plain `T`/`F` swaps. Line 66 also drops the now-redundant comparison, so it becomes `!exclude` rather than `exclude==FALSE`.

```
# commit 1 (previous)
# line 60
make_saf <- function(ann_file, exclude=F){

# line 66
  if (exclude==F){

# line 117
                           excl_rna = T,

# line 155
                      allowMultiOverlap = T,

# line 156
                      fraction = T,
```

```
# commit 2 (this commit)
# line 60
make_saf <- function(ann_file, exclude=FALSE){

# line 66
  if (!exclude){

# line 117
                           excl_rna = TRUE,

# line 152
                      allowMultiOverlap = TRUE,

# line 153
                      fraction = TRUE,
```

#### Documentation comments

The roxygen documented defaults are brought in line with the corrected code.

```
# commit 1 (previous)
# line 100
#' @param excl_rna A boolean indicating if misc RNA features (rRNA, tRNA) are excluded from quantification. (Defaults=T)

# line 101
#' @param ... Optional parameters passed on to featureCounts()  Default: allowMultiOverlap = T, fraction = T
```

```
# commit 2 (this commit)
# line 100
#' @param excl_rna A boolean indicating if misc RNA features (rRNA, tRNA) are excluded from quantification. (Defaults=TRUE)

# line 101
#' @param ... Optional parameters passed on to featureCounts()  Default: allowMultiOverlap = TRUE, fraction = TRUE
```

### tpm_norm_flagging.R

#### Executable sites

Line 18 carries two instances. Line 25 also drops the now-redundant comparison, so it becomes `is_gff` rather than `is_gff == TRUE`.

```
# commit 1 (previous)
# line 18 (two instances: is_gff = T and excl_rna = T)
tpm_normalisation <- function(count_table, complete_ann, feature_type = c("putative_sRNA", "putative_UTR"), is_gff = T, output_file = NA, excl_rna = T) {

# line 25
  if (is_gff == T){
```

```
# commit 2 (this commit)
# line 18 (both instances)
tpm_normalisation <- function(count_table, complete_ann, feature_type = c("putative_sRNA", "putative_UTR"), is_gff = TRUE, output_file = NA, excl_rna = TRUE) {

# line 25
  if (is_gff){
```

#### Documentation comments

```
# commit 1 (previous)
# line 8
#' @param is_gff A boolean indicating whether annotation is gff file, default=T

# line 9
#' @param excl_rna A boolean indicating if misc RNA features (rRNA, tRNA) are excluded from normalisation. (Defaults=T)
```

```
# commit 2 (this commit)
# line 8
#' @param is_gff A boolean indicating whether annotation is gff file, default=TRUE

# line 9
#' @param excl_rna A boolean indicating if misc RNA features (rRNA, tRNA) are excluded from normalisation. (Defaults=TRUE)
```

---

## `paste0()` and `fixed = TRUE` modernisation

**Issue Summary:** Some `paste(..., sep = "")` calls should use `paste0()`, and some `grepl()` calls on a literal pattern should pass `fixed = TRUE`.

**Solution Summary:** Replace each `paste(..., sep = "")` with `paste0()` carrying the same arguments, and add `fixed = TRUE` to the literal `grepl()` calls.

**Note:** Behaviour-preserving and chiefly an idiom and modernisation fix, with only a marginal speed gain; output is identical to the originals.

### count_features.R

```
# commit 1 (previous)
# line 134
  count_file_name <- paste(output_file, "_Counts.csv", sep = "")

# line 135
  summary_file_name <- paste(output_file, "_Count_summary.csv", sep="")
```

```
# commit 2 (this commit)
# line 134
  count_file_name <- paste0(output_file, "_Counts.csv")

# line 135
  summary_file_name <- paste0(output_file, "_Count_summary.csv")
```

### tpm_norm_flagging.R

#### paste0 modernisation

```
# commit 1 (previous)
# line 114
      new_line <- paste(i, ";expression_flag=", flags[feature_name], sep = "")
```

```
# commit 2 (this commit)
# line 114
      new_line <- paste0(i, ";expression_flag=", flags[feature_name])
```

#### grepl literal to fixed = TRUE

This line also has its redundant `==TRUE` dropped, so the updated form both adds `fixed = TRUE` and removes the comparison.

```
# commit 1 (previous)
# line 174
  while (grepl("#",f[i])==TRUE) {
```

```
# commit 2 (this commit)
# line 174
  while (grepl("#", f[i], fixed = TRUE)) {
```

### feature_file_editor.R

#### paste0 modernisation

```
# commit 1 (previous)
# line 133
    ori_sRNA_biotype <- paste("biotype=", original_sRNA_annotation, sep = "")

# line 176
    names(IGR_sRNAs) <- apply(as.data.frame(IGR_sRNAs),1, function(x) paste("ID=putative_sRNA:p", x[1], "_", x[2], ";", sep = ''))

# line 178
    names(IGR_sRNAs) <- apply(as.data.frame(IGR_sRNAs),1, function(x) paste("ID=putative_sRNA:m", x[1], "_", x[2], ";", sep = ''))

# line 217
    names(UTRs) <- apply(as.data.frame(UTRs),1, function(x) paste("ID=putative_UTR:p", x[1], "_", x[2],";", sep = ''))

# line 219
    names(UTRs) <- apply(as.data.frame(UTRs),1, function(x) paste("ID=putative_UTR:m", x[1], "_", x[2],";", sep = ''))

# line 286
        feature_attribute <- paste(cmp_strand[i,9],"upstream_feature=", previous_feature_name, ";downstream_feature=", next_feature_name, sep = "")

# line 288
        feature_attribute <- paste(cmp_strand[i,9],"upstream_feature=", next_feature_name, ";downstream_feature=", previous_feature_name, sep = "")
```

```
# commit 2 (this commit)
# line 133
    ori_sRNA_biotype <- paste0("biotype=", original_sRNA_annotation)

# line 176
    names(IGR_sRNAs) <- apply(as.data.frame(IGR_sRNAs),1, function(x) paste0("ID=putative_sRNA:p", x[1], "_", x[2], ";"))

# line 178
    names(IGR_sRNAs) <- apply(as.data.frame(IGR_sRNAs),1, function(x) paste0("ID=putative_sRNA:m", x[1], "_", x[2], ";"))

# line 217
    names(UTRs) <- apply(as.data.frame(UTRs),1, function(x) paste0("ID=putative_UTR:p", x[1], "_", x[2],";"))

# line 219
    names(UTRs) <- apply(as.data.frame(UTRs),1, function(x) paste0("ID=putative_UTR:m", x[1], "_", x[2],";"))

# line 286
        feature_attribute <- paste0(cmp_strand[i,9],"upstream_feature=", previous_feature_name, ";downstream_feature=", next_feature_name)

# line 288
        feature_attribute <- paste0(cmp_strand[i,9],"upstream_feature=", next_feature_name, ";downstream_feature=", previous_feature_name)
```

#### grepl literal to fixed = TRUE

This line also has its redundant `==TRUE` dropped, so the updated form both adds `fixed = TRUE` and removes the comparison.

```
# commit 1 (previous)
# line 367
    while (grepl("#",f[i])==TRUE) {
```

```
# commit 2 (this commit)
# line 367
    while (grepl("#", f[i], fixed = TRUE)) {
```
