# Commit 14: Parameter Scout Merge and Source Annotation Audit

## Commit Abstract: Merge the two parameter scout modules into one packageable file, and correct the annotation of the three production source files, without changing one executable line anywhere.

### Commit Summary:
Merge parameter scout and audit production sources

### Commit Description:

```
refactor(parameter-scout): merge the pair and make it packageable

Merge parameter_scout_paired_end.R and parameter_scout_single_end.R
into one parameter_scout.R exporting parameter_scout(). The library
type becomes the paired_end_data argument, defaulting to FALSE to
match feature_file_editor(). Delete the top-level attach block and
declare the imports as roxygen instead; qualify read.delim; remove
the module's own version stamp. Reduce the annotation from 768 lines
to 405.

Audit feature_file_editor.R, count_features.R and tpm_norm_flagging.R
to the same annotation standard and remove the commit markers.

No computed result changes. Every definition deparses identically to
its Commit 11 version but for the four lines the merge and the
qualification alter.
```

## Merge the parameter scout into one packageable file

**Issue Summary:** The parameter scout ships as two modules, one for paired-end and one for single-end libraries, 99.2% identical at 1,767 of 1,782 lines. Both carry a top-level attach block, which a package must not have, and both are annotated as working files rather than as source.

**Solution Summary:** Merge the pair into one `parameter_scout.R` exporting `parameter_scout()`, with the library type becoming a `paired_end_data` argument defaulting to `FALSE`. Delete the attach block and declare the imports as roxygen, qualify `read.delim`, and reduce the annotation to the standard the source files are audited to below.

**Note:** A packaging change, not a behaviour change. Every definition appears once and deparses identically to its Commit 11 version but for the four lines the merge and the qualification alter.

### parameter_scout.R

`parameter_scout_paired_end.R` and `parameter_scout_single_end.R` are deleted and replaced by one file. The before column is one module, the two being near-identical.

| measure | one module before | merged file after |
|---|---|---|
| lines | 1,782 | 1,403 |
| executable lines | 854 | 842 |
| annotation lines | 768 | 405 |

Executable lines fall by twelve: the attach block takes fourteen with it, and `read.delim` gains two through wrapping at its two call sites. `documentation/parameter_scout_instructions.md` is rewritten to match.

---

## Remove the module's own version stamp

**Issue Summary:** The scout carries a version string of its own in its header block, which would report a number unrelated to the package version once the module is part of the package.

**Solution Summary:** Remove the version rather than align it, so the package version is the only version.

**Note:** Done in the same pass as the merge, since the stamp sits in the header block the merge rewrites.

---

## Documentation true-up

**Issue Summary:** The commit summaries, the integration notes and the scout instructions all sit together in `documentation/`, so a reader cannot tell the project record from the user documentation without opening a file. Two of the files in there describe a repository that no longer exists.

**Solution Summary:** Create `commit_notes/` for the project record and move every commit summary into it, leaving `documentation/` for the scout instructions alone. Delete the files that are out of date.

**Note:** Documentary only. No source file is touched by this and nothing in the package reads any of it.

### commit_notes/

The eleven existing summaries, `commit1_summary.md` through `commit11_summary.md`, are moved from `documentation/` into the new folder. Three are added with them:

* `commit12_summary.md` — the repository configuration commit
* `commit13_summary.md` — the bundled example data commit
* `commit14_summary.md` — this commit

### Deletions

Three items are removed from `documentation/`, all describing the repository as it was before it became an installable package:

* `integration.md`
* `notes.md`
* the `integration/` folder

---

## Correct the annotation of the three source files

**Issue Summary:** The three source files carry eleven commits' worth of accreted annotation and a block of `#commitN completed` markers at the foot of each. Several comments name commits or code the reader cannot see, and several are wrong against the code beside them.

**Solution Summary:** Correct each comment that is wrong, opinionated or restates the code, to one standard: a comment describes the code as it currently is, with no commits, no issue numbers, no dates and no comparison against code the reader cannot see. Then remove the marker block from the foot of each file.

**Note:** Comment lines and trailing marker lines only. Not one executable line, in any of the three files.

### count_features.R

#### Lines 1 to 2, the file header block. Names a parameter the function does not have.

```
# commit 13 (previous)
## Structural-RNA vocabulary for make_saf(): the feature types treated as tRNA
## or rRNA when excl_rna = TRUE. Transcript-level types, gene-level types and
```

```
# commit 14 (this commit)
## Structural-RNA vocabulary for make_saf(): the feature types treated as tRNA
## or rRNA when it is called with exclude = TRUE. Transcript-level types, gene-level types and
```

#### Line 11, `read_annotation_file` description. Wrong against the code.

```
# commit 13 (previous)
#' This function pastes together the path and filename of annotation file (if not in current directory) and tests file for existence.
```

```
# commit 14 (this commit)
#' This function joins the directory and filename of the annotation file and tests the result for existence.
```

#### Line 93, `@param strand_param`. Strandedness keyword, first site. Wrong separator, incomplete set.

```
# commit 13 (previous)
#' @param strand_param user input 'stranded' or 'reversely-stranded'
```

```
# commit 14 (this commit)
#' @param strand_param user input 'unstranded', 'stranded' or 'reversely_stranded'
```

#### Line 117, `@param ann_file`. Understates what the argument accepts.

```
# commit 13 (previous)
#' @param ann_file A GFF3 annotation file.
```

```
# commit 14 (this commit)
#' @param ann_file Either a GFF3 annotation file or a pre-built GFF cache (see \code{load_gff_cache}).
```

#### Line 118, `@param exclude`. Sense inverted against the argument name.

```
# commit 13 (previous)
#' @param exclude A boolean to indicate whether or not to include rRNA/tRNA features
```

```
# commit 14 (this commit)
#' @param exclude A boolean to indicate whether or not to exclude rRNA/tRNA features
```

#### After line 120, `make_saf()` roxygen. Filter policy asymmetry. New text, additive.

```
# commit 14 (this commit)
#' @note By design this quantification filter excludes only tRNA and rRNA when
#'   \code{exclude = TRUE}, keeping all other annotated RNAs so they are counted.
#'   This is deliberately more permissive than \code{major_features()}
#'   (prediction), which excludes all pre-annotated ncRNAs to avoid
#'   re-discovering them.
```

#### Line 164, `@param output_dir`. Wrong format.

```
# commit 13 (previous)
#' @param output_dir The full directory path for CSV output files to be written
```

```
# commit 14 (this commit)
#' @param output_dir The full directory path for the output files to be written. The files carry a .csv suffix but are tab-separated.
```

#### Line 167, `@param strandedness`. Strandedness keyword, second site. Wrong separator.

```
# commit 13 (previous)
#' @param strandedness A string outlining the type of the sequencing library: unstranded, stranded, or reversely stranded.
```

```
# commit 14 (this commit)
#' @param strandedness A string outlining the type of the sequencing library: unstranded, stranded, or reversely_stranded.
```

#### Line 170, `@param largest_overlap`. Markup inconsistent with the file.

```
# commit 13 (previous)
#' @param largest_overlap A boolean; if TRUE, assigns each read to the feature with the largest number of overlapping bases. Maps to featureCounts largestOverlap. Combined with the package's `fraction = TRUE`, Rsubread >= 2.14.0 is recommended, since earlier versions silently miscount that combination. (Default: FALSE)
```

```
# commit 14 (this commit)
#' @param largest_overlap A boolean; if TRUE, assigns each read to the feature with the largest number of overlapping bases. Maps to featureCounts largestOverlap. Combined with the package's \code{fraction = TRUE}, Rsubread >= 2.14.0 is recommended, since earlier versions silently miscount that combination. (Default: FALSE)
```

#### Line 173, `@param count_multi_mapping_reads`. Compares against code the reader cannot see.

```
# commit 13 (previous)
#' @param count_multi_mapping_reads A boolean; if FALSE, reads mapping to multiple locations are excluded from counts. Maps to featureCounts countMultiMappingReads. (Default: FALSE, reproducing 2019 behaviour)
```

```
# commit 14 (this commit)
#' @param count_multi_mapping_reads A boolean; if FALSE, reads mapping to multiple locations are excluded from counts. Maps to featureCounts countMultiMappingReads. (Default: FALSE)
```

#### Line 177, `@return` of `count_features()`. Silent on what is returned.

```
# commit 13 (previous)
#' @return Count tables for each feature are written into separate files, as well as the result summary.
```

```
# commit 14 (this commit)
#' @return NULL, invisibly. Count tables for each feature are written into separate files, as well as the result summary.
```

#### Lines 231 to 232, inline in the `featureCounts` call. Names a commit.

```
# commit 13 (previous)
                      # countReadPairs requires Rsubread >= 2.4.3 (present in the 2.4.3 reference manual, RELEASE_3_12, 30 March 2021);
                      # forwarded unconditionally, with the matching DESCRIPTION Imports floor added at the documentation commit.
```

```
# commit 14 (this commit)
                      # countReadPairs requires Rsubread >= 2.4.3.
```

#### Lines 249 to 258, the commit markers. Taken last within the file.

---

### feature_file_editor.R

#### Line 4, `peak_union_calc()` description. Contradicts its own `@return`.

```
# commit 13 (previous)
#' The function goes over each BAM file in the directory and finds the expression peaks that satisfy the coverage boundary and length criteria in each file. Then it unifies the peak information to obtain a single set of peak genomic coordinates.
```

```
# commit 14 (this commit)
#' The function goes over each BAM file in the directory and finds the expression peaks that satisfy the coverage boundary and length criteria in each file. Then it unifies the peak information to obtain one set of peak genomic coordinates for each strand.
```

#### Lines 8 to 10, the three threshold `@param` entries of `peak_union_calc()`. Threshold inclusivity.

```
# commit 13 (previous)
#' @param low_coverage_cutoff An integer indicating the low coverage threshold value.
#' @param high_coverage_cutoff An integer indicating the high coverage threshold value.
#' @param peak_width An integer indicating the minimum peak width.
```

```
# commit 14 (this commit)
#' @param low_coverage_cutoff An integer indicating the low coverage threshold value. Inclusive: coverage equal to this value qualifies.
#' @param high_coverage_cutoff An integer indicating the high coverage threshold value. Exclusive: coverage must be greater than this value.
#' @param peak_width An integer indicating the minimum peak width. Exclusive: a peak must carry more than this many positions above the high cut-off.
```

#### Line 12, `@param strandedness`. Strandedness keyword, third site. Commit 10 touched this line.

```
# commit 13 (previous)
#' @param strandedness A string outlining the type of the sequencing library: stranded, or reversely stranded. Defaults to "stranded"; "unstranded" is rejected with an error.
```

```
# commit 14 (this commit)
#' @param strandedness A string outlining the type of the sequencing library: stranded, or reversely_stranded. Defaults to "stranded"; "unstranded" is rejected with an error.
```

#### Lines 43 to 44, `@return` of `peak_union_calc()`. Markup inconsistent with the file.

```
# commit 13 (previous)
#' @return A named list with two IRanges objects, `plus` and `minus`, holding the
#'   unified peak coordinates for each strand.
```

```
# commit 14 (this commit)
#' @return A named list with two IRanges objects, \code{plus} and \code{minus}, holding the
#'   unified peak coordinates for each strand.
```

#### Lines 82 to 83, in `compute_strand_peaks()`. Compares against code the reader cannot see.

```
# commit 13 (previous)
    ## Count, per peak, the positions above the high cut-off, then keep the peaks
    ## whose count exceeds the minimum width. This reproduces the count-based test.
```

```
# commit 14 (this commit)
    ## Count, per peak, the positions above the high cut-off, then keep the peaks
    ## whose count exceeds the minimum width.
```

#### Line 166, `peak_analysis()` description. Describes a test the code does not perform.

```
# commit 13 (previous)
#' This is a helper function that is used to examine if the peak had a continuous stretch of a given width that has coverage above the high cut-off value.
```

```
# commit 14 (this commit)
#' This is a helper function that is used to examine whether the number of positions with coverage above the high cut-off value exceeds the given width. The positions need not be contiguous; contiguity is imposed by the low cut-off slice that produces the RleViews line.
```

#### Lines 169 to 170, `peak_analysis()`'s two threshold `@param` entries. Threshold inclusivity.

```
# commit 13 (previous)
#' @param high_cutoff An integer indicating the high coverage threshold value.
#' @param min_sRNA_length An integer indicating the minimum sRNA length (peak width).
```

```
# commit 14 (this commit)
#' @param high_cutoff An integer indicating the high coverage threshold value. Exclusive: coverage must be greater than this value.
#' @param min_sRNA_length An integer indicating the minimum sRNA length (peak width). Exclusive: the peak must carry more than this many positions above the high cut-off.
```

#### Line 177, inline in `peak_analysis()`. Verbatim restatement of the roxygen. Deleted, not rewritten.

```
# commit 13 (previous)
  ## This is a helper function that is used to examine if the peak had a continuous stretch of a given width that has coverage above the high cut-off value.
```

```
# commit 14 (this commit)
  (deleted, no replacement)
```

#### Line 190, `@param annotation_file` of `major_features()`. Understates what the argument accepts.

```
# commit 13 (previous)
#' @param annotation_file  GFF3 genome annotation file.
```

```
# commit 14 (this commit)
#' @param annotation_file  Either a GFF3 genome annotation file or a pre-built GFF cache (see \code{load_gff_cache}).
```

#### After line 195, `major_features()` roxygen. Filter policy asymmetry. New text, additive.

```
# commit 14 (this commit)
#' @note By design this prediction filter excludes pre-annotated ncRNAs so they
#'   are not re-discovered as novel features, tRNA and rRNA being retained for
#'   masking. This is deliberately more aggressive than \code{make_saf()}
#'   (quantification), which excludes only tRNA and rRNA so other annotated RNAs
#'   are still counted.
```

#### Lines 432 to 434, the three threshold `@param` entries of `feature_file_editor()`. Threshold inclusivity.

```
# commit 13 (previous)
#' @param low_coverage_cutoff An integer indicating the low coverage threshold value.
#' @param high_coverage_cutoff An integer indicating the high coverage threshold value.
#' @param min_sRNA_length An integer indicating the minimum peak width/sRNA length.
```

```
# commit 14 (this commit)
#' @param low_coverage_cutoff An integer indicating the low coverage threshold value. Inclusive: coverage equal to this value qualifies.
#' @param high_coverage_cutoff An integer indicating the high coverage threshold value. Exclusive: coverage must be greater than this value.
#' @param min_sRNA_length An integer indicating the minimum peak width/sRNA length. Exclusive: a peak must carry more than this many positions above the high cut-off.
```

#### Line 437, `@param strandedness`. Strandedness keyword, fourth site. Commit 10 touched this line.

```
# commit 13 (previous)
#' @param strandedness A string outlining the type of the sequencing library: stranded, or reversely stranded. Defaults to "stranded"; "unstranded" is rejected with an error.
```

```
# commit 14 (this commit)
#' @param strandedness A string outlining the type of the sequencing library: stranded, or reversely_stranded. Defaults to "stranded"; "unstranded" is rejected with an error.
```

#### Line 468, `@return` of `feature_file_editor()`. Silent on what is returned.

```
# commit 13 (previous)
#' @return Outputs a new GFF3 file populated with predicted sRNAs and UTRs.
```

```
# commit 14 (this commit)
#' @return The path to the output GFF3 file, returned invisibly. The written file is a new GFF3 populated with predicted sRNAs and UTRs.
```

#### Lines 541 to 550, the commit markers. Taken last within the file.

---

### tpm_norm_flagging.R

No annotation change. The only edit is the marker removal below.

#### Lines 150 to 159, the commit markers

```
# commit 13 (previous)
#commit1 completed
#commit2 completed
#commit3 completed
#commit4 completed
#commit5 completed
#commit6 completed
#commit7 completed
#commit8 completed
#commit9 completed
#commit10 completed
```

```
# commit 14 (this commit)
# the ten marker lines and the blank separator above them are deleted; no replacement
```
