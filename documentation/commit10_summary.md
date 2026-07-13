# Commit 10: Feature-type handling and GFF robustness

## Commit Abstract: Make baerhunter's feature-type and GFF handling robust across annotation sources, correcting the neighbour-name parse that silently captured the wrong token on RefSeq, honouring rRNA and tRNA exclusion and masking on non-Ensembl annotations, and validating strandedness at the entry to feature prediction.

### Commit Summary:
Feature-type handling and GFF robustness

### Commit Description:
fix(annotation): bound the ID= capture so neighbour names parse
on any source; recognise tRNA/rRNA by biotype rather than by one
source's type names; keep tRNA/rRNA masked whatever the source;
validate strandedness at feature-prediction entry

---

## Bound the `ID=` capture so neighbour names parse on any annotation source

**Issue Summary:** The name parse in `strand_feature_editor()` was an unbounded capture whose lazy `.*?` ran to the first colon anywhere in column 9, so the token captured depended entirely on the annotation's attribute layout and the failure was silent: RefSeq `ID=gene-Rv0001;Dbxref=GeneID:885041;` wrong-matched to the Dbxref number `885041`, a colon-less RefSeq ID carried its whole attribute string, and a later-attribute colon such as Prokka's `inference=...:Prodigal` captured an arbitrary token. All of these corrupted the `upstream_feature` and `downstream_feature` fields of the augmented GFF, and where the whole attribute string is carried its own semicolons inject spurious `Name=`, `gbkey=` and `locus_tag=` pairs into column 9, malforming 3,182 predicted rows on *M. bovis* and 1,876 on SL1344; Cortes shows 0 such rows, because a numeric capture is a clean token, so it is silently wrong inside a structurally valid GFF. Quantification was unaffected, because `make_saf()` reads a separately bounded `ID=`.

**Solution Summary:** Bound the capture to the `ID=` value with `ID=([^;]*)`, then take the part after the first colon for an Ensembl `type:name` token or the whole token otherwise, so the capture can never run on into a later attribute. A purely-numeric neighbour-name guard is folded in alongside the retained no-match accumulator.

**Note:** A correctness fix, not a speed optimisation. Byte-identical on Ensembl gene rows and on every baerhunter predicted row, with the type parse byte-identical throughout since it runs only on colon-form predicted rows; on RefSeq the neighbour fields change to the gene ID, and the retained warning now fires only on a row carrying no `ID=` at all.

### feature_file_editor.R

#### The parse helper block in `strand_feature_editor()`

One bounded token function replaces the general-purpose helper, with two thin wrappers deriving the name and the type from it, and a second accumulator for purely-numeric names.

```
# commit 9 (previous)
# lines 305-314
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

```
# commit 10 (this commit)
# lines 315-341
  ## Collect GFF attributes that carry no parseable ID, and any neighbour name that
  ## comes back purely numeric (a wrong-match signature), to report each once after
  ## the loop.
  unparsed_attrs <- character(0)
  numeric_ids    <- character(0)

  ## Isolate the ID= value, bounded to the first ';', so the capture can never run
  ## on into a later attribute such as Dbxref=GeneID:885041. Records the attribute
  ## when no ID= token is present.
  parse_id_token <- function(attr) {
    token <- sub("ID=([^;]*).*", "\\1", attr)
    if (token == attr) {
      unparsed_attrs <<- c(unparsed_attrs, attr)
    }
    token
  }
  ## Feature name: the part after the first ':' for an Ensembl type:name token, or
  ## the whole token for a hyphen or bare ID (RefSeq gene-Rv0001, Prokka PROKKA_0001).
  parse_name <- function(attr) {
    name <- sub("^[^:]*:", "", parse_id_token(attr))
    if (grepl("^[0-9]+$", name)) {
      numeric_ids <<- c(numeric_ids, name)
    }
    name
  }
  ## Feature type: the part before the first ':' within the ID token.
  parse_type <- function(attr) sub(":.*", "", parse_id_token(attr))
```

#### The five parse sites

Each call loses its pattern argument and splits into the name or type wrapper. Line 346 of the previous commit is the only type parse.

```
# commit 9 (previous)
# line 340
  previous_feature_name <- parse_id(cmp_strand[nrow(cmp_strand),9], "ID=.*?:(.*?);.*")

# line 344
    feature_name <- parse_id(cmp_strand[i,9], "ID=.*?:(.*?);.*")

# line 346  (type pattern: captures before the colon)
      feature_type <- parse_id(cmp_strand[i,9], "ID=(.*?):.*?;.*")

# line 351
        next_feature_name <- parse_id(cmp_strand[i+1,9], "ID=.*?:(.*?);.*")

# line 353
        next_feature_name <- parse_id(cmp_strand[1,9], "ID=.*?:(.*?);.*")
```

```
# commit 10 (this commit)
# line 367
  previous_feature_name <- parse_name(cmp_strand[nrow(cmp_strand),9])

# line 371
    feature_name <- parse_name(cmp_strand[i,9])

# line 373
      feature_type <- parse_type(cmp_strand[i,9])

# line 378
        next_feature_name <- parse_name(cmp_strand[i+1,9])

# line 380
        next_feature_name <- parse_name(cmp_strand[1,9])
```

#### The summary warning, and the numeric guard added after it

The warning is retained but its message no longer claims IDs must take the form `ID=type:name;`, since a colon-less ID now parses. The numeric guard is new and belt-and-braces.

```
# commit 9 (previous)
# lines 372-379
  ## Report unparsable feature IDs once, as a single deduplicated summary.
  if (length(unparsed_attrs) > 0L) {
    n_failed <- length(unique(unparsed_attrs))
    warning(paste0(n_failed, " of ", nrow(cmp_strand),
                   " feature IDs could not be parsed from the GFF attribute column (e.g. ",
                   unparsed_attrs[1], "); IDs should have the form ID=type:name;."),
            call. = FALSE, immediate. = TRUE)
  }
```

```
# commit 10 (this commit)
# lines 399-413
  ## Report rows with no parseable ID once, as a single deduplicated summary.
  if (length(unparsed_attrs) > 0L) {
    n_failed <- length(unique(unparsed_attrs))
    warning(paste0(n_failed, " of ", nrow(cmp_strand),
                   " feature rows had no parseable ID= attribute (e.g. ",
                   unparsed_attrs[1], "); each feature row should carry an ID= tag."),
            call. = FALSE, immediate. = TRUE)
  }
  ## Report purely-numeric neighbour names once: a parse wrong-match signature.
  if (length(numeric_ids) > 0L) {
    warning(paste0(length(unique(numeric_ids)),
                   " feature name(s) parsed to a purely-numeric value (e.g. ",
                   numeric_ids[1], "); this is a signature of an ID-parsing wrong-match."),
            call. = FALSE, immediate. = TRUE)
  }
```

---

## Guarantee tRNA and rRNA are masked across annotation sources

**Issue Summary:** The contract of `major_features()` is that pre-annotated tRNA and rRNA are masked, but it is honoured only by accident of annotation structure: the `[^tr]RNA` regex, built only on the `original_sRNA_annotation = "unknown"` path, does not itself do the masking, since tRNA and rRNA fail to match it and are therefore merely not removed from `major_f`. They are masked instead by their retained `gene` parents on hierarchical annotations such as RefSeq and Ensembl Bacteria, and by direct retention on flat ones such as Prokka and Bakta, so a tRNA or rRNA leaks into the sRNA candidate output only if it is a child feature whose covering parent is absent or removed.

**Solution Summary:** Retain any column 3 `tRNA` or `rRNA` feature in the masking set, but filter those same types out of the neighbour-annotation set passed to `strand_feature_editor()`, so masking is guaranteed without changing the upstream and downstream fields. The column 3 test matches the child transcript rows, not the gene-level `tRNA_gene` or `rRNA_gene` wrappers, so a gene-level-only annotation is not masked by this fix; it is still excluded from quantification by the biotype recognition below, and full masking coverage awaits the configurable allow-list left as future work.

**Note:** A correctness fix, and monotone for the sRNA candidate set: it can only remove tRNA and rRNA candidates, never suppress a genuine novel sRNA. Confirmed byte-identical on *M. bovis* and SL1344, where the mask gains all 48 and 108 tRNA and rRNA rows respectively and no neighbour field changes; the only behavioural effect is on a leak annotation, a child tRNA or rRNA whose covering parent is absent, which occurs on no arm tested; there it removes a spurious candidate and may shift an adjacent UTR, with the neighbour fields still unchanged.

### feature_file_editor.R

#### The compound filter in `major_features()`

The existing filter chain is unchanged but parenthesised, so the new retention applies as a disjunction.

```
# commit 9 (previous)
# line 204
  major_f <- gff[!grepl("Parent", gff[,9], ignore.case = TRUE) & gff[,3]!='chromosome' & gff[,3]!='biological_region' & !grepl(ori_sRNA_biotype, gff[,9], ignore.case = TRUE) & gff[,3]!='region' & gff[,3]!='sequence_feature',]
```

```
# commit 10 (this commit)
# lines 213-214
  is_trna_rrna <- gff[,3] %in% c("tRNA", "rRNA")
  major_f <- gff[(!grepl("Parent", gff[,9], ignore.case = TRUE) & gff[,3]!='chromosome' & gff[,3]!='biological_region' & !grepl(ori_sRNA_biotype, gff[,9], ignore.case = TRUE) & gff[,3]!='region' & gff[,3]!='sequence_feature') | is_trna_rrna,]
```

#### The two `strand_feature_editor()` calls

Each call now receives a neighbour set with tRNA and rRNA filtered out. This is what keeps the neighbour fields byte-identical while the mask grows.

```
# commit 9 (previous)
# line 451
    plus_annot_dataframe <- strand_feature_editor("+", plus_sRNA, plus_UTR, maj_plus_features)

# line 459
    minus_annot_dataframe <- strand_feature_editor("-", minus_sRNA, minus_UTR, maj_minus_features)
```

```
# commit 10 (this commit)
# lines 495-497
    ## tRNA/rRNA are masked (above) but dropped from the neighbour set, so neighbour fields are unchanged.
    plus_neighbour_features <- maj_plus_features[!(maj_plus_features[,3] %in% c("tRNA", "rRNA")), ]
    plus_annot_dataframe <- strand_feature_editor("+", plus_sRNA, plus_UTR, plus_neighbour_features)

# lines 505-506
    minus_neighbour_features <- maj_minus_features[!(maj_minus_features[,3] %in% c("tRNA", "rRNA")), ]
    minus_annot_dataframe <- strand_feature_editor("-", minus_sRNA, minus_UTR, minus_neighbour_features)
```

#### Documentation comments

The retention changes what `major_features()` returns, so the roxygen contract and the inline filter comment are corrected to match.

```
# commit 9 (previous)
# line 179 (description), line 186 (@return), line 203 (inline comment)
#' The function extracts parent features only; it also excludes all non-coding RNAs that are already annotated in the file.
#' @return A dataframe with the major features from a set strand.
  ## Select only the major genomic features: remove all child features (like CDS, mRNA etc.), previously annotated sRNAs and extra features
```

```
# commit 10 (this commit)
# line 188 (description), line 195 (@return), line 212 (inline comment)
#' The function extracts parent features, plus any tRNA/rRNA rows (retained for masking even when child features); it also excludes non-coding RNAs already annotated in the file.
#' @return A dataframe with the major features for a set strand, plus any tRNA/rRNA rows retained for masking.
  ## Select the major genomic features: remove child features (CDS, mRNA etc.), previously annotated sRNAs and extra features, but retain tRNA/rRNA for masking.
```

---

## Recognise structural RNA by biotype rather than by one source's type names

**Issue Summary:** `make_saf()` honoured `excl_rna = TRUE` by matching column 3 against the literal strings `tRNA_gene` and `rRNA_gene`, which belong to legacy Ensembl and the Sequence Ontology alone, so on RefSeq (`gene` plus `gene_biotype=`) and modern Ensembl Bacteria (`ncRNA_gene` plus `biotype=`) the match found nothing and the function silently excluded nothing. Because `excl_rna` defaults to `TRUE`, this was the default path.

**Solution Summary:** Key the exclusion on a prefix-agnostic `biotype=(tRNA|rRNA)(;|$)` tag, which catches both the RefSeq `gene_biotype=` and the Ensembl `biotype=` forms in a single test, and keep a type-name list alongside it for flat Prokka and Bakta rows and for the legacy Sequence Ontology terms. A zero-match warning is added for the case where `excl_rna = TRUE` excludes nothing at all.

**Note:** A correctness fix. Byte-identical wherever the existing code already worked, including baerhunter's own bundled EnsemblBacteria r40 vignette annotation, where 48 features are excluded under both the old and the new form; it is a behavioural change on the default path for RefSeq, modern Ensembl and flat Prokka, and on Cortes H37Rv, where rRNA and tRNA carry about 60 per cent of the TPM budget, every retained gene's TPM rescales by roughly 2.5 times. `tmRNA` is out of scope and the `(tRNA|rRNA)` pattern correctly does not match it.

### count_features.R

#### Structural-RNA vocabulary constants

Two internal seed constants are added at the top of the file, consumed by `make_saf()` alone. `feature_file_editor.R` carries its own literals, because its three sites perform different operations on a narrower set; unifying them is future work.

```
# commit 9 (previous)
# no equivalent; the type names were inline literals inside make_saf()
```

```
# commit 10 (this commit)
# lines 1-7
## Structural-RNA vocabulary for make_saf(): the feature types treated as tRNA
## or rRNA when excl_rna = TRUE. Transcript-level types, gene-level types and
## the biotype tag are all recognised, so the exclusion holds across annotation
## sources (Ensembl tRNA_gene / rRNA_gene, flat Prokka/Bakta tRNA / rRNA,
## and RefSeq/Ensembl biotype=tRNA / biotype=rRNA attributes).
.bh_structural_rna_types      <- c("tRNA", "rRNA")
.bh_structural_rna_gene_types <- c("tRNA_gene", "rRNA_gene")
```

#### The exclude branch of `make_saf()`

Only the `else` branch changes. The shared base filter is lifted into `base_keep` so the exclusion mask can be tested against it, which is what the zero-match warning keys on.

```
# commit 9 (previous)
# lines 122-133
  if (!exclude){
    major_f <- gff[!grepl("Parent", gff[,9], ignore.case = TRUE) &
                     gff[,3]!='chromosome' & gff[,3]!='biological_region' &
                     gff[,3]!='region' & gff[,3]!='sequence_feature',]
  }else{ # if you want to exclude rRNA and tRNA features)
    major_f <- gff[!grepl("Parent", gff[,9], ignore.case = TRUE) &
                     gff[,3]!='chromosome' & gff[,3]!='biological_region' &
                     gff[,3]!='region' &
                     gff[,3]!='sequence_feature' &
                     gff[,3]!='tRNA_gene' &
                     gff[,3]!='rRNA_gene',]
  }
```

```
# commit 10 (this commit)
# lines 130-146
  if (!exclude){
    major_f <- gff[!grepl("Parent", gff[,9], ignore.case = TRUE) &
                     gff[,3]!='chromosome' & gff[,3]!='biological_region' &
                     gff[,3]!='region' & gff[,3]!='sequence_feature',]
  }else{ # exclude structural RNA (tRNA/rRNA) across annotation sources
    base_keep <- !grepl("Parent", gff[,9], ignore.case = TRUE) &
                   gff[,3]!='chromosome' & gff[,3]!='biological_region' &
                   gff[,3]!='region' & gff[,3]!='sequence_feature'
    excluded  <- gff[,3] %in% c(.bh_structural_rna_types, .bh_structural_rna_gene_types) |
                   grepl("biotype=(tRNA|rRNA)(;|$)", gff[,9])
    if (!any(base_keep & excluded)) {
      warning("excl_rna = TRUE but no tRNA or rRNA features were found to exclude; ",
              "the annotation may contain none, or may use RNA type names or biotype ",
              "tags this function does not recognise.", call. = FALSE, immediate. = TRUE)
    }
    major_f <- gff[base_keep & !excluded,]
  }
```

---

## Validate `strandedness` at feature-prediction entry

**Issue Summary:** After the Commit 3 rewrite, `peak_union_calc()` silently treated any unrecognised `strandedness` value as ordinary forward-stranded data, and its own default, `"unstranded"`, was itself such a value, so the default was an unhandled one.

**Solution Summary:** Validate `strandedness` at function entry, raise an explicit error naming the unrecognised value and the valid set, and align the `peak_union_calc()` default from `"unstranded"` to `"stranded"`. The guard sits at both `feature_file_editor()`, the primary entry point, and `peak_union_calc()`, which is exported and so can be called directly, with an identical message and accepted set in each.

**Note:** A robustness fix. Byte-identical on valid input, with only the invalid path changing: unstranded was never a working prediction mode, since in the baseline it matched no branch and aborted on the first `coverage()` call with an opaque error. `find_strandedness()` in the counting stage legitimately keeps `"unstranded"`, a principled stage difference the error message states explicitly.

### feature_file_editor.R

#### `peak_union_calc()`: the default, and the entry guard

The default moves to a valid value, and the guard sits between `match.arg()` and the BAM-file listing, so all entry validation fires before any file access.

```
# commit 9 (previous)
# line 52; the default is itself an unhandled value
peak_union_calc <- function(bam_location = ".", bam_txt_list = "", low_coverage_cutoff, high_coverage_cutoff, peak_width, paired_end_data = FALSE, strandedness = "unstranded", scanbamparam = NULL, mapqFilter = 10, coverage_model = c("fragment", "footprint")) {
  coverage_model <- match.arg(coverage_model)
  ## Find all BAM files in the directory.
```

```
# commit 10 (this commit)
# lines 52-63
peak_union_calc <- function(bam_location = ".", bam_txt_list = "", low_coverage_cutoff, high_coverage_cutoff, peak_width, paired_end_data = FALSE, strandedness = "stranded", scanbamparam = NULL, mapqFilter = 10, coverage_model = c("fragment", "footprint")) {
  coverage_model <- match.arg(coverage_model)
  ## Validate strandedness at entry, before any BAM is read.
  valid_strandedness <- c("stranded", "reversely_stranded")
  if (!strandedness %in% valid_strandedness) {
    stop("Invalid 'strandedness' value: '", strandedness,
         "'. Feature prediction requires stranded data; must be one of: ",
         paste(valid_strandedness, collapse = ", "),
         ". ('unstranded' is valid for read counting but not for prediction.)",
         call. = FALSE)
  }
  ## Find all BAM files in the directory.
```

#### `feature_file_editor()`: the second guard

Purely additive: the signature already defaulted to `"stranded"`, so only the guard block is new. It is the first statement in the function, ahead of the BAM-directory listing.

```
# commit 9 (previous)
# line 438 signature (default already valid); no guard, the function body begins
# directly with the BAM-directory listing
feature_file_editor <- function(bam_directory = ".", bam_list = "", original_annotation_file, annot_file_dir = ".", output_file, original_sRNA_annotation, low_coverage_cutoff, high_coverage_cutoff, min_sRNA_length, min_UTR_length, paired_end_data = FALSE, strandedness  = "stranded", scanbamparam = NULL, mapqFilter = 10, coverage_model = c("fragment", "footprint")) {
  test <- list.files(path = bam_directory, pattern = "\\.BAM$", full.names = TRUE, ignore.case = TRUE)
```

```
# commit 10 (this commit)
# lines 472-483 (signature unchanged; the guard block is inserted)
feature_file_editor <- function(bam_directory = ".", bam_list = "", original_annotation_file, annot_file_dir = ".", output_file, original_sRNA_annotation, low_coverage_cutoff, high_coverage_cutoff, min_sRNA_length, min_UTR_length, paired_end_data = FALSE, strandedness  = "stranded", scanbamparam = NULL, mapqFilter = 10, coverage_model = c("fragment", "footprint")) {
  ## Validate strandedness at entry, before any filesystem work. peak_union_calc()
  ## validates too, so a direct caller of that function is guarded as well.
  valid_strandedness <- c("stranded", "reversely_stranded")
  if (!strandedness %in% valid_strandedness) {
    stop("Invalid 'strandedness' value: '", strandedness,
         "'. Feature prediction requires stranded data; must be one of: ",
         paste(valid_strandedness, collapse = ", "),
         ". ('unstranded' is valid for read counting but not for prediction.)",
         call. = FALSE)
  }
  test <- list.files(path = bam_directory, pattern = "\\.BAM$", full.names = TRUE, ignore.case = TRUE)
```

#### Documentation comments

The roxygen `@param strandedness` line is identical in both functions and is updated in both, so the documented contract states the new default and records that `"unstranded"` is now rejected.

```
# commit 9 (previous)
# line 12 (peak_union_calc) and line 403 (feature_file_editor), identical text
#' @param strandedness A string outlining the type of the sequencing library: stranded, or reversely stranded.
```

```
# commit 10 (this commit)
# line 12 (peak_union_calc) and line 437 (feature_file_editor), identical text
#' @param strandedness A string outlining the type of the sequencing library: stranded, or reversely stranded. Defaults to "stranded"; "unstranded" is rejected with an error.
```
