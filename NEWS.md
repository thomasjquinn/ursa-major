# baerhunter (ursa-major development version)

Development codename for the optimised version of baerhunter, tracking changes made
relative to the `v0-baseline` tag (upstream v0.9.1).

* Each GFF annotation is now read and parsed once per run and the result
  reused, via a new exported helper `load_gff_cache()`. The speed gain grows
  with annotation size. Predicted features, counts, TPM values, flags, and
  filtered output are unchanged for valid input.

* For paired-end data, `peak_union_calc()` and `feature_file_editor()` gain a
  `coverage_model` argument. The default, `"fragment"`, reproduces the existing
  behaviour: each read pair contributes coverage across the whole fragment,
  including the unsequenced gap between the mates. The alternative,
  `"footprint"`, counts coverage over the two mate alignments only and excludes
  the gap, the model most standard RNA-seq coverage tools use. The argument has
  no effect on single-end data, so the predicted sRNAs and UTRs are unchanged
  unless `coverage_model = "footprint"` is requested.

* Peak detection now filters reads by mapping quality and alignment flag
  before building coverage. `peak_union_calc()` and `feature_file_editor()`
  gain a `scanbamparam` argument and a `mapqFilter` argument. When `scanbamparam`
  is not supplied, a default `ScanBamParam` is built that keeps only primary,
  mapped, quality-passing reads with a mapping quality of at least `mapqFilter`
  (default 10), dropping secondary and supplementary alignments, and for
  paired-end data additionally requiring a properly paired read with a mapped
  mate. Multi-mapping reads, which aligners commonly mark with a mapping quality
  of zero, are therefore excluded by default. This changes the predicted sRNAs
  and UTRs for any run whose BAM files contain low-quality or multi-mapping
  alignments, on both single-end and paired-end data. Set `mapqFilter` no
  higher than the mapping quality your aligner assigns a unique read (see
  `?feature_file_editor` for per-aligner values; Bowtie2 users should set
  `mapqFilter = 1`). Pass `mapqFilter = NA` to disable the quality filter, or
  `scanbamparam = ScanBamParam()` to restore the previous unfiltered behaviour.
  If too high a threshold leaves a BAM file with no reads, a warning now names
  that file rather than the run continuing on empty coverage.

* `count_features()` gains two optional arguments, `count_multi_mapping_reads`
  and `count_read_pairs`, which set the matching featureCounts options
  explicitly so counts are reproducible across Rsubread versions. Multi-mapping
  reads are now excluded by default (`count_multi_mapping_reads = FALSE`),
  restoring the behaviour of the original analyses. Current Rsubread versions
  count multi-mapping reads by default, which the package's `fraction = TRUE`
  setting then weights fractionally, so excluding them lowers counts at
  repetitive and paralogous loci relative to the current Rsubread default.
  `count_read_pairs` defaults to `TRUE`, the
  existing fragment-counting behaviour, and affects paired-end data only. Using
  `count_read_pairs` requires Rsubread 2.4.3 or later, the first release to
  provide the argument; it exposes the `--countReadPairs` option introduced in
  the Subread tool at version 2.0.2.

* Predicted UTRs on the plus strand are now filtered at the minimum UTR
  length, the same cut-off already applied to the minus strand. A parameter
  mix-up previously filtered plus-strand UTRs at the minimum sRNA length
  instead, so when the two lengths differ, short plus-strand UTRs between the
  two cut-offs were retained on one strand only. This correction removes that
  asymmetry, so the predicted UTR set will differ from earlier versions for any
  run where the minimum sRNA and UTR lengths are not equal. The minus-strand
  output is unchanged.

* `tpm_flagging()` now returns the path of the file it writes, invisibly,
  instead of `NULL`, so it composes more cleanly when called from scripts. It
  continues to write the same flagged annotation.

* The count table, count summary, and TPM table are now written without
  quotation marks around text labels such as feature IDs and sample names. A
  reader using `read.delim()` parses the files to the same data as before; only
  a strict byte comparison with older output differs.

* `tpm_normalisation()` no longer accepts a `feature_type` argument. It was
  documented as selecting feature types from the annotation but had no effect on
  the result, so any value passed was silently ignored. Removing it is a
  breaking change for a call that supplied `feature_type` by name or relied on
  argument position; the computed TPM table is unchanged.

* `count_features()` gains three optional arguments, `largest_overlap`,
  `frac_overlap_feature`, and `read_to_pos`, exposing the matching featureCounts
  overlap settings. Each defaults to the value featureCounts already uses, so
  counts are unchanged unless set. These replace the previous undocumented route
  of passing the camelCase names through `...`, which now stops with a
  duplicate-argument error. Setting `largest_overlap = TRUE` with the package's
  fractional counting is reliable only on Rsubread 2.14.0 or later, which fixed
  an earlier silent miscount of that combination.

* BAM files are discovered using a pattern that matches a literal `.bam`/`.BAM`
  extension, so a file whose name ends in `bam` without a separating dot, such
  as `sampleBAM`, is no longer picked up by mistake. Directories of normally
  named BAM files are unaffected.

* Feature IDs in the annotation are parsed more robustly. Expression flagging
  now matches a feature ID whether or not it carries a `type:` prefix, so
  annotations with bare or hyphen-delimited IDs, such as RefSeq-style GFFs, are
  flagged correctly rather than being silently skipped; for Ensembl-style IDs
  the flagged output is unchanged. Two previously silent failures are now
  reported: a warning is emitted if no feature ID matches the count table during
  flagging, and a single summary warning if any feature IDs cannot be parsed
  during strand annotation.

* Filtering flagged features now matches the chosen flag as a literal,
  case-sensitive string rather than a case-insensitive regular expression. The
  standard flag names match exactly as before; the difference is visible only
  for a flag supplied in a different case, or one containing regular-expression
  metacharacters.

* TPM normalisation, expression flagging, and flag-based filtering are faster,
  particularly with many samples or large annotations. The TPM values, flags,
  and filtered features are unchanged for valid input.

* `peak_union_calc()` now computes both strands in a single call: it returns
  a named list with `plus` and `minus` peak coordinates and no longer takes a
  `target_strand` argument. This is a breaking change for any code that calls
  `peak_union_calc()` directly; the package's own annotation pipeline is
  updated to match.

* Peak detection is substantially faster, with predicted sRNAs and UTRs
  unchanged for valid input.

* Functions whose job is to write an output file now return invisibly rather
  than the literal string `"Done!"`, so they no longer auto-print at the
  console and compose more cleanly when called from scripts.

* Status messages shown while predicting and annotating sRNAs and UTRs now
  use `message()` instead of `print()`, so they are written to standard
  error and can be silenced with `suppressMessages()`. The message text is
  unchanged.

* Functions that detect invalid input now halt with `stop()` rather than
  returning a message string that was passed downstream.

* Removed the `assertthat` dependency; input checks now use base R's
  `stopifnot()`. User-facing error messages are unchanged.
