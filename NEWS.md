# baerhunter (ursa-major development version)

Development codename for the optimised version of baerhunter, tracking changes made
relative to the `v0-baseline` tag (upstream v0.9.1).

* New optional module, Parameter Scout, which reports the coverage percentiles
  of the intergenic regions for each BAM file so that `low_coverage_cutoff` and
  `high_coverage_cutoff` can be chosen from the data rather than guessed. It is
  run before the rest of the pipeline and ships as two files,
  `parameter_scout_paired_end.R` and `parameter_scout_single_end.R`, exporting
  `parameter_scout_paired_end()` or `parameter_scout_single_end()` alongside
  `suggest_cutoffs()` and `plot_scout_distribution()`. Run the one that matches
  the library: reading paired-end data through the single-end module counts each
  mate as an independent record and inflates every number reported, and reading
  single-end data through the paired-end module discards every read. A scan
  writes a percentile table, a run-length summary of the coverage distribution,
  two suggested cutoff pairs and a figure, and returns them invisibly. Coverage
  is built under the same read filter and coverage model `peak_union_calc()`
  uses, so the percentiles are on the scale the pipeline later applies; the
  arguments controlling that must be given the same values in both places. The
  module is advisory: it writes no annotation, sets no parameter and passes
  nothing to the rest of the pipeline, so the chosen values are typed into
  `feature_file_editor()` by hand. Nothing existing behaves differently, and
  deleting the two files leaves the package unchanged. Until the packaging pass
  they are loaded with `source()` rather than attached with the package. Full
  instructions are in `parameter_scout_instructions.md`.

* Feature prediction now rejects an `"unstranded"` or otherwise unrecognised
  `strandedness` value with an explicit error naming the offending value and the
  accepted set, in place of the opaque failure that occurred when no strand
  branch matched. The `peak_union_calc()` default changes from `"unstranded"` to
  `"stranded"`, since unstranded was never a working prediction mode; read
  counting still accepts `"unstranded"`, which remains a valid counting mode.

* On RefSeq-style annotations the `upstream_feature` and `downstream_feature`
  fields in the augmented GFF previously carried the wrong value (a `Dbxref`
  number where the neighbouring gene row had one, otherwise that gene's whole
  attribute block) and now carry the gene ID. Ensembl-style `type:name` IDs and
  baerhunter's own predicted features are unchanged, and counts and TPM were
  never affected, since quantification uses a separately bounded ID parse. Two
  previously silent conditions are now reported: a warning if a feature name
  parses to a purely numeric value, the signature of the old wrong-match, and
  the existing unparsable-ID warning is reworded, now firing only on a feature
  row carrying no `ID=` attribute at all.

* `excl_rna = TRUE`, the default, now excludes tRNA and rRNA on all annotation
  sources by recognising them from their biotype rather than only the legacy
  Ensembl `tRNA_gene`/`rRNA_gene` type names. On RefSeq, modern Ensembl, and
  flat Prokka or Bakta annotations, where the exclusion previously matched
  nothing and silently kept those genes, the structural-RNA rows are now removed
  from the count matrix and the TPM denominator, so every retained gene's TPM
  rises (roughly 2.5-fold on the Cortes H37Rv arm, where tRNA and rRNA carry
  about 60 per cent of the TPM budget). Legacy Ensembl annotations,
  including the bundled vignette data, are unchanged, and a warning is now
  emitted if `excl_rna = TRUE` finds nothing to exclude.

* Pre-annotated tRNA and rRNA are now retained in the set of features masked
  from the sRNA candidate search, rather than being masked only as a side effect
  of the annotation's structure. Previously a tRNA or rRNA could be predicted as
  a novel sRNA where its covering parent feature was absent or removed, as on
  flat annotations with no gene-level parent. Those rows are filtered out of the
  neighbour set passed to the annotation step, so the `upstream_feature` and
  `downstream_feature` fields are unchanged. The effect on the predicted set is
  monotone: it can remove a tRNA or rRNA candidate and may shift an adjacent
  UTR, but it can never suppress a genuine novel sRNA.

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
