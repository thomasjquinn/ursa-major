# Notes

This document highlights important project details for the updated R files.

## Summary Files Notes

The per commit summary files (`commit1_summary.md` through `commit9_summary.md`) are located in the `docs` folder of the repository. These files detail the code changes by line number for each commit and the reason these changes were implemented.

## Performance Notes

Across Commits 1 to 9, most changes gave only a small speed gain on their own. The bulk of the overall speedup comes from a handful of changes that are detailed below.

**`coverage()` recomputed across strand calls (Commit 3).** This change reads each BAM once and splits by strand rather than reading every file twice, and most of the pipeline speedup is attributable to it. On the single-end Cortes _M. tuberculosis_ H37Rv profiling run (E-MTAB-1616) it accounts for roughly 80 per cent of the time saving in `feature_file_editor`, the most resource-heavy stage of the pipeline. This is a textbook example of the Pareto Principle.

**Vector-growing loop in `tpm_flagging` (Commit 4).** This change, together with the row-filter change below, accounts for most of the speedup in the `tpm_norm_flagging.R` file, which is now far faster than before, although its total runtime is very short. It is faster because the old code added each annotation line to a growing list one at a time, whereas the new code builds the whole list in a single step via vectorisation.

**Vectorise the row filter in `tpm_flag_filtering` (Commit 4).** This change shares that `tpm_norm_flagging.R` speedup with the loop change detailed above. It is faster for the same reason: instead of checking the table one row at a time in a loop, the new code tests every row at once across whole columns, again using vectorisation.

## Behavioural Change Notes

Unlike the speed and robustness changes elsewhere in the package, the following updates may alter the output from the baseline version of baerhunter.

**TPM-to-GFF matching regex requires a colon in the feature ID (Commit 5).** In `tpm_flagging()` the ID-extraction regex required a colon, so on the project's RefSeq H37Rv annotation (GCF\_000195955.2, sequence NC\_000962.3), whose IDs are colon-less (`ID=gene-Rv0001;`), the gene set was silently passed through unflagged; the fix makes the capture colon-agnostic so it matches how `make_saf()` builds the count-table IDs. On the single-end Cortes _M. tuberculosis_ H37Rv test data (E-MTAB-1616, six BAMs) this newly flagged roughly 4060 gene-level features that the baseline and Commit 4 had left unflagged (4060 of 10848 lines in `flagged.gff3` and 4060 of 8448 in `filtered.gff3`), while the four non-flagging outputs stayed byte-identical.

**Plus-strand UTR length filter mismatch (Commit 7).** In current baerhunter the plus-strand `UTR_calc()` call filters at `min_sRNA_length` (40) instead of `min_UTR_length` (50), so the updated version removes the plus-strand putative UTRs whose length falls in the 40 to 49 nt gap and brings the plus strand into line with the minus-strand cut-off. This was proven out on a paired-end _M. tuberculosis_ H37Rv dataset (run ERR2103718, project PRJEB65014 / E-MTAB-6011), where the updated version dropped exactly 58 plus-strand `putative_UTR` features of 40 to 49 bp (58 plus, 0 minus) that the baseline version did not, matching the predicted behaviour.

**Set `countMultiMappingReads` and `countReadPairs` explicitly (Commit 7).** When baerhunter was first published, the Bioconductor package it uses for counting (Rsubread) had a particular default for handling multi-mapping reads, but that default later changed, so the unmodified code could produce different counts depending on which Rsubread version happened to be installed. This fix hardcodes the original settings (`countMultiMappingReads = FALSE` and `countReadPairs = TRUE`) so baerhunter once again behaves as originally intended and returns the same counts on every Rsubread version.

**Read-quality filter for paired-end peak detection (Commit 8).** This adds a read-quality filter to the BAM-reading step, controlled by the `mapqFilter` setting (default 10), so low-mapping-quality and multi-mapping reads no longer inflate the coverage that peaks are called from. Because different aligners report mapping quality on different scales, the user needs to set `mapqFilter` to a value appropriate for the aligner that produced their BAM files earlier in the pipeline, since a threshold set too high for a given aligner can discard every read.
