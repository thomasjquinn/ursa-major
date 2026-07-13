# Notes

This document highlights important project details for the updated R files.

## Summary Files Notes

The per commit summary files (`commit1_summary.md` through `commit10_summary.md`) are located in the `documentation` folder of the repository. These files detail the code changes by line number for each commit and the reason these changes were implemented.

## Performance Notes

Across Commits 1 to 10, most changes gave only a small speed gain on their own. The bulk of the overall speedup comes from a handful of changes that are detailed below.

**`coverage()` recomputed across strand calls (Commit 3).** This change reads each BAM once and splits by strand rather than reading every file twice, and most of the pipeline speedup is attributable to it. On the single-end Cortes _M. tuberculosis_ H37Rv profiling run (E-MTAB-1616) it accounts for roughly 80 per cent of the time saving in `feature_file_editor`, the most resource-heavy stage of the pipeline. This is a textbook example of the Pareto Principle.

**Vector-growing loop in `tpm_flagging` (Commit 4).** This change, together with the row-filter change below, accounts for most of the speedup in the `tpm_norm_flagging.R` file, which is now far faster than before, although its total runtime is very short. It is faster because the old code added each annotation line to a growing list one at a time, whereas the new code builds the whole list in a single step via vectorisation.

**Vectorise the row filter in `tpm_flag_filtering` (Commit 4).** This change shares that `tpm_norm_flagging.R` speedup with the loop change detailed above. It is faster for the same reason: instead of checking the table one row at a time in a loop, the new code tests every row at once across whole columns, again using vectorisation.

## Behavioural Change Notes

Unlike the speed and robustness changes elsewhere in the package, the following updates may alter the output from the baseline version of baerhunter.

**TPM-to-GFF matching regex requires a colon in the feature ID (Commit 5).** In `tpm_flagging()` the ID-extraction regex required a colon, so on the project's RefSeq H37Rv annotation (GCF\_000195955.2, sequence NC\_000962.3), whose IDs are colon-less (`ID=gene-Rv0001;`), the gene set was silently passed through unflagged; the fix makes the capture colon-agnostic so it matches how `make_saf()` builds the count-table IDs. On the single-end Cortes _M. tuberculosis_ H37Rv test data (E-MTAB-1616, six BAMs) this newly flagged roughly 4060 gene-level features that the baseline and Commit 4 had left unflagged (4060 of 10848 lines in `flagged.gff3` and 4060 of 8448 in `filtered.gff3`), while the four non-flagging outputs stayed byte-identical.

**Plus-strand UTR length filter mismatch (Commit 7).** In current baerhunter the plus-strand `UTR_calc()` call filters at `min_sRNA_length` (40) instead of `min_UTR_length` (50), so the updated version removes the plus-strand putative UTRs whose length falls in the 40 to 49 nt gap and brings the plus strand into line with the minus-strand cut-off. This was proven out on a paired-end _M. tuberculosis_ H37Rv dataset (run ERR2103718, project PRJEB65014 / E-MTAB-6011), where the updated version dropped exactly 58 plus-strand `putative_UTR` features of 40 to 49 bp (58 plus, 0 minus) that the baseline version did not, matching the predicted behaviour.

**Set `countMultiMappingReads` and `countReadPairs` explicitly (Commit 7).** When baerhunter was first published, the Bioconductor package it uses for counting (Rsubread) had a particular default for handling multi-mapping reads, but that default later changed, so the unmodified code could produce different counts depending on which Rsubread version happened to be installed. This fix hardcodes the original settings (`countMultiMappingReads = FALSE` and `countReadPairs = TRUE`) so baerhunter once again behaves as originally intended and returns the same counts on every Rsubread version.

**Read-quality filter for paired-end peak detection (Commit 8).** This adds a read-quality filter to the BAM-reading step, controlled by the `mapqFilter` setting (default 10), so low-mapping-quality and multi-mapping reads no longer inflate the coverage that peaks are called from. Because different aligners report mapping quality on different scales, the user needs to set `mapqFilter` to a value appropriate for the aligner that produced their BAM files earlier in the pipeline, since a threshold set too high for a given aligner can discard every read.

**Neighbour-name parsing on non-Ensembl annotations (Commit 10).** The parse that fills the `upstream_feature` and `downstream_feature` fields of the augmented GFF captured an unbounded stretch of the attribute column, so the token it picked up depended on how the annotation happened to lay out its attributes. On RefSeq it wrong-matched to the `Dbxref` number rather than the gene ID, and where no colon was present it carried the whole attribute block, injecting spurious tags into column 9. The capture is now bounded to the `ID=` value, so those two fields carry the gene ID on any annotation source. Ensembl-style `type:name` IDs and baerhunter's own predicted features are unchanged, and counts and TPM were never affected, because quantification parses the ID separately.

**Structural RNA excluded by biotype rather than by type name (Commit 10).** `excl_rna = TRUE`, the default, recognised tRNA and rRNA only by the legacy Ensembl type names `tRNA_gene` and `rRNA_gene`, so on RefSeq, modern Ensembl and flat Prokka or Bakta annotations it matched nothing and silently kept those genes in the count matrix. The exclusion is now keyed on the biotype tag as well as the type names, so structural RNA is removed on all of these sources. Because tRNA and rRNA sit in the TPM denominator, removing them raises the TPM of every retained gene, by roughly 2.5 times on the Cortes _M. tuberculosis_ H37Rv arm, where they carry about 60 per cent of the TPM budget. Legacy Ensembl annotations, including baerhunter's own vignette data, are unchanged, and a warning is now raised if the exclusion finds nothing to exclude.

**tRNA and rRNA masking guaranteed across annotation sources (Commit 10).** Pre-annotated tRNA and rRNA were masked from the sRNA candidate search only as a side effect of annotation structure, so a tRNA or rRNA could leak into the predictions where its covering parent feature was absent. Those rows are now retained in the masking set directly, and filtered out of the neighbour set so the `upstream_feature` and `downstream_feature` fields do not change. The effect is monotone: it can only remove tRNA and rRNA candidates, never suppress a genuine novel sRNA. No arm tested produced such a leak, so the predictions on the project's datasets are unchanged.

**Strandedness now validated before prediction (Commit 10).** Feature prediction silently treated any unrecognised `strandedness` value as ordinary forward-stranded data, and the `peak_union_calc()` default, `"unstranded"`, was itself such a value. The value is now checked at the entry to prediction and an unrecognised one raises an explicit error naming the accepted set, and the `peak_union_calc()` default changes from `"unstranded"` to `"stranded"`. Read counting still accepts `"unstranded"`, which remains a valid counting mode; unstranded was never a working prediction mode.
