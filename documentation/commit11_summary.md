# Commit 11: Parameter Scout

## Commit Abstract: Add Parameter Scout, an optional advisory tool that reports the coverage percentiles of the intergenic regions per BAM file so that the two coverage cutoffs can be read off the data rather than guessed, shipped as a paired-end and a single-end sibling and modifying no existing file.

### Commit Summary:

Add Parameter Scout coverage cutoff advisor

### Commit Description:

feat(parameter-scout): report intergenic coverage percentiles  
per BAM and offer two cutoff pairs, as paired-end and single-end  
siblings; advisory only, writes no annotation and sets no  
parameter; no existing file is modified

---

## Parameter Scout Summary

baerhunter requires the user to supply `low_coverage_cutoff` and `high_coverage_cutoff` and offers no way to derive them, and because they are raw read depths, unscaled to library or genome, a pair that suits one dataset can be a poor choice for another. Parameter Scout restricts per-base coverage to the intergenic regions, where baerhunter predicts features, and reports its percentiles for each BAM file, so the cutoffs can be chosen from the data. Parameter Scout is designed to be run as an optional step 0, before the baerhunter pipeline. It is advisory only: it writes no annotation, sets no parameter and passes nothing to the rest of the pipeline, so the user records the values and types them in by hand.

Full details and complete end-user instructions are in `documentation/parameter_scout_instructions.md`.

---

## Parameter Scout Paired End

The two files are siblings that share every line of their logic, and the difference that matters is which reads each one is built to read. `parameter_scout_paired_end.R` must be the file run on paired-end data: it requires properly paired alignments and measures coverage across the whole fragment, so pointing it at single-end data discards every read and returns `NA` percentiles, while pointing the single-end sibling at paired-end data counts each mate as an independent record and inflates every number reported.

The two files differ in fifteen lines and in nothing else. Twelve of the fifteen are comment or roxygen lines: the four header lines, the marker on the entry point's own section, two roxygen lines on the entry point itself, and five roxygen cross-references naming the scan in the documentation of `.write_scout_outputs()`, `suggest_cutoffs()` and `plot_scout_distribution()`. All but two of those twelve name something different rather than say something different, so they are not shown here. Only three lines of executable code differ: the name of the exported function, the default of `paired_end_data`, and the entry-point name inside one error message. The comment-stripped executable code is 854 lines in each file, and those three are the only lines in it that differ.

### parameter\_scout\_paired\_end.R

#### The three executable lines, 828, 831 and 869

The default of `paired_end_data` is the whole of the behavioural difference: every branch on library type sits inside `.strand_split_reads()`, which takes the value as an argument. The error message names the function the caller actually called.

```
# parameter_scout_single_end.R (sibling)
# line 828
parameter_scout_single_end <- function(bam_location = ".", bam_txt_list = "",
# line 831
                            paired_end_data = FALSE,
# line 869
      stop("parameter_scout_single_end: could not create out_dir: ", out_dir, call. = FALSE)
```

```
# parameter_scout_paired_end.R (this file)
# line 828
parameter_scout_paired_end <- function(bam_location = ".", bam_txt_list = "",
# line 831
                            paired_end_data = TRUE,
# line 869
      stop("parameter_scout_paired_end: could not create out_dir: ", out_dir, call. = FALSE)
```

---

## Parameter Scout Single End

The two files are siblings that share every line of their logic, and the difference that matters is which reads each one is built to read. `parameter_scout_single_end.R` must be the file run on single-end data: reading paired-end data through it counts each mate as an independent record, with no proper-pair requirement and no fragment span, which inflates every percentile it reports.

The two files differ in fifteen lines and in nothing else. Twelve of the fifteen are comment or roxygen lines: the four header lines, the marker on the entry point's own section, two roxygen lines on the entry point itself, and five roxygen cross-references naming the scan in the documentation of `.write_scout_outputs()`, `suggest_cutoffs()` and `plot_scout_distribution()`. All but two of those twelve name something different rather than say something different, so they are not shown here. Only three lines of executable code differ: the name of the exported function, the default of `paired_end_data`, and the entry-point name inside one error message. The comment-stripped executable code is 854 lines in each file, and those three are the only lines in it that differ.

### parameter\_scout\_single\_end.R

#### The three executable lines, 828, 831 and 869

The exported name differs because a package collates every file in `R/` into one namespace, and two files defining one name leave only the later-collated one standing, silently. The default of `paired_end_data` is the whole of the behavioural difference: every branch on library type sits inside `.strand_split_reads()`, which takes the value as an argument, and in this file that branch also inverts the strand assignment when the library is reversely stranded. The error message names the function the caller actually called.

```
# parameter_scout_paired_end.R (sibling)
# line 828
parameter_scout_paired_end <- function(bam_location = ".", bam_txt_list = "",
# line 831
                            paired_end_data = TRUE,
# line 869
      stop("parameter_scout_paired_end: could not create out_dir: ", out_dir, call. = FALSE)
```

```
# parameter_scout_single_end.R (this file)
# line 828
parameter_scout_single_end <- function(bam_location = ".", bam_txt_list = "",
# line 831
                            paired_end_data = FALSE,
# line 869
      stop("parameter_scout_single_end: could not create out_dir: ", out_dir, call. = FALSE)
```

