# ursa-major

**Improving the speed and usability of the baerhunter software for identifying unannotated expressed regions in bacterial transcriptomes**

## About

`ursa-major` is the working codename for an MSc Bioinformatics thesis project at Birkbeck, University of London, which aims to improve the performance and functionality of [baerhunter](https://github.com/irilenia/baerhunter).

baerhunter is an R package that uses a coverage-based method to predict, annotate, and filter unannotated expressed regions, such as small RNAs (sRNAs) and untranslated regions (UTRs), from bacterial RNA-seq data. This project revisits the package to make it faster, more robust, and easier to use, while preserving its existing scientific output.

The work is supervised by Dr Irilenia Nobeli at Birkbeck, University of London. The updates have been validated against _Mycobacterium tuberculosis_ H37Rv RNA-seq data. The improved pipeline has since been applied to _Mycobacterium bovis_ AF2122/97 RNA-seq data as a biological test case. Runtime has additionally been measured on paired-end _M. tuberculosis_ H37Rv RNA-seq (run ERR2103718, from E-MTAB-6011 / PRJEB65014).

## Installation

```r
install.packages(c("remotes", "BiocManager"))
options(repos = BiocManager::repositories())
remotes::install_github("thomasjquinn/ursa-major", upgrade = "never")
```

The package installs as `baerhunter`, the same name as upstream:

```r
library(baerhunter)
packageVersion("baerhunter")
#> '0.9.2.0000'
```

The package requires R 4.5 or later and Bioconductor 3.22 or later. `install_github` resolves CRAN dependencies on its own but not Bioconductor ones, which is why the install block above sets `repos`. The dependency versions are declared in `DESCRIPTION`.

## What this is

**A fork of [`irilenia/baerhunter`](https://github.com/irilenia/baerhunter), not a replacement for it.** The modified files keep their original package paths, so the changes can be read as a well-scoped diff against upstream and merged back if the maintainers wish.

Three source files were revised across fourteen commits, a fourth was added, and the fifteenth turns the result back into an installable R package with generated documentation. What changed, and why, is recorded in `NEWS.md` and in the commit summaries under `commit_notes/`.

**In-scope source files:** `feature_file_editor.R`, `count_features.R`, `tpm_norm_flagging.R`, and the new `parameter_scout.R`.

## What this is not

**This fork does not carry `differential_expression.R`.** Upstream's package includes a wrapper around DESeq2 for differential expression testing; that file was out of scope for this project and is not part of this package. Anyone needing that step should use upstream baerhunter, or call DESeq2 directly on the count matrix `count_features()` produces.

This is a scope decision taken at the outset of the project, not a judgement about the function.

## Project goals

1. Test the package more thoroughly for bugs and improve its robustness.
2. Speed up the code where possible, without changing its scientific results.
3. Lay the groundwork for functional improvements, such as parameter selection.

## Parameter Scout

baerhunter requires the user to supply two important parameters, `low_coverage_cutoff` and `high_coverage_cutoff`, but offers no way to derive them. Since they are raw read depths, unscaled to library or genome, the coverage cutoffs that suit one dataset can be poor choices for another.

Parameter Scout is an optional advisory utility that reports the coverage percentiles of the intergenic regions, per BAM file, so the cutoffs can be read off the data rather than guessed. It runs as an optional step 0, before `feature_file_editor()`. It writes no annotation and sets no parameter: it produces a percentile table, a figure and two suggested cutoff pairs, and the user types the chosen values into the pipeline by hand.

**Full instructions:** [`documentation/parameter_scout_instructions.md`](documentation/parameter_scout_instructions.md).

## Repository structure

Two kinds of document live here, and they have different audiences.

**`documentation/`** is for anyone using the package. **`commit_notes/`** is a point-in-time record of how the package came to be in the state it is in, written for the readers of a thesis; it is not maintained documentation. Neither directory is installed with the package, so the scout instructions are read here on GitHub rather than through `help()`.

```
ursa-major/
├── README.md
├── LICENSE
├── NEWS.md                 # changes in this fork, by version
├── DESCRIPTION             # package metadata and dependencies
├── NAMESPACE               # generated; never edited by hand
├── .Rbuildignore           # what is not part of the installed package
├── R/                      # package source
│   ├── feature_file_editor.R
│   ├── count_features.R
│   ├── tpm_norm_flagging.R
│   ├── parameter_scout.R
│   └── baerhunter-package.R
├── man/                    # generated help pages
├── inst/
│   └── extdata/            # example annotation and BAM subsets, from upstream
├── documentation/          # user documentation
│   └── parameter_scout_instructions.md
└── commit_notes/           # project record: one summary per commit, and an index
    ├── commit1_summary.md
    ├── commit2_summary.md
    ├── commit3_summary.md
    ├── commit4_summary.md
    ├── commit5_summary.md
    ├── commit6_summary.md
    ├── commit7_summary.md
    ├── commit8_summary.md
    ├── commit9_summary.md
    ├── commit10_summary.md
    ├── commit11_summary.md
    ├── commit12_summary.md
    ├── commit13_summary.md
    ├── commit14_summary.md
    ├── commit15_summary.md
    ├── commit16_summary.md
    └── index.md
```

## Not yet included

**There is no vignette.** Upstream ships one, and reconciling it against the current API is deliberately held over rather than declined: the read-quality filter added in this fork did not exist when the vignette was written, so it has to be rebuilt against a package that installs before its output can be trusted.

**The exported functions carry no runnable examples.** Neither do upstream's, so this is an inherited gap rather than one introduced here. It is recorded as future work alongside the vignette.

## References

Ozuna, A., Liberto, D., Joyce, R. M., Arnvig, K. B., & Nobeli, I. (2020). baerhunter: An R package for the discovery and analysis of expressed non-coding regions in bacterial RNA-seq data. _Bioinformatics_, _36_(3), 966-969. [https://doi.org/10.1093/bioinformatics/btz643](https://doi.org/10.1093/bioinformatics/btz643)

Original package: [https://github.com/irilenia/baerhunter](https://github.com/irilenia/baerhunter)

## License

baerhunter is released under the MIT License. This repository preserves the original copyright (© 2019 irilenia) and licenses the present author's modifications under the same terms. See the `LICENSE` file for details.

## Authors

baerhunter was written by Alina Ozuna, and is maintained by Dr Irilenia Nobeli, with contributions from Jennifer Stiens. The modifications in this fork are by Thomas Quinn.

Email: [tquinn04@student.bbk.ac.uk](mailto:tquinn04@student.bbk.ac.uk) until 30 September 2026, and [thomquinn@gmail.com](mailto:thomquinn@gmail.com) from 1 October 2026.
