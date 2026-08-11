# ursa-major

**Improving the speed and usability of the baerhunter software for identifying unannotated expressed regions in bacterial transcriptomes**

> ⚠️ **WORK IN PROGRESS**  
> This project is currently under active development. Some features may be incomplete.

## About

`ursa-major` is the working codename for an MSc Bioinformatics thesis project at Birkbeck, University of London, which aims to improve the performance and functionality of [baerhunter](https://github.com/irilenia/baerhunter).

baerhunter is an R package that uses a coverage-based method to predict, annotate, and filter unannotated expressed regions, such as small RNAs (sRNAs) and untranslated regions (UTRs), from bacterial RNA-seq data. This project revisits the package to make it faster, more robust, and easier to use, while preserving its existing scientific output.

The work is supervised by Dr Irilenia Nobeli at Birkbeck, University of London. So far, the updates have been validated against _Mycobacterium tuberculosis_ H37Rv RNA-seq data, with additional paired-end testing on _Salmonella_ Typhimurium. The improved pipeline has since been applied to _Mycobacterium bovis_ AF2122/97 RNA-seq data as a biological test case. Runtime has additionally been measured on paired-end _M. tuberculosis_ H37Rv RNA-seq (run ERR2103718, from E-MTAB-6011 / PRJEB65014). Additional tests are forthcoming.

## Project goals

The aims of this project are to:

1. Test the package more thoroughly for bugs and improve its robustness.
2. Speed up the code where possible, without changing its scientific results.
3. Lay the groundwork for functional improvements, such as parameter selection.

## Parameter Scout

baerhunter requires the user to supply two important parameters: `low_coverage_cutoff` and `high_coverage_cutoff`, but it offers no way to derive them. Since they are raw read depths, unscaled to library or genome, the coverage cutoffs that suit one dataset can be poor choices for another.

Parameter Scout is an optional advisory utility that reports the coverage percentiles of the intergenic regions, per BAM file, so the cutoffs can be read off the data rather than guessed. It is run as an optional step 0, before `feature_file_editor()`. It writes no annotation and sets no parameter: it produces a percentile table, a figure and two suggested cutoff pairs, and the user types the chosen values into the pipeline by hand.

Parameter Scout ships as two files, one for paired-end and one for single-end reads, which differ only in the name of the function they export and in one default. The method is adapted from Sivasankaran (2024) and rebuilt to work with the updated baerhunter code base. It has been tested mainly on three mycobacteria datasets: six _M. bovis_ AF2122/97 BAMs, three _M. tuberculosis_ H37Rv BAMs from E-MTAB-6011, and six single-end _M. tuberculosis_ H37Rv BAMs from E-MTAB-1616 (Cortes et al., 2013).

Further details and complete instructions for Parameter Scout can be found here: `documentation/parameter_scout_instructions.md`.

## What this repository contains

This repository is a _contribution overlay_, not a full copy of baerhunter. It holds only the files created or modified as part of this thesis:

* the R source files under active development, at their original package paths,
* the thesis project documentation.

All other parts of the baerhunter package (the unchanged source files, including `differential_expression.R`, plus the generated help pages, the vignette, example data, and package metadata) are not redistributed here.

Because it is an overlay rather than a complete package, this repository will not install with `R CMD INSTALL` on its own. Instructions for applying these changes onto an upstream checkout, and for rebuilding the package documentation, are given in `documentation/integration/integration.md`.

In-scope baerhunter source files: `feature_file_editor.R`, `count_features.R`, and `tpm_norm_flagging.R`.

In-scope Parameter Scout source files: `parameter_scout_paired_end.R` and `parameter_scout_single_end.R`. Both are new files, added at commit 11; no existing file was modified to accommodate them.

## Repository structure



```
ursa-major/
├── README.md
├── LICENSE
├── NEWS.md                # summary of changes
├── .gitignore
├── .gitattributes         # enforce LF line endings; .rds treated as binary
├── R/                     # in-scope source files (modified)
│   ├── feature_file_editor.R
│   ├── count_features.R
│   ├── tpm_norm_flagging.R
│   ├── parameter_scout_paired_end.R
│   └── parameter_scout_single_end.R
└── documentation/         # project notes and records
    ├── notes.md
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
    ├── parameter_scout_instructions.md
    └── integration/        # applied at the end of the project
        └── integration.md
```



## Integrating these changes upstream

The modified files keep their original package paths so they can be merged back into baerhunter as a small, well-scoped change. `documentation/integration/integration.md` records the upstream base the changes apply against, which files they replace, the dependency version requirements they introduce, and the steps to regenerate the package documentation.

The two Parameter Scout files keep their package paths like the rest, but they are written to be evaluated one at a time, with `source()`, which is how they are used and tested here. Evaluating both into a single namespace, as `devtools::load_all()` and a package build do, needs a conversion first: that is a separate, planned piece of work and is recorded in the project documentation.

## References

Cortes, T., Schubert, O. T., Rose, G., Arnvig, K. B., Comas, I., Aebersold, R., & Young, D. B. (2013). Genome-wide mapping of transcriptional start sites defines an extensive leaderless transcriptome in _Mycobacterium tuberculosis_. _Cell Reports_, _5_(4), 1121-1131. [https://doi.org/10.1016/j.celrep.2013.10.031](https://doi.org/10.1016/j.celrep.2013.10.031)

Ozuna, A., Liberto, D., Joyce, R. M., Arnvig, K. B., & Nobeli, I. (2020). baerhunter: An R package for the discovery and analysis of expressed non-coding regions in bacterial RNA-seq data. _Bioinformatics_, _36_(3), 966-969. [https://doi.org/10.1093/bioinformatics/btz643](https://doi.org/10.1093/bioinformatics/btz643)

Sivasankaran, A. (2024). _Detection and annotation of non-coding RNAs in bacteria using baerhunter: Developing a user-friendly RShiny application with intergenic coverage analysis for informed parameter selection_ [Unpublished MSc dissertation]. Birkbeck, University of London.

Original package: [https://github.com/irilenia/baerhunter](https://github.com/irilenia/baerhunter)

## License

baerhunter is released under the MIT License. This repository preserves the original copyright (© 2019 irilenia) and licenses the present author's modifications under the same terms. See the `LICENSE` file for details.

## Author

Thomas Quinn, MSc Bioinformatics, Birkbeck, University of London.

Email: [tquinn04@student.bbk.ac.uk](mailto:tquinn04@student.bbk.ac.uk) until 30 September 2026, and [thomquinn@gmail.com](mailto:thomquinn@gmail.com) from 1 October 2026.
