# ursa-major

**Improving the speed and usability of baerhunter for identifying unannotated expressed regions in bacterial transcriptomes**

> ⚠️ **WORK IN PROGRESS**
> This project is currently under active development. Some features may be incomplete.

## About

`ursa-major` is the working codename for an MSc Bioinformatics thesis project at Birkbeck, University of London, which aims to improve the performance and functionality of [baerhunter](https://github.com/irilenia/baerhunter).

baerhunter is an R package that uses a coverage-based method to predict, annotate, and filter unannotated expressed regions, such as small RNAs (sRNAs) and untranslated regions (UTRs), from bacterial RNA-seq data. This project revisits the package to make it faster, more robust, and easier to use, while preserving its existing scientific output.

The work is supervised by Dr Irilenia Nobeli at Birkbeck, University of London. So far, the updates have been validated against *Mycobacterium tuberculosis* H37Rv RNA-seq data, with additional paired-end testing on *Salmonella* Typhimurium. Additional tests are forthcoming.

## Project goals

The aims of this project are to:

1. Test the package more thoroughly for bugs and improve its robustness.
2. Speed up the code where possible, without changing its scientific results.
3. Lay the groundwork for functional improvements, such as parameter selection.

## What this repository contains

This repository is a *contribution overlay*, not a full copy of baerhunter. It holds only the files created or modified as part of this thesis:

- the R source files under active development, at their original package paths,
- the thesis project documentation.

All other parts of the baerhunter package (the unchanged source files, including `differential_expression.R`, plus the generated help pages, the vignette, example data, and package metadata) are not redistributed here.

Because it is an overlay rather than a complete package, this repository will not install with `R CMD INSTALL` on its own. Instructions for applying these changes onto an upstream checkout, and for rebuilding the package documentation, are given in `documentation/integration/integration.md`.

In-scope source files: `feature_file_editor.R`, `count_features.R`, and `tpm_norm_flagging.R`. `differential_expression.R` is not modified in this project and is recorded as future work.

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
│   └── tpm_norm_flagging.R
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
    └── integration/        # applied at the end of the project
        └── integration.md
```

## Integrating these changes upstream

The modified files keep their original package paths so they can be merged back into baerhunter as a small, well-scoped change. `documentation/integration/integration.md` records the upstream base the changes apply against, which files they replace, the dependency version requirements they introduce, and the steps to regenerate the package documentation.

## References

Ozuna, A., Liberto, D., Joyce, R. M., Arnvig, K. B., & Nobeli, I. (2020). baerhunter: An R package for the discovery and analysis of expressed non-coding regions in bacterial RNA-seq data. *Bioinformatics*, *36*(3), 966-969. https://doi.org/10.1093/bioinformatics/btz643

Original package: https://github.com/irilenia/baerhunter

## License

baerhunter is released under the MIT License. This repository preserves the original copyright (© 2019 irilenia) and licenses the present author's modifications under the same terms. See the [`LICENSE`](LICENSE) file for details.

## Author

Thomas Quinn, MSc Bioinformatics, Birkbeck, University of London. Email: tquinn04@student.bbk.ac.uk
