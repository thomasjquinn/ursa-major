# Integrating these changes into baerhunter

> ⚠️ **PLACEHOLDER — NOT YET APPLICABLE**
> This project is still in progress, so integration is not yet possible. This note is a forward-looking reminder of best practice for the end of the project, once the in-scope changes are complete and approved. The placeholders below (upstream base commit, final version numbers) are filled in at that point.

## Purpose

This repository is a contribution overlay: it carries only the files created or modified during the thesis, at their original package paths, rather than a full copy of baerhunter. This note records how those changes should be folded back into the upstream package as a small, self-contained update. It is a checklist of good practice, not a record of work already merged.

## Provenance and base

The changes are authored against upstream baerhunter (`irilenia/baerhunter`, package version 0.9.1.0000, `master`). Record the exact upstream base commit here once work begins:

- Upstream base commit: `________`

Good practice is to make the first commit of the working history the unmodified upstream baseline, so that every later change reads as a clean diff against it. That diff is the authoritative record of what is, and is not, the author's work.

## What gets integrated

Three source files replace their upstream namesakes and two are new; nothing else in `R/` is touched.

| Overlay file | Replaces in upstream |
| --- | --- |
| `R/feature_file_editor.R` | `R/feature_file_editor.R` |
| `R/count_features.R` | `R/count_features.R` |
| `R/tpm_norm_flagging.R` | `R/tpm_norm_flagging.R` |

| New overlay file | Upstream counterpart |
| --- | --- |
| `R/parameter_scout_paired_end.R` | none, added by this work |
| `R/parameter_scout_single_end.R` | none, added by this work |

The two Parameter Scout files are an optional advisory module: it reports the coverage percentiles of the intergenic regions so the two coverage cutoffs can be chosen from the data, and it writes no annotation and sets no parameter. It replaces nothing, and no existing file was modified to accommodate it, so it can be integrated independently of the three replacements above, or left out. **It is not yet package-ready**, which is the one thing to know before copying it across; see Regenerating package documentation below.

`differential_expression.R` is not modified and should be left as upstream. The package vignette gains a short addition for the GFF-cache helper, and `DESCRIPTION`, `NAMESPACE`, and the `man/` help pages are regenerated (see below).

## Behaviour-affecting changes

Some changes deliberately alter scientific output relative to the current upstream package, so integration is not byte-identical. The authoritative list lives in `NEWS.md` and the per-commit summaries; in outline:

- Read counting is made deterministic across Rsubread versions by setting `countMultiMappingReads` and `countReadPairs` explicitly. Multi-mapper counts at repetitive loci change relative to the modern Rsubread default.
- The feature-identifier regex is made delimiter-agnostic, so NCBI RefSeq style identifiers are parsed correctly. This is a correctness fix that affects the annotation output.
- The default `mapqFilter` is raised to 10 (on the supervisor's instruction). This value suits most aligners; Bowtie2 is the documented exception and needs its own value.

Integrate with these in view, and carry the matching `NEWS.md` entries across.

## Dependency requirements

The changes introduce a hard dependency on a recent enough Rsubread, and Parameter Scout adds a further set of packages. Add an `Imports` floor to `DESCRIPTION` at the documentation pass:

- Absolute minimum: `Rsubread (>= 2.4.3)`, the lowest version confirmed to carry `countReadPairs`.
- Recommended: `Rsubread (>= 2.16)`, per the dependency audit. This also covers the `largest_overlap = TRUE` with `fraction = TRUE` combination, which was silently miscounted before Rsubread 2.14.0.

The requirement is documented rather than enforced at runtime: on an older Rsubread the call fails with a loud, self-explaining error rather than silently changing results.

Parameter Scout needs `S4Vectors`, `IRanges`, `GenomicRanges`, `GenomicAlignments`, `Rsamtools`, `Seqinfo`, `ggplot2` and `grid`, of which the first five are already package dependencies. `grDevices`, `stats` and `utils` are used qualified and also need declaring, and `RColorBrewer` belongs in `Suggests`, being reached only through `requireNamespace()` with a base-R fallback. One point needs a decision rather than a copy: the `seqinfo` generics moved from `GenomeInfoDb` into their own `Seqinfo` package, and the module currently chooses between the two at run time, which a static `NAMESPACE` import cannot do. The package imports `Seqinfo` and takes its Bioconductor floor from there. The full measured list is in `issue068.md`.

## Regenerating package documentation

**Do the Parameter Scout package conversion first.** The two scout files belong in `R/` alongside the other three; the problem is not where they sit but how they are evaluated. They are written to be evaluated one at a time, which is how they were built and tested, and three properties do not survive being evaluated together into a single namespace, which is what `devtools::load_all()` and a real package build both do: the module attaches its packages at top level, which a package must not do; the two siblings define the same 54 objects with identical bodies, and R keeps only the last one evaluated, silently; and two of those duplicates carry `@export`, so `NAMESPACE` would take duplicate entries and roxygen would be asked to write two help pages for one topic. **`devtools::load_all()` on this overlay would already misbehave**, which is easy to reach for by habit. The conversion is specified in `issue068.md` and is not optional: step 1 below will not produce correct output until it is done.

After copying the five files onto an upstream checkout:

1. Run `devtools::document()` to regenerate `NAMESPACE` and `man/*.Rd` from the roxygen tags in the source files. The overlay does not ship these; they are rebuilt here.
2. Apply the vignette addition and update `NEWS.md`.
3. Run `devtools::check()`. It should pass, or show only the same notes as the upstream baseline. Expect the `ggplot2` aesthetic names to raise `no visible binding for global variable` on the first run; that is the ordinary non-standard-evaluation note and is answered with `utils::globalVariables()` or `.data$` pronouns, not by changing the plotting code.
4. Ship `parameter_scout_instructions.md` alongside the module. It is the user-facing guide, written for someone who has the two files and nothing else.

## Workflow reminders

Best practice when the changes are folded in:

- Work through feature branches and pull requests, one logical change per branch.
- Use Conventional Commits messages (`type(scope): description`), imperative mood, subject under about 72 characters.
- Keep the maintainer (`cre`) in `DESCRIPTION` as Dr Irilenia Nobeli.
- Use `NEWS.md` as the changelog, one bullet per change.
- Bump the version with `usethis::use_version("minor")` (0.9.1.0000 to 0.10.0) and tag the release when the work is merged.
- Locate edits by code content, not line number, since lines shift between changes.
- Keep LF line endings throughout.

## Out of scope

`differential_expression.R` is left for future work and is recorded as a recommendation for the next maintainer rather than implemented here.

*Parameter selection was listed here as future work until 10 August 2026. It is no longer: Parameter Scout is implemented, tested on three datasets, and integrated as described above. What remains of it is the package conversion, which is scoped work with a ticket rather than a recommendation.*
