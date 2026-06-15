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

Three source files replace their upstream namesakes; nothing else in `R/` is touched.

| Overlay file | Replaces in upstream |
| --- | --- |
| `R/feature_file_editor.R` | `R/feature_file_editor.R` |
| `R/count_features.R` | `R/count_features.R` |
| `R/tpm_norm_flagging.R` | `R/tpm_norm_flagging.R` |

`differential_expression.R` is not modified and should be left as upstream. The package vignette gains a short addition for the GFF-cache helper, and `DESCRIPTION`, `NAMESPACE`, and the `man/` help pages are regenerated (see below).

## Behaviour-affecting changes

Some changes deliberately alter scientific output relative to the current upstream package, so integration is not byte-identical. The authoritative list lives in `NEWS.md` and the per-commit summaries; in outline:

- Read counting is made deterministic across Rsubread versions by setting `countMultiMappingReads` and `countReadPairs` explicitly. Multi-mapper counts at repetitive loci change relative to the modern Rsubread default.
- The feature-identifier regex is made delimiter-agnostic, so NCBI RefSeq style identifiers are parsed correctly. This is a correctness fix that affects the annotation output.
- The default `mapqFilter` is raised to 10 (on the supervisor's instruction). This value suits most aligners; Bowtie2 is the documented exception and needs its own value.

Integrate with these in view, and carry the matching `NEWS.md` entries across.

## Dependency requirements

The changes introduce a hard dependency on a recent enough Rsubread. Add an `Imports` floor to `DESCRIPTION` at the documentation pass:

- Absolute minimum: `Rsubread (>= 2.4.3)`, the lowest version confirmed to carry `countReadPairs`.
- Recommended: `Rsubread (>= 2.16)`, per the dependency audit. This also covers the `largest_overlap = TRUE` with `fraction = TRUE` combination, which was silently miscounted before Rsubread 2.14.0.

The requirement is documented rather than enforced at runtime: on an older Rsubread the call fails with a loud, self-explaining error rather than silently changing results.

## Regenerating package documentation

After copying the three files onto an upstream checkout:

1. Run `devtools::document()` to regenerate `NAMESPACE` and `man/*.Rd` from the roxygen tags in the source files. The overlay does not ship these; they are rebuilt here.
2. Apply the vignette addition and update `NEWS.md`.
3. Run `devtools::check()`. It should pass, or show only the same notes as the upstream baseline.

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

`differential_expression.R` is left for future work, as is the parameter-selection functionality. These are recorded as recommendations for the next maintainer rather than implemented here.
