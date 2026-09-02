# Commit 12: Repository Configuration

## Commit Abstract: Configure Git's handling of binary data and build output, so that the alignment files bundled in the next commit are stored byte-exact rather than line-ending normalised, and so that generated build artefacts are not tracked.

### Commit Number:

`78c263edff1b2c2b9d72a90985eb8d73655e8727`, short form `78c263e`. Parent `325c327`. Landed 31 August 2026.

### Commit Summary:
Configure Git handling for data and build output

### Commit Description:

```
chore(repo): add binary and build-output rules

Mark *.bam and *.bai binary so Git applies no line-ending
conversion to the bundled alignment data added in the next
commit. A normalised BAM does not error; it produces different
coverage. Pin LF explicitly for *.gff3, *.txt, *.Rd and *.Rmd.

Add package build and check output to .gitignore: inst/doc/,
doc/, Meta/, the knitr intermediates, and output_BH_5_10/.

No existing rule is removed from either file.

Part of the repository file structure work.
```

### Commit Details:

Two files changed, 27 insertions and 1 deletion: 13 insertions and 1 deletion in `.gitattributes`, 14 insertions in `.gitignore`.

---

## Mark alignment data binary and pin line endings for the remaining text formats

**Issue Summary:** `.gitattributes` sets `* text=auto eol=lf` and names only `*.R` and `*.md` explicitly, so the alignment files about to be bundled would have their type inferred rather than declared, and the annotation and plain-text files would rely on auto-detection.

**Solution Summary:** Declare `*.bam` and `*.bai` binary alongside the existing `*.rds`, and pin LF explicitly for `*.Rd`, `*.Rmd`, `*.gff3` and `*.txt`, grouping the additions under comment headings that say what each group is for.

**Note:** This is a data-integrity fix rather than a style change. Auto-detection would probably have classified a BAM correctly, but "probably" is not a property to rest bundled alignment data on: a normalised BAM does not error, it produces different coverage, and nothing downstream would report it.

### .gitattributes

The header block and the first three rules are unchanged and are shown for context. The additions are the three commented groups below them.

```
# commit 11 (previous)
# lines 1-14
# .gitattributes
# Line-ending and binary handling for the ursa-major repository.
# Pinned here so behaviour is identical for every contributor, regardless of
# each machine's local Git configuration (this overrides core.autocrlf).

# Default: auto-detect text vs binary; store and check out all text as LF.
* text=auto eol=lf

# Be explicit about the text file types in this repository.
*.R    text eol=lf
*.md   text eol=lf

# Snapshot objects are binary: never apply line-ending conversion.
*.rds  binary
```

```
# commit 12 (this commit)
# lines 1-26
# .gitattributes
# Line-ending and binary handling for the ursa-major repository.
# Pinned here so behaviour is identical for every contributor, regardless of
# each machine's local Git configuration (this overrides core.autocrlf).

# Default: auto-detect text vs binary; store and check out all text as LF.
* text=auto eol=lf

# Be explicit about the text file types in this repository.
*.R    text eol=lf
*.md   text eol=lf

# Generated package documentation.
*.Rd   text eol=lf
*.Rmd  text eol=lf

# Bundled data: pin the endings rather than let them be inferred.
*.gff3 text eol=lf
*.txt  text eol=lf

# Binary: never apply line-ending conversion. Alignment files and their
# indexes are compressed, and a normalised BAM does not error, it produces
# different coverage.
*.rds  binary
*.bam  binary
*.bai  binary
```

**The single deletion is a comment, not a rule.** The line `# Snapshot objects are binary: never apply line-ending conversion.` is removed and its content absorbed into the fuller three-line note above the binary block, which now covers the alignment files as well. `*.rds binary` itself is untouched, and no rule is removed from either file.

---

## Ignore package build and check output

**Issue Summary:** `.gitignore` covers R session data, RStudio files, the `R CMD build` tarball and `.Rcheck` directory, but not the directories a package build and a vignette build leave behind, so building the package in place would offer generated artefacts for staging.

**Solution Summary:** Add the build output directories, the knitr and rmarkdown intermediates, and the output directory a vignette run produces when it is run in place rather than in `tempdir()`.

**Note:** This is anticipatory rather than corrective: none of these paths exists in the repository yet, and each becomes possible in a later commit. Adding them now costs nothing and removes the chance of a generated file being committed by accident during the packaging work.

### .gitignore

Inserted after the existing `R CMD build / check output` block and before `# Script backups and scratch files`; nothing already in the file changes.

```
# commit 11 (previous)
# lines 11-15
# R CMD build / check output
/*.tar.gz
/*.Rcheck/
*-Ex.R

```

```
# commit 12 (this commit)
# lines 11-29
# R CMD build / check output
/*.tar.gz
/*.Rcheck/
*-Ex.R
inst/doc/
doc/
Meta/

# knitr and rmarkdown intermediates
vignettes/*.html
vignettes/*.pdf
vignettes/*_cache/
vignettes/*_files/
*.knit.md
*.utf8.md

# Vignette run output, if run in place rather than in tempdir()
output_BH_5_10/

```
