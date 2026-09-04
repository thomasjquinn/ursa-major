# Commit 16: Installation Instructions and News Dependencies

## Commit Abstract: Correct the installation instructions so that a reader without the package's dependencies can follow them to a working install, and declare the two packages R's markdown news reader requires so that the documented `news()` call is honest on every machine.

### Commit Summary:
Fix installation instructions and news deps

### Commit Description:

```
docs(install): correct the README install block and declare news deps

The published installation instructions do not reach a working
install on a machine that lacks the package's dependencies. The
remotes install line was commented out, install_github does not
resolve Bioconductor dependencies without the repositories being
set, and the upgrade prompt recommends updating packages a reader
may wish to keep pinned.

Declare commonmark and xml2 in Suggests, and qualify the news()
recommendation in the package-level topic, so that the documented
route to the changelog is accurate whether or not those packages
are present.
```

---

## Why this commit exists

**Parent `669d3d5`, Commit 15**, landed 3 September 2026. That commit made the package installable; this one makes the instructions for installing it correct.

**Both defects were found by installing the package on a second machine for the first time**, within hours of Commit 15 being pushed. Neither is a fault in the code. Both are faults in what the package tells a user to do, and neither could have been found by `R CMD check`, which runs where everything is already present.

---

## C1 — the installation block

**Issue Summary:** The published block does not reach a working install on a machine lacking the dependencies, failing in three ways.

**Solution Summary:** Make both install lines executable, set the repositories so `install_github` can resolve Bioconductor, and decline package updates explicitly.

`remotes` is not installed and its install line was commented out. **Measured absent on three of three environments**, including one carrying `devtools` 2.5.2, which does not pull it in. `install_github` resolves against `getOption("repos")`, CRAN alone by default, and none of the seven Bioconductor imports is on CRAN; on a machine without them it reports `Skipping 7 packages not available` and installs a package that cannot load. The upgrade prompt invites a reader to update packages they may wish to keep pinned.

### README.md

```
# commit 15 (previous)
    # install.packages("remotes")
    remotes::install_github("thomasjquinn/ursa-major")
```

```
# commit 16 (this commit)
    install.packages(c("remotes", "BiocManager"))
    options(repos = BiocManager::repositories())
    remotes::install_github("thomasjquinn/ursa-major", upgrade = "never")
```

`BiocManager::repositories()` returns CRAN plus the Bioconductor repositories for the current release, so `install_github` resolves the seven itself. Listing them by hand was declined: a hard-coded list in a README goes stale and `DESCRIPTION` already carries the authoritative version.

The surrounding prose gains two sentences, replacing the existing sentence about requiring R 4.5 and a current Bioconductor installation, which is true but does not tell a reader what to do about it:

```
The package requires R 4.5 or later and Bioconductor 3.22 or later.
install_github resolves CRAN dependencies on its own but not
Bioconductor ones, which is why the install block above sets repos.
```

**Two further README edits.** `commit16_summary.md` is added to the `commit_notes/` tree listing, and the author line is shortened to `The modifications in this fork are by Thomas Quinn.`

---

## C2 — `Suggests` gains `commonmark` and `xml2`

**Issue Summary:** `news(package = "baerhunter")` fails without both. Neither is a base or recommended package, and neither arrives with any declared dependency of this package.

**Solution Summary:** Declare both in `Suggests`.

```
# commit 16 (this commit)
Suggests: 
    knitr,
    rmarkdown,
    RColorBrewer,
    commonmark,
    xml2,
    testthat (>= 3.0.0)
```

`Suggests` is the right field. The package's own code calls neither, so an `Imports` entry would be false as well as forcing both onto every user for a call most will never make.

---

## C3 — qualify the `news()` sentence

**Issue Summary:** `man/baerhunter-package.Rd` tells the reader the changes are readable with `news(package = "baerhunter")`. On a machine without the two packages that call fails, so the instruction is false for a subset of users.

**Solution Summary:** Qualify the sentence. Every edited line is a `#'` comment; `man/baerhunter-package.Rd` is regenerated rather than edited.

```
# commit 15 (previous) -- R/baerhunter-package.R
#' \code{NEWS.md} and are readable after installation with
#' \code{news(package = "baerhunter")}.
```

```
# commit 16 (this commit)
#' \code{NEWS.md}, readable after installation with
#' \code{news(package = "baerhunter")} where the \pkg{commonmark} and
#' \pkg{xml2} packages are available, and on the repository otherwise.
```

---

## Documents

`commit16_summary.md` and the updated `index.md` are filed in `commit_notes/`. `index.md` gains row 16 and records `669d3d5` against row 15, whose number was unknown when it was written.
