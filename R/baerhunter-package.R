#' @details
#' baerhunter predicts, annotates and filters expressed intergenic regions,
#' such as small RNAs and untranslated regions, from bacterial RNA-seq data,
#' using a coverage-based method that needs no reference set of known
#' features. It also supports the stages that follow, producing a count
#' matrix and TPM values, and it carries an optional advisory module that
#' reports the coverage percentiles from which the two coverage cutoffs can
#' be chosen.
#'
#' @section Fork provenance:
#' This is \code{ursa-major}, a fork of \code{irilenia/baerhunter} made as an
#' MSc Bioinformatics thesis project at Birkbeck, University of London. It
#' installs under the upstream package name, and its modified files keep
#' their original paths, so the changes read as a scoped difference against
#' upstream rather than as a new package. The changes are listed in
#' \code{NEWS.md} and are readable after installation with
#' \code{news(package = "baerhunter")}.
#'
#' This fork does not include \code{differential_expression()}. Upstream's
#' wrapper around DESeq2 was outside the scope of the work and is not part of
#' this package; run that step with upstream baerhunter, or call DESeq2
#' directly on the count matrix \code{count_features()} produces.
#'
#' @references
#' Ozuna, A., Liberto, D., Joyce, R. M., Arnvig, K. B., & Nobeli, I. (2020).
#' baerhunter: An R package for the discovery and analysis of expressed
#' non-coding regions in bacterial RNA-seq data. \emph{Bioinformatics},
#' \emph{36}(3), 966-969. \doi{10.1093/bioinformatics/btz643}
#'
#' @seealso
#' Upstream: \url{https://github.com/irilenia/baerhunter}.
#' This fork: \url{https://github.com/thomasjquinn/ursa-major}.
#'
#' @keywords internal
"_PACKAGE"
