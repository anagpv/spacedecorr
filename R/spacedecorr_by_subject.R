#' Run SpaceDecorr separately for each subject
#'
#' Splits cells by subject and runs \code{spacedecorr()} within each subject.
#'
#' @param assay_matrix Numeric matrix, ideally cells x genes.
#' @param metadata Data.frame with one row per cell.
#' @param subject_col Column name in metadata indicating subject ID.
#' @param return One of "combined", "list", or "both".
#' @param min_cells Minimum number of cells required per subject; smaller subjects are skipped.
#' @param ... Additional arguments passed to \code{spacedecorr()}.
#'
#' @return Depending on \code{return}:
#' \itemize{
#'   \item "combined": list with \code{residuals} (matrix) and \code{metadata} (data.frame)
#'   \item "list": named list of \code{spacedecorr_fit} objects
#'   \item "both": list with \code{combined} and \code{fits}
#' }
#'
#' @export
#'
spacedecorr_by_subject <- function(
    assay_matrix,
    metadata,
    subject_col = "subject",
    return = c("combined", "list", "both"),
    min_cells = 50,
    ...
) {
  return <- match.arg(return)

  assay_matrix <- as.matrix(assay_matrix)
  metadata <- as.data.frame(metadata)

  if (!subject_col %in% names(metadata)) {
    stop("metadata must contain subject_col: ", subject_col)
  }

  subjects <- unique(metadata[[subject_col]])
  subjects <- subjects[!is.na(subjects)]

  fits <- vector("list", length(subjects))
  names(fits) <- as.character(subjects)

  for (sj in subjects) {
    idx <- which(metadata[[subject_col]] == sj)
    if (length(idx) < min_cells) next

    fit <- spacedecorr(
      assay_matrix = assay_matrix[idx, , drop = FALSE],
      metadata = metadata[idx, , drop = FALSE],
      ...
    )

    fits[[as.character(sj)]] <- fit
  }

  fits <- Filter(Negate(is.null), fits)

  if (return == "list") return(fits)

  residuals <- do.call(rbind, lapply(fits, `[[`, "residuals"))
  meta_out  <- do.call(rbind, lapply(fits, `[[`, "metadata"))

  combined <- list(
    residuals = residuals,
    metadata  = meta_out
  )

  if (return == "both") return(list(combined = combined, fits = fits))
  combined
}
