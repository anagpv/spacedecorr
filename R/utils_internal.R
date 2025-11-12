#' @keywords internal
#' @noRd
clipRes <- function(vec, n){
  b <- sqrt(n)
  vec[vec < -b] <- -b
  vec[vec >  b] <-  b
  vec
}

#' @keywords internal
#' @noRd
getResiduals <- function(target, data, formula, assay_matrix, num_cells, family, highdim, verbose){
  if (isTRUE(verbose)) message(target)
  data$y <- assay_matrix[, target]

  mod <- if (isTRUE(highdim)) {
    mgcv::bam(formula, family = family, data = data,
              method = "fREML", discrete = TRUE, gc.level = 2, select = TRUE)
  } else {
    mgcv::gam(formula, family = family, data = data,
              method = "REML", select = TRUE)
  }

  res   <- clipRes(mgcv::residuals.gam(mod, type = "pearson"), num_cells)
  space <- mgcv::predict.gam(mod, newdata = data, type = "terms")

  kp <- try(mgcv::k.check(mod), silent = TRUE)
  kpvalue <- if (!inherits(kp, "try-error") && "p-value" %in% colnames(kp)) {
    suppressWarnings(min(kp[, "p-value"], na.rm = TRUE))
  } else NA_real_

  list(
    target    = target,
    residuals = res,
    splines   = space,
    pvalue    = as.numeric(kpvalue)
  )
}
