#' Spatial de-correlation and covariate-adjusted residuals for single-cell data
#'
#' Models each gene or protein independently to remove spatially structured variation
#' and adjust for covariates, returning Pearson residuals that represent decorrelated
#' expression values. Diagnostic summaries of the spatial adjustment are also provided.
#'
#' @section Overview:
#' \code{spacedecorr()} takes an expression matrix with cells in rows and features in
#' columns, together with metadata containing spatial coordinates and, optionally,
#' library size and other covariates. Advanced users may override the default model
#' specification via a custom formula.
#'
#' @param assay_matrix Numeric matrix (n x p). Rows are cells/spots; columns are features.
#' @param metadata Data frame with at least the location and library-size columns.
#'   Can include any additional covariates to adjust for.
#' @param loc_cols Character length-2. Names of the x and y columns in \code{metadata}.
#'   Default: \code{c("sdimx","sdimy")}.
#' @param libsize_col Column name for library size in metadata, or NULL.
#'   If NULL and family="nb", library size is computed as rowSums(assay_matrix)
#'   and used as offset(log(libsize)).
#' @param covariates Either \code{NULL} (no extra covariates), \code{"all"} (use all
#'   columns in \code{metadata} except location/library-size), or a character vector of
#'   column names to include.
#' @param kprop Numeric in (0,1]. Proportion of cells used to set the spline basis
#'   dimension \code{k}. Ignored if \code{k} is provided. Default: \code{NULL}.
#' @param k Integer. Basis dimension for the spatial smooth (overrides \code{kprop}).
#'   Default: \code{NULL}.
#' @param m Integer. Derivative order for Duchon splines when \code{basis} is \code{"d"} or \code{"ds"}.
#'   Ignored otherwise. Default: \code{2}.
#' @param nCores Integer. Number of parallel workers. Use 1 for sequential. Default: \code{1}.
#' @param family Either \code{"nb"} (negative binomial; default), \code{"gaussian"},
#'   or a \link[stats]{family} object. For \code{"nb"}, \code{offset(log(libsize))} is added.
#' @param basis Character. Spatial smoother basis for \code{mgcv::s}: one of
#'   \code{"tp"}, \code{"ts"}, \code{"d"}, \code{"ds"}, \code{"te"}. Default: \code{"ts"}.
#' @param verbose Logical. Print feature names as they are processed. Default: \code{FALSE}.
#' @param formula_override Optional formula to fully override the internally constructed model.
#'   Must reference standardized column names \code{libsize}, \code{sdimx}, \code{sdimy}, and any
#'   covariates present in \code{metadata}.
#' @param use_bam Logical; if TRUE uses \code{mgcv::bam()} (faster for large n),
#'   if FALSE uses \code{mgcv::gam()}. If NULL, defaults to TRUE when n >= 5000.
#' @param return_splines Logical; if TRUE, return the estimated spatial smooth term per gene.
#'   Can be large for many genes; set FALSE to save memory.
#'
#' @return A named list with:
#' \itemize{
#'   \item \code{Residuals}: data.frame (n x p) of Pearson residuals (one column per feature).
#'   \item \code{Splines}:   data.frame (n x p) of the spatial-term predictions (first term column).
#'   \item \code{k.check.p}: data.frame (1 x p) with a single numeric p-value per feature.
#'   \item \code{metadata}:  data.frame used for model fitting (standardized columns included).
#' }
#'
#' @details
#' If \code{nrow(assay_matrix) != nrow(metadata)} but \code{ncol(assay_matrix) == nrow(metadata)},
#' the matrix is transposed. If both \code{assay_matrix} and \code{metadata} have rownames,
#' they must match; rows are reordered to align.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' # simulate toy spatial data
#' n <- 500   # number of cells
#' p <- 5     # number of genes
#' metadata <- data.frame(x = runif(n), y = runif(n))
#' Sigma <- telefit::maternCov(as.matrix(dist(coords)), smoothness = 0.5, range = 0.25)
#' Delta <- diag(1, nrow = p)
#'
#' # Simulate Matrix-Variate
#' tmp <- MASS::mvrnorm(n = 1, mu = rep(0, num_cells * p), Sigma = Sigma %x% Delta)
#' tmp <- as.matrix(matrix(tmp, ncol = p, byrow = T))
#' gamma <- qgamma(pnorm(tmp), shape = 1, rate = 1)
#'
#' # Counts
#' Y <- matrix(rpois(num_cells * p, lambda = as.numeric(gamma)),
#'             nrow = num_cells, ncol = p)
#' colnames(Y) <- paste0("gene", 1:p)
#'
#' # Fit spacedecorr
#' fit <- spacedecorr(
#'   assay_matrix = Y,
#'   metadata = coords,
#'   loc_cols = c("x", "y"),
#'   libsize_col = NA,
#'   family = "nb",
#'   k = 200,
#'   return_splines = TRUE
#' )
#' }
#' @importFrom stats as.formula gaussian
#' @importFrom mgcv gam bam nb predict.gam residuals.gam k.check
#' @export


spacedecorr <- function(assay_matrix,
                        metadata,
                        loc_cols    = c("sdimx","sdimy"),
                        libsize_col = NULL,
                        covariates  = NULL,
                        kprop = NULL,
                        k = 200,
                        m = 2,
                        nCores = 1,
                        family = "nb",
                        basis  = "ts",
                        verbose = FALSE,
                        formula_override = NULL,
                        use_bam = NULL,
                        return_splines = TRUE) {

  assay_matrix <- as.matrix(assay_matrix)
  metadata <- as.data.frame(metadata)

  if (length(loc_cols) != 2) stop("loc_cols must be length-2: c(x_col, y_col).")
  if (!all(loc_cols %in% names(metadata)))
    stop("metadata must contain location columns: ", paste(loc_cols, collapse = ", "))
  xcol <- loc_cols[1]; ycol <- loc_cols[2]
  rn_meta <- rownames(metadata)

  # Standardized columns
  metadata$sdimx   <- metadata[[xcol]]
  metadata$sdimy   <- metadata[[ycol]]

  # library size
  use_libsize <- FALSE
  if (identical(family, "nb")){
    if(is.na(libsize_col)){
      use_libsize <- FALSE
      metadata$libsize <- NA
    } else if (is.null(libsize_col) | (!libsize_col %in% names(metadata))) {
      metadata$libsize <- rowSums(assay_matrix)
      use_libsize <- TRUE
    } else {
      ls <- metadata[[libsize_col]]
    }
  }


  if (nrow(assay_matrix) != nrow(metadata)) {
    if (ncol(assay_matrix) == nrow(metadata)) {
      assay_matrix <- t(assay_matrix)
    } else {
      stop("assay_matrix must be n x p with n = nrow(metadata).")
    }
  }
  rna <- rownames(assay_matrix)
  if (!is.null(rna) && !is.null(rn_meta)) {
    if (!setequal(rna, rn_meta)) stop("Row names of assay_matrix and metadata do not match.")
    assay_matrix <- assay_matrix[rn_meta, , drop = FALSE]
  }

  # family ---
  family <- if (identical(family, "gaussian")) stats::gaussian() else mgcv::nb()

  #  k / kprop ---
  n <- nrow(assay_matrix)
  ksplines <- if (!is.null(k)) as.integer(k) else as.integer(round((if (is.null(kprop)) 0.1 else kprop) * n))

  # smoother ---
  if (basis %in% c("d","ds")) {
    smoother <- paste0('s(sdimx, sdimy, bs="', basis, '", k=', ksplines, ', xt=list(m=', m, '))')
  } else {
    smoother <- paste0('s(sdimx, sdimy, bs="', basis, '", k=', ksplines, ')')
  }

  # choose covariates ---
  reserved_std <- c("sdimx", "sdimy", "libsize")  # always exclude these from covariates
  reserved_orig <- unique(c(loc_cols, libsize_col))
  reserved_orig <- reserved_orig[!is.na(reserved_orig)]  # handles libsize_col = NULL

  if (is.null(covariates)){
    covnames <- character(0)
  }else if (identical(covariates, "all")){
    covnames <- setdiff(names(metadata), c(reserved_std, reserved_orig))
  }else if (is.character(covariates)){
    missing_cov <- setdiff(covariates, names(metadata))
    if (length(missing_cov)) warning("Dropping missing covariates: ", paste(missing_cov, collapse = ", "))
    covnames <- setdiff(intersect(covariates, names(metadata)), reserved_std)
  }else{
    stop("covariates must be NULL, 'all', or a character vector.")
  }

  # build model data & formula ---
  dat <- metadata[, unique(c("sdimx", "sdimy", "libsize", covnames)), drop = FALSE]
  if (!use_libsize) dat$libsize <- NULL

  if (!is.null(formula_override)) {
    form <- stats::as.formula(formula_override)
  } else {
    rhs <- paste(c(covnames, smoother), collapse = " + ")
    if (identical(family$family, "gaussian") || !use_libsize) {
      form <- stats::as.formula(paste0("y ~ ", rhs))
    } else {
      form <- stats::as.formula(paste0("y ~ offset(log(libsize)) + ", rhs))
    }
  }

  # labels & dims ---
  if (is.null(colnames(assay_matrix))) colnames(assay_matrix) <- paste0("X", seq_len(ncol(assay_matrix)))
  targets   <- colnames(assay_matrix)
  num_cells <- nrow(assay_matrix)
  if(is.null(use_bam)) use_bam <- (num_cells >= 5000)

  # run ---
  if (nCores > 1 && .Platform$OS.type != "windows") {
    pieces <- parallel::mclapply(
      targets,
      getResiduals,
      data = dat, formula = form, assay_matrix = assay_matrix,
      num_cells = num_cells, family = family, use_bam = use_bam, verbose = verbose,
      return_splines = return_splines,
      mc.cores = nCores
    )

  } else if (nCores > 1) {
    cl <- parallel::makeCluster(getOption("cl.cores", nCores))
    parallel::clusterEvalQ(cl, { library(mgcv) })
    pieces <- parallel::parLapply(
      cl,targets, getResiduals,
      data = dat, formula = form, assay_matrix = assay_matrix,
      num_cells = num_cells, family = family, use_bam = use_bam, verbose = verbose,
      return_splines = return_splines
    )
    suppressWarnings(parallel::stopCluster(cl))
  } else {
    pieces <- lapply(
      targets,getResiduals,
      data = dat, formula = form, assay_matrix = assay_matrix,
      num_cells = num_cells, family = family, use_bam = use_bam, verbose = verbose,
      return_splines = return_splines
    )

  }

  # combine ---
  Residuals <- do.call(cbind, lapply(pieces, `[[`, "residuals"))
  colnames(Residuals) <- vapply(pieces, `[[`, character(1), "target")
  Residuals <- as.data.frame(Residuals)

  Splines <- NULL
  if (isTRUE(return_splines)){
    Splines <- do.call(cbind, lapply(pieces, `[[`, "splines"))
    colnames(Splines) <- vapply(pieces, `[[`, character(1), "target")
    Splines <- as.data.frame(Splines)
  }

  Pvalues <- as.data.frame(t(vapply(pieces, function(x) x$pvalue, numeric(1))))
  names(Pvalues) <- vapply(pieces, `[[`, "", "target")
  rownames(Pvalues) <- "pvalue"

  structure(
    list(
      residuals = Residuals,
      splines = Splines,
      k_check_p = Pvalues,
      metadata = dat,
      call = match.call(),
      config = list(family = family$family, basis = basis, k = ksplines, m = m, use_bam = use_bam,
                    covariates = covariates, loc_cols = loc_cols, libsize_col = libsize_col)
    ),
    class = "spacedecorr_fit"
  )
}

#' @noRd
getResiduals <- function(target, data, formula, assay_matrix,
                         num_cells, family, use_bam, verbose,return_splines){
  if (isTRUE(verbose)) message(target)
  data$y <- assay_matrix[, target]

  mod <- if (isTRUE(use_bam)) {
    mgcv::bam(formula, family = family, data = data,
              method = "fREML", discrete = TRUE, gc.level = 2, select = TRUE)
  } else {
    mgcv::gam(formula, family = family, data = data,
              method = "REML", select = TRUE)
  }

  # Residuals
  res   <- clip_residuals_(mgcv::residuals.gam(mod, type = "pearson"), num_cells)

  # Splines
  space <- NULL
  if (isTRUE(return_splines)) {
    space <- mgcv::predict.gam(mod, data, type = "terms")[, "s(sdimx,sdimy)"]
  }

  # Test for eigen dimension
  kp <- try(mgcv::k.check(mod), silent = TRUE)
  kpvalue <- if (!inherits(kp, "try-error") && "p-value" %in% colnames(kp)) {
    suppressWarnings(min(kp[, "p-value"], na.rm = TRUE))
  } else NA_real_

  list(
    target    = target,
    residuals = as.numeric(res),
    splines   = as.numeric(space),
    pvalue    = as.numeric(kpvalue)
  )
}

#' @noRd
clip_residuals_ <- function(vec, n) {
  b <- sqrt(n)
  vec[vec < -b] <- -b
  vec[vec >  b] <-  b
  vec
}

