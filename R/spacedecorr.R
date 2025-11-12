#' Spatial de-correlation and covariate-adjusted residuals for single-cell data
#'
#' Fits a GAM/BAM per feature (gene/protein) to adjust for library size,
#' optional covariates, and a 2D spatial smooth. Returns Pearson residuals,
#' the spatial term, and one scalar p-value per feature (from \code{mgcv::k.check}).
#'
#' @section Overview:
#' \code{spacedecorr()} expects an expression matrix with cells in rows (n) and
#' features in columns (p), and a \code{metadata} data.frame that contains x/y
#' coordinates and a library-size column. You can (i) select which covariates to
#' adjust (\code{covariates}), or (ii) fully override the model with
#' \code{formula_override} (advanced).
#'
#' @param assay_matrix Numeric matrix (n x p). Rows are cells/spots; columns are features.
#' @param metadata Data frame with at least the location and library-size columns.
#'   Can include any additional covariates to adjust for.
#' @param loc_cols Character length-2. Names of the x and y columns in \code{metadata}.
#'   Default: \code{c("sdimx","sdimy")}.
#' @param libsize_col Character scalar. Name of the library-size column in \code{metadata}.
#'   If missing, it is computed as \code{rowSums(assay_matrix)}. Default: \code{"libsize"}.
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
#'   covariates present in \code{metadata}. Use only if you know what you're doing.
#'
#' @returns A named list with:
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
#' P-values are currently the minimum \code{k.check} p-value across smooths for each feature
#' (fast, simple heuristic). You can swap this strategy inside \code{getResiduals()} if desired.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 100; p <- 5
#' X <- matrix(rpois(n*p, 5), n, p)
#' meta <- data.frame(
#'   sdimx = runif(n), sdimy = runif(n),
#'   libsize = rowSums(X),
#'   batch = sample(letters[1:3], n, TRUE)
#' )
#'
#' out <- spacedecorr(
#'   assay_matrix = X,
#'   metadata     = meta,
#'   covariates   = c("batch"),
#'   family       = "nb",
#'   basis        = "ts",
#'   nCores       = 1
#' )
#' str(out)
#'
#' # Using all covariates (except coords/libsize):
#' out2 <- spacedecorr(X, meta, covariates = "all", family = "nb")
#'
#' # Advanced: full formula override (must reference libsize/sdimx/sdimy):
#' out3 <- spacedecorr(
#'   X, meta, covariates = NULL, family = "nb",
#'   formula_override = y ~ offset(log(libsize)) + batch + s(sdimx, sdimy, bs = "ts", k = 80)
#' )
#' }
#'
#' @importFrom stats as.formula gaussian
#' @importFrom mgcv gam bam nb predict.gam residuals.gam k.check
#' @export


spacedecorr <- function(assay_matrix,
                        metadata,
                        loc_cols    = c("sdimx","sdimy"),  # names of x,y cols in metadata
                        libsize_col = "libsize",           # name of library size col in metadata
                        covariates  = NULL,        # NULL | "all" | character()
                        kprop = NULL,
                        k = NULL,
                        m = 2,
                        nCores = 1,
                        family = "nb",             # "gaussian" or "nb"
                        basis  = "ts",             # "tp","ts","d","ds","te"
                        verbose = FALSE,
                        formula_override = NULL) { # optional: full mgcv formula (uses standardized names)

  assay_matrix <- as.matrix(assay_matrix)
  metadata <- as.data.frame(metadata)

  # ---- check/standardize required columns ------------------------------------
  if (length(loc_cols) != 2) stop("loc_cols must be length-2: c(x_col, y_col).")
  if (!all(loc_cols %in% names(metadata)))
    stop("metadata must contain location columns: ", paste(loc_cols, collapse = ", "))
  xcol <- loc_cols[1]; ycol <- loc_cols[2]

  # if libsize missing, compute it; then standardize names
  if (!libsize_col %in% names(metadata)) {
    metadata[[libsize_col]] <- rowSums(assay_matrix)
  }

  # Keep rownames for alignment
  rn_meta <- rownames(metadata)
  # Make standardized columns used by the model
  metadata$sdimx   <- metadata[[xcol]]
  metadata$sdimy   <- metadata[[ycol]]
  metadata$libsize <- metadata[[libsize_col]]

  # ---- orientation/alignment --------------------------------------------------
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

  # ---- family -----------------------------------------------------------------
  family <- if (identical(family, "gaussian")) stats::gaussian() else mgcv::nb()

  # ---- k / kprop --------------------------------------------------------------
  n <- nrow(assay_matrix)
  ksplines <- if (!is.null(k)) as.integer(k) else as.integer(round((if (is.null(kprop)) 0.1 else kprop) * n))

  # ---- smoother ---------------------------------------------------------------
  if (basis %in% c("d","ds")) {
    smoother <- paste0('s(sdimx, sdimy, bs="', basis, '", k=', ksplines, ', xt=list(m=', m, '))')
  } else {
    smoother <- paste0('s(sdimx, sdimy, bs="', basis, '", k=', ksplines, ')')
  }

  # ---- choose covariates ------------------------------------------------------
  reserved_std <- c("sdimx","sdimy","libsize")
  if (is.null(covariates)) {
    covnames <- character(0)
  } else if (identical(covariates, "all")) {
    # everything except the *original* loc/libsize columns and the standardized ones
    reserved_orig <- unique(c(loc_cols, libsize_col))
    covnames <- setdiff(colnames(metadata), c(reserved_std, reserved_orig))
  } else if (is.character(covariates)) {
    missing_cov <- setdiff(covariates, colnames(metadata))
    if (length(missing_cov)) warning("Dropping missing covariates: ", paste(missing_cov, collapse = ", "))
    covnames <- intersect(covariates, setdiff(colnames(metadata), reserved_std))
  } else {
    stop("covariates must be NULL, 'all', or a character vector.")
  }

  # ---- build model data & formula --------------------------------------------
  dat <- metadata[, unique(c("libsize","sdimx","sdimy", covnames)), drop = FALSE]

  if (!is.null(formula_override)) {
    # NOTE: formula_override must reference standardized names: libsize, sdimx, sdimy, and any covariates by name.
    form <- stats::as.formula(formula_override)
  } else {
    fprefix <- if (identical(family$family, "gaussian")) "y ~ " else "y ~ offset(log(libsize)) + "
    rhs <- paste(c(covnames, smoother), collapse = " + ")
    form <- stats::as.formula(paste0(fprefix, rhs))
  }

  # ---- labels & dims ----------------------------------------------------------
  if (is.null(colnames(assay_matrix))) colnames(assay_matrix) <- paste0("X", seq_len(ncol(assay_matrix)))
  targets   <- colnames(assay_matrix)
  num_cells <- nrow(assay_matrix)
  highdim   <- (num_cells >= 5000)

  # ---- run --------------------------------------------------------------------
  if (nCores > 1) {
    cl <- parallel::makeCluster(getOption("cl.cores", nCores))
    parallel::clusterEvalQ(cl, { library(mgcv) })
    pieces <- parallel::parLapply(
      cl, targets, getResiduals,
      data = dat, formula = form, assay_matrix = assay_matrix,
      num_cells = num_cells, family = family, highdim = highdim, verbose = verbose
    )
    parallel::stopCluster(cl)
  } else {
    pieces <- lapply(
      targets, getResiduals,
      data = dat, formula = form, assay_matrix = assay_matrix,
      num_cells = num_cells, family = family, highdim = highdim, verbose = verbose
    )
  }

  # ---- combine ---------------------------------------------------------------
  Residuals <- do.call(cbind, lapply(pieces, `[[`, "residuals"))
  colnames(Residuals) <- vapply(pieces, `[[`, character(1), "target")
  Residuals <- as.data.frame(Residuals)

  Splines <- do.call(cbind, lapply(pieces, function(x) as.numeric(x$splines[, 1])))
  colnames(Splines) <- vapply(pieces, `[[`, character(1), "target")
  Splines <- as.data.frame(Splines)

  Pvalues <- as.data.frame(t(vapply(pieces, function(x) x$pvalue, numeric(1))))
  names(Pvalues) <- vapply(pieces, `[[`, "", "target")
  rownames(Pvalues) <- "pvalue"

  list(Residuals = Residuals,
       Splines   = Splines,
       k.check.p = Pvalues,
       metadata  = dat)
}


