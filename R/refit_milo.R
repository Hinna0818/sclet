#' @title Refit Milo DA after dropping separated neighbourhoods
#'
#' @description 
#' This function removes neighbourhoods (nhoods) that show perfect separation 
#' (as detected by \code{checkSeparation}) from a Milo object and refits the 
#' differential abundance (DA) model. It supports both GLM and GLMM modes.
#'
#' @param milo A \code{Milo} object containing neighbourhood counts and metadata.
#' @param design_df Data.frame of experimental design information (samples x covariates).
#' @param bad_hood Integer vector of indices of nhoods to drop (e.g. from \code{checkSeparation}).
#' @param formula Model formula used in \code{testNhoods}.
#' @param contrasts Character string or vector of model contrasts.
#' @param glmm Logical; if TRUE fit a NB-GLMM, otherwise fit a NB-GLM. Default = TRUE.
#' @param glmm.solver Solver for GLMM variance components. Default = "Fisher".
#' @param REML Logical; use restricted maximum likelihood for GLMM. Default = TRUE.
#' @param max.iters Maximum number of iterations for GLMM fitting. Default = 50.
#' @param fail.on.error Logical; continue on errors if FALSE. Default = FALSE.
#' @param fdr.weighting Method for FDR weighting, passed to \code{testNhoods}. Default = "graph-overlap".
#' @param norm.method Normalisation method. Default = "TMM".
#' @param BPPARAM BiocParallel parameter. Default = \code{SerialParam()}.
#' @param old_results Optional; previous DA results for comparison of significant nhood counts.
#' @param seed Random seed for reproducibility. Default = 2025.
#' @param ... Additional arguments passed to \code{testNhoods}.
#'
#' @importFrom miloR testNhoods
#' @importFrom BiocParallel SerialParam
#' @return A list with:
#' \item{milo_clean}{Milo object with separated nhoods removed.}
#' \item{da_clean}{Data.frame of refitted DA results.}
#' \item{dropped}{Indices of dropped nhoods.}
#'
#' @export

refit_milo <- function(
    milo,
    design_df,
    bad_hood,
    formula,        
    contrasts,     
    glmm = TRUE,       
    glmm.solver = "Fisher",
    REML = TRUE,
    max.iters = 50,
    fail.on.error = FALSE,
    fdr.weighting = "graph-overlap",
    norm.method = "TMM",
    BPPARAM = BiocParallel::SerialParam(),
    old_results = NULL,
    seed = 2025,
    ...
) {
  set.seed(seed)
  
  # update new milo object
  milo2 <- if (length(bad_hood) > 0) milo[-bad_hood, ] else milo
  
  # refit
  if (glmm) {
    da_new <- miloR::testNhoods(
      milo2,
      design          = as.formula(formula),
      design.df       = design_df,
      model.contrasts = contrasts,
      fdr.weighting   = fdr.weighting,
      norm.method     = norm.method,
      glmm.solver     = glmm.solver,
      REML            = REML,
      max.iters       = max.iters,
      fail.on.error   = fail.on.error,
      BPPARAM         = BPPARAM,
      ...
    )
  } else {
    da_new <- miloR::testNhoods(
      milo2,
      design          = as.formula(formula),
      design.df       = design_df,
      model.contrasts = contrasts,
      fdr.weighting   = fdr.weighting,
      norm.method     = norm.method,
      BPPARAM         = BPPARAM,
      ...
    )
  }
  
  # comparison (optional)
  if (!is.null(old_results) && "SpatialFDR" %in% names(old_results)) {
    before <- sum(old_results$SpatialFDR < 0.1, na.rm = TRUE)
    after  <- sum(da_new$SpatialFDR < 0.1, na.rm = TRUE)
    message(sprintf("Significant nhoods (FDR<0.1): before=%s, after=%s", before, after))
  }
  
  return(list(milo_clean = milo2, da_clean = da_new, dropped = bad_hood))
}
