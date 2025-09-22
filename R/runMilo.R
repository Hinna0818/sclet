#' @title Run Milo Differential Abundance Analysis
#'
#' @description This function performs a complete Milo differential abundance (DA) analysis workflow on a
#' SingleCellExperiment (SCE) object. It constructs a k-nearest neighbor (KNN) graph, defines representative
#' neighborhoods, performs DA testing using edgeR, and returns both the modified SCE and Milo object.
#'
#' @param sce A \code{SingleCellExperiment} object containing single-cell data. Must have PCA and UMAP in \code{reducedDims}.
#' @param sample_col The column name in \code{colData(sce)} that indicates the sample label.
#' @param condition_col The column name in \code{colData(sce)} that indicates the experimental condition to test.
#' @param covariates optional character vector of additional columns as covariates for differential analysis in colData(sce)
#' @param contrasts character vector of contrast expressions
#' @param k1 Number of neighbors used in \code{buildGraph()}, default = 10.
#' @param k2 Number of neighbors used in \code{makeNhoods()}, default = 10.
#' @param d1 Number of PCA dimensions used in \code{buildGraph()}, default = 30.
#' @param d2 Number of PCA dimensions used in \code{makeNhoods()}, default = 30.
#' @param d3 Number of PCA dimensions used in \code{calcNhoodDistance()}, default = 30.
#' @param prop Proportion of cells to sample as neighborhood centers, default = 0.1.
#' @param per_contrast Logical; test contrasts separately (\code{TRUE}) or jointly (\code{FALSE}).
#' @param is_refined Logical; whether to apply refinement algorithm when sampling neighborhoods, default = TRUE.
#' @param refinement_scheme Refinement scheme passed to \code{makeNhoods}.
#' @param adjustment Method for multiple testing correction within neighborhoods (default \code{"BH"}).
#' @param seed Random seed for reproducibility, default = 2025.
#'
#' @return A list with:
#' \item{sce}{Input \code{SingleCellExperiment} with DA results stored in metadata.}
#' \item{milo}{Milo object with graph, neighborhoods, and DA results.}
#' \item{da_results}{Data.frame of DA testing results.}
#'
#' @importFrom SingleCellExperiment reducedDim reducedDim<- reducedDimNames colData metadata<-
#' @importFrom miloR Milo, buildGraph, makeNhoods, countCells, calcNhoodDistance, testNhoods, buildNhoodGraph
#' @importFrom dplyr distinct
#' @export
#' 
runMilo <- function(
    sce,
    sample_col,
    condition_col,
    covariates = NULL,
    contrasts = NULL,
    k1 = 10,
    k2 = 10,
    d1 = 30,
    d2 = 30,
    d3 = 30,
    prop = 0.1,
    per_contrast = TRUE, 
    is_refined = TRUE,
    refinement_scheme = "graph",
    adjustment = "BH",
    seed = 2025
) {
  
  set.seed(seed)
  
  if (!("PCA"  %in% reducedDimNames(sce)))  stop("PCA not found in reducedDims(sce).")
  if (!("UMAP" %in% reducedDimNames(sce)))  stop("UMAP not found in reducedDims(sce).")
  if (!(sample_col %in% colnames(colData(sce)))) stop("`sample_col` not found in colData(sce).")
  if (!(condition_col %in% colnames(colData(sce)))) stop("`condition_col` not found in colData(sce).")
  
  if (is.null(contrasts)) {
    stop("A valid contrasts vector should be supplied, e.g. 'ConditionB - ConditionA'.")
  }
  
  ## create miloR object
  milo <- miloR::Milo(sce)
  reducedDim(milo, "UMAP") <- reducedDim(sce, "UMAP")
  
  ## build KNN graph and find neighborhoods for each cell
  milo <- miloR::buildGraph(milo, k = k1, d = d1)
  milo <- miloR::makeNhoods(milo, prop = prop, k = k2, d = d2, refined = is_refined, refinement_scheme = refinement_scheme)
  
  ## count cells by sample
  milo <- miloR::countCells(milo, meta.data = data.frame(colData(milo)), samples = sample_col)
  
  ## create experiment-design dataframe
  need_cols <- unique(c(sample_col, condition_col, covariates))
  design_df <- data.frame(colData(milo))[, need_cols, drop = FALSE]
  design_df <- dplyr::distinct(design_df)
  rownames(design_df) <- design_df[[sample_col]]
  design_df <- design_df[colnames(nhoodCounts(milo)), , drop = FALSE]
  
  ## enforce factor and drop unused levels
  design_df[[condition_col]] <- droplevels(factor(design_df[[condition_col]]))
  
  ## create design matrix(~ 0 + condition + covariates)
  rhs <- c(paste0("0 + ", condition_col), covariates)
  fml <- as.formula(paste("~", paste(rhs, collapse = " + ")))
  
  ## run differential analysis
  if (!per_contrast) {
    ## combine QL F test
    da_results <- miloR::testNhoods(
      milo,
      design = fml,
      design.df = design_df,
      model.contrasts = contrasts,
      fdr.weighting = "graph-overlap",
      norm.method = "TMM"
    )
  } else {
    ## compare group by group
    res_list <- lapply(contrasts, function(ct) {
      r <- miloR::testNhoods(
        milo,
        design = fml,
        design.df = design_df,
        model.contrasts = ct,
        fdr.weighting = "graph-overlap",
        norm.method = "TMM"
      )
      r$contrast <- ct
      r
    })
    
    da_results <- do.call(rbind, res_list)
    
    ## BH adjustment(optional)
    da_results$SpatialFDR_adj_withinNhood <- ave(
      da_results$SpatialFDR,
      da_results$Nhood,
      FUN = function(x) p.adjust(x, method = adjustment)
    )
  }
  
  ## build cell neighborhood graph
  milo <- miloR::calcNhoodDistance(milo, d = d3)
  milo <- miloR::buildNhoodGraph(milo)
  
  ## store results
  metadata(sce)$milo_results <- list(
    da_results = da_results,
    design_df  = design_df,
    formula    = deparse(fml),
    contrasts  = contrasts,
    params     = list(k1 = k1, k2 = k2, d1 = d1, d2 = d2, d3 = d3,
                      prop = prop, is_refined = is_refined, seed = seed))
  
  return(list(
    sce = sce,
    milo = milo,
    da_results = da_results
  ))
  
}
