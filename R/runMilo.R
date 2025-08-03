#' @title Run Milo Differential Abundance Analysis
#'
#' @description This function performs a complete Milo differential abundance (DA) analysis workflow on a
#' SingleCellExperiment (SCE) object. It constructs a k-nearest neighbor (KNN) graph, defines representative
#' neighborhoods, performs DA testing using edgeR, and returns both the modified SCE and Milo object.
#'
#' @param sce A \code{SingleCellExperiment} object containing single-cell data. Must have PCA and UMAP in \code{reducedDims}.
#' @param sample_col The column name in \code{colData(sce)} that indicates the sample label.
#' @param condition_col The column name in \code{colData(sce)} that indicates the experimental condition to test.
#' @param k1 Number of neighbors used in \code{buildGraph()}, default = 10.
#' @param k2 Number of neighbors used in \code{makeNhoods()}, default = 10.
#' @param d1 Number of PCA dimensions used in \code{buildGraph()}, default = 30.
#' @param d2 Number of PCA dimensions used in \code{makeNhoods()}, default = 30.
#' @param d3 Number of PCA dimensions used in \code{calcNhoodDistance()}, default = 30.
#' @param prop Proportion of cells to sample as neighborhood centers, default = 0.1.
#' @param is_refined Logical; whether to apply refinement algorithm when sampling neighborhoods, default = TRUE.
#' @param seed Random seed for reproducibility, default = 2025.
#'
#' @return A list with the following components:
#' \item{sce}{The input \code{SingleCellExperiment} object with DA results added to \code{metadata()}.}
#' \item{milo}{The resulting \code{Milo} object containing graph, neighborhood, and DA information.}
#' \item{da_results}{A data.frame of differential abundance testing results for each neighborhood.}
#'
#' @importFrom SingleCellExperiment reducedDim reducedDim<- reducedDimNames colData metadata<-
#' @importFrom miloR Milo, buildGraph, makeNhoods, countCells, calcNhoodDistance, testNhoods, buildNhoodGraph
#' @importFrom dplyr distinct
#' @importFrom S4Vectors metadata
#' @export
#' 
runMilo <- function(
    sce,
    sample_col,
    condition_col,
    k1 = 10,
    k2 = 10,
    d1 = 30,
    d2 = 30,
    d3 = 30,
    prop = 0.1,
    is_refined = TRUE,
    seed = 2025
) {
  set.seed(seed)
  
  stopifnot("PCA" %in% reducedDimNames(sce))
  message("Please run PCA before runMilo.")
  
  stopifnot("UMAP" %in% reducedDimNames(sce))
  message("Please run UMAP before runMilo.")
  
  stopifnot(sample_col %in% colnames(colData(sce)))
  message("Please add the sample information into the colData of sce.")
  
  stopifnot(condition_col %in% colnames(colData(sce)))
  message("Please add the condition information into the colData of sce.")
  
  milo <- miloR::Milo(sce)
  reducedDim(milo, "UMAP") <- reducedDim(sce, "UMAP")
  
  milo <- miloR::buildGraph(milo, k = k1, d = d1)
  
  milo <- miloR::makeNhoods(milo, prop = prop, k = k2, d = d2, refined = is_refined)
  
  milo <- miloR::countCells(milo, meta.data = data.frame(colData(milo)), samples = sample_col)
  
  design_df <- data.frame(colData(milo))[, c(sample_col, condition_col)]
  design_df <- dplyr::distinct(design_df)
  rownames(design_df) <- design_df[[sample_col]]
  design_df <- design_df[colnames(nhoodCounts(milo)), , drop = FALSE]
  
  milo <- miloR::calcNhoodDistance(milo, d = d3)
  rownames(design_df) <- design_df[[sample_col]]
  da_results <- miloR::testNhoods(
    milo, 
    design = as.formula(paste0("~ ", condition_col)), 
    design.df = design_df)
  
  milo <- miloR::buildNhoodGraph(milo)
  
  S4Vectors::metadata(sce)$milo_DA_results <- da_results
  
  return(list(
    sce = sce,
    milo = milo,
    da_results = da_results
  ))
}