# Example dataset preprocessing: KotliarovPBMCData (RNA)
#
# Reference:
# Kotliarov et al. (2020). 
# Broad immune activation underlies shared set point signatures 
# for vaccine responsiveness in healthy individuals and disease 
# activity in patients with lupus. Nat. Med. 26, 618–629
#
# Description:
# - The KotliarovPBMCData dataset contains PBMC (peripheral blood 
#   mononuclear cells) single-cell RNA-seq data from healthy donors 
#   and patients with lupus. 
# - We use this dataset as a reproducible example for testing MiloR 
#   differential abundance (DA) analysis workflows.
#

library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(scRNAseq)
library(scuttle)
library(irlba)
library(BiocParallel)
library(ggplot2)
library(sparseMatrixStats)

## data preprocessing
# download raw data
bpparam <- SerialParam()
register(bpparam)

pbmc.sce <- KotliarovPBMCData(mode="rna", ensembl=TRUE) 

# downsample cells to reduce compute time
prop.cells <- 0.75
n.cells <- floor(ncol(pbmc.sce) * prop.cells)
set.seed(42)
keep.cells <- sample(colnames(pbmc.sce), size=n.cells)
pbmc.sce <- pbmc.sce[, colnames(pbmc.sce) %in% keep.cells]

# downsample the number of samples
colData(pbmc.sce)$ObsID <- paste(colData(pbmc.sce)$tenx_lane, colData(pbmc.sce)$sampleid, sep="_")
n.samps <- 80
keep.samps <- sample(unique(colData(pbmc.sce)$ObsID), size=n.samps)
keep.cells <- rownames(colData(pbmc.sce))[colData(pbmc.sce)$ObsID %in% keep.samps]
pbmc.sce <- pbmc.sce[, colnames(pbmc.sce) %in% keep.cells]

pbmc.sce

# data preprocessing and normalization
set.seed(42)
# remove sparser cells
drop.cells <- colSums(counts(pbmc.sce)) < 1000
pbmc.sce <- computePooledFactors(pbmc.sce[, !drop.cells], BPPARAM=bpparam)
pbmc.sce <- logNormCounts(pbmc.sce)

pbmc.hvg <- modelGeneVar(pbmc.sce)
pbmc.hvg$FDR[is.na(pbmc.hvg$FDR)] <- 1

pbmc.sce <- runPCA(pbmc.sce, subset_row=rownames(pbmc.sce)[pbmc.hvg$FDR < 0.1], scale=TRUE, ncomponents=50, assay.type="logcounts", name="PCA", BPPARAM=bpparam)
pbmc.sce

set.seed(42)
pbmc.sce <- runUMAP(pbmc.sce, n_neighbors=30, pca=50, BPPARAM=bpparam)

colData(pbmc.sce)$ObsID <- paste(colData(pbmc.sce)$tenx_lane, colData(pbmc.sce)$sampleid, sep="_")

## save
save(pbmc.sce, file = "./pbmc.sce.rda")

