# Example dataset preprocessing: MouseSMARTseqData (Thymus)
#
# Reference:
# MouseThymusAgeing package (Single-cell RNA-seq data resource)
# The dataset consists of SMART-seq2 profiles of mouse thymic 
# cells sampled across different ages and sorting days, designed 
# to study immune ageing.
#
# Description:
# - The dataset captures thymic cell populations across age groups, 
#   allowing exploration of age-related changes in cell abundance 
#   and transcriptional state.
# - We use this as a smaller example dataset to demonstrate MiloR 
#   DA analysis on mouse thymus single-cell data.

library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(MouseThymusAgeing)
library(scuttle)

## data preprocessing
# load data
thy.sce <- MouseSMARTseqData()
colData(thy.sce)$Sample <- paste(colData(thy.sce)$SortDay, colData(thy.sce)$Age, sep="_")

# Normalization and reduction
thy.sce <- logNormCounts(thy.sce)
thy.sce <- runUMAP(thy.sce)

# save
save(thy.sce, file = "./thy.sce.rda")


