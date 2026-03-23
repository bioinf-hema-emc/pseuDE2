
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pseuDE2

<!-- badges: start -->

<!-- badges: end -->

Perform pseudobulk differential expression analysis for single cell RNA-seq data using Seurat.

## Installation

You can install pseuDE2 using the following devtools command:

``` r
devtools::install_github("bioinf-hema-emc/pseuDE2")
```

## Introduction

Normal differential expression analysis (e.g. using Seurat::FindMarkers)
treats all cells as independent measurements. Strictly speaking this is
not fair as individual cells from the same sample are not independent.
This also explains the inflated p-values one typically obtains when
performing differential expression analysis with single cell data. A
common method to make the measurements independent is by summing the
counts of individual cells per sample per group (e.g. cell type). This
creates a single independent ‘measurement’ where the measurement is now
representing a bulk RNA experiment per sample per group (called
pseudobulk). Now we can compare the independent measurements between
conditions (e.g. before and after treatment). For a proper statistical
test this requires at least 5 independent samples per condition (thus at
least 5 single cell experiments per condition). Pseudobulk also allows
for performing a paired differential analysis. This can be done when the
same samples are sequenced between conditions. The between sample
variability can then be accounted for when performing the differential
expression analysis between conditions.

Design equations:

- Unpaired analysis: ‘~ + com_group’
- Paired analysis: ‘~ + pair_group + com_group’

## Example

We have a Seurat object with cells from two conditions (A and B). For 5
patients single cell sequencing is performed in both conditions. We want
to compare condition A with condition B in T-cells using 5 patients.
There are two approaches we can take.

Approach 1:

We create a metadata entry where we concatenate the metadata that
defines the conditions and the cell types. We would then get, for each
cell, something like cellType_condition (e.g. Tcell_conditionA;
Tcell_conditionB, neutrophil_conditionA, etc.). This allows us to tell
the function to compare conditionA with conditionB only for the T-cells.
This new metadata entry is going to be our choice for
`compare_groups_metadata` together with
`comp_group1 = 'Tcell_conditionA'` and
`comp_group2 = 'Tcell_conditionB'`.

Approach 2:

We subset our Seurat object to only include T-cells, for example with
`subset(object, subset = cellTypes == 'Tcells')`. Now we can make the
comparison between conditionA and conditionB directly without adding a
new metadata entry. The entry for ‘compare_groups_metadata’ is now
‘Condition’ with comp_group1 = ‘conditionA’ and comp_group2 =
‘conditionB’.

For either approach, the metadata that defines the patients is going to
be our entry for ‘aggregate_groups_metadata’. Since the same patients
are present in both condition A and condition B, when can perform a
paired comparison by setting `paired = T`. This corrects for any
biological variability within the patients between the conditions. When
`paired = TRUE`, the function creates a table with two columns called
‘com_group’ and ‘pair_group’:

|                       | com_group        | pair_group |
|-----------------------|------------------|------------|
| pat1-Tcell_conditionA | Tcell_conditionA | pat1       |
| pat1-Tcell_conditionB | Tcell_conditionB | pat1       |
| pat2-Tcell_conditionA | Tcell_conditionA | pat2       |
| pat2-Tcell_conditionB | Tcell_conditionB | pat2       |
| pat3-Tcell_conditionA | Tcell_conditionA | pat3       |
| pat3-Tcell_conditionB | Tcell_conditionB | pat3       |
| pat4-Tcell_conditionA | Tcell_conditionA | pat4       |
| pat4-Tcell_conditionB | Tcell_conditionB | pat4       |
| pat5-Tcell_conditionA | Tcell_conditionA | pat5       |
| pat5-Tcell_conditionB | Tcell_conditionB | pat5       |

In case `paired = FALSE`, the ‘pair_group’ column is not created and the
comparison is only made based on the ‘com_group’ column.

``` r
library(pseuDE2)
library(Seurat)
library(EnhancedVolcano)

pseude_res <- pseuDE2(seurat_object,
aggregate_groups_metadata = 'Patients',
compare_group_metadata = 'cellType_condition',
outputDir = 'output/',
comp_group1 = 'Tcell_conditionA',
comp_group2 = 'Tcell_conditionB',
paired = T,
cores = 3)

res_counts <- pseude_res[[1]] # Table with log2FC, p-values and normalized counts.
de_out <- pseude_res[[2]] # DESeq2 object

# CREATE VOLCANOPLOT ----
EnhancedVolcano(res_counts,
                rownames(res_counts),
                x ="log2FoldChange",
                y ="padj",
                drawConnectors = TRUE,
                cutoffLineWidth = 0.3,
                widthConnectors = 0.2,
                max.overlaps = 60,
                xlab = bquote(~Log[2]~ "fold change"),
                title = paste0("Pseudobulk DEG"),
                pCutoff = 0.05,
                FCcutoff = 0.499,
                gridlines.major = FALSE,
                gridlines.minor = FALSE
)

# CREATE FEATUREPLOTS ----
DimPlot(seurat_object,
        reduction = 'umap',
        group.by = 'Patients',
        raster = F) +
    DimPlot(seurat_object,
            reduction = 'umap',
            group.by = 'cellType_condition',
            raster = F) +
    FeaturePlot(seurat_object,
                features = head(rownames(res_counts), 5),
                raster = F,
                order = T)
```
