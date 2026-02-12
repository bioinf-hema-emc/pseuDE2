################################################################################
### Author:             Gregory van Beek
### Date:               02-01-2026
### Email:              g.vanbeek@erasmusmc.nl
###
### Version:            2.1
###
### Description:        Pseudobulk analysis for single cell analysis using DEseq2.
###                     Both unpaired as well as paired sample analysis.
###                     Based on the algorithm from PseuDE [https://github.com/weversMJW/pseude]
###                     It loads the following function(s):
###                         pseuDE2()
################################################################################

#' @title pseuDE2
#'
#' @description
#' Perform differential expression analysis using pseudobulk.
#' It requires a Seurat object as main input.
#' The different metadata entries in the Seurat object can be used to define how the data needs to be aggregated and which comparison should be made.
#' A positive log2FC means an up-regulation in comp_group1 relative to comp_group2 (or all other cells when comp_group2 is not defined).
#'
#' @param object Seurat object. Different layers in the RNA assay will be merged. Only the raw RNA-counts are used in this function.
#' @param aggregate_groups_metadata Metadata entry in object. Aggregation will be based on grouping in this metadata.
#' @param compare_groups_metadata Metadata entry in object. Should contain at least two different groups for which the DE analysis is performed.
#' @param outputDir Path to a folder where results are stored. It will be created if it does not exist.
#' @param comp_group1 Value present in compare_groups_metadata. This will be the reference group in the comparison.
#' @param comp_group2 Default = NULL. Value present in compare_groups_metadata. This will be the test group in the comparison. If not set, it automatically takes all cells not in comp_group1. If there are more than 2 groups and this is not set, all other cells will be assigned to a group 'x'.
#' @param do_prefilter Default = TRUE. Either TRUE (T) or FALSE (F). When set to TRUE, then only the genes are kept that have more than 2 fragments per million (fpm) in more than half of the samples in either of the two groups.
#' @param paired Default = FALSE. Either TRUE (T) or FALSE (F). Whether to perform pseudobulk DE analysis using paired samples.
#' @param pct Default = 0.1. Value between 0.0 and 1.0. Ratio of cells in which a gene should be present in either of the two groups in order to be taken into account in the DE analysis.
#' @param order_results Default = padj. Should be one of the following strings: 'padj', 'log2FoldChange', 'pvalue', 'baseMean' or 'lfcSE'. Defines how the final results should be ordered. It does not influence the results, only the ordering with which it is saved.
#' @param saveRDS Default = FALSE. Either TRUE (T) or FALSE (F). Whether to save the DESeq2 object in outputDir.
#' @param return_results Default = TRUE. Either TRUE (T) or FALSE (F). Whether to return the results (when TRUE) or only save the results in outputDir (when FALSE).
#' @param saveName Default = NULL. String. Name to add to the filenames. If NULL, then the output files will be named 'log2fc.csv' and, if saveRDS == TRUE, 'DESeq2object.RDS'. For example, if saveName = 'sample1', then the filenames will become 'log2fc_sample1.csv' and 'DESeq2object_sample1.RDS'.
#' @param cores Default = 1. Positive integer value. How many cores to be used when running DESeq2. When value is > 1, then BiocParallel is used for parallelization.
#'
#' @return When return_results = TRUE: List containing two objects:
#' \itemize{
#' 	\item A matrix with the DESeq2 lfcShrink output and the normalized counts
#' 	\item The DESeq2 output object.
#' }
#'
#' @details
#' Normal differential expression analysis (e.g. using Seurat::FindMarkers) treats all cells as independent measurements.
#' Strictly speaking this is not fair as individual cells from the same sample are not independent.
#' This also explains the inflated p-values one typically obtains when performing differential expression analysis with single cell data.
#' A common method to make the measurements independent is by summing the counts of individual cells per sample per group (e.g. cell type).
#' This creates a single independent 'measurement' where the measurement is now representing a bulk RNA experiment per sample per group (called pseudobulk).
#' Now we can compare the independent measurements between conditions (e.g. before and after treatment).
#' For a proper statistical test this requires at least 5 independent samples per condition (thus at least 5 single cell experiments per condition).
#' Pseudobulk also allows for performing a paired differential analysis.
#' This can be done when the same samples are sequenced between conditions.
#' The between sample variability can then be accounted for when performing the differential expression analysis between conditions.
#'
#' Design equations:
#' \itemize{
#' \item Unpaired analysis: '~ + com_group'
#' \item Paired analysis: '~ + pair_group + com_group'
#' }
#'
#'
#' EXAMPLE
#'
#' We have a Seurat object with cells from two conditions (A and B).
#' For 5 patients single cell sequencing is performed in both conditions.
#' We want to compare condition A with condition B in T-cells using 5 patients.
#' There are two approaches we can take.
#'
#' Approach 1:
#'
#' We create a metadata entry where we concatenate the metadata that defines the conditions and the cell types.
#' We would then get, for each cell, something like cellType_condition (e.g. Tcell_conditionA; Tcell_conditionB, neutrophil_conditionA, etc.).
#' This allows us to tell the function to compare conditionA with conditionB only for the T-cells.
#' This new metadata entry is going to be our choice for 'compare_groups_metadata' together with comp_group1 = 'Tcell_conditionA' and comp_group2 = 'Tcell_conditionB'.
#'
#' Approach 2:
#'
#' We subset our Seurat object to only include T-cells, for example with `subset(object, subset = cellTypes == 'Tcells')`.
#' Now we can make the comparison between conditionA and conditionB directly without adding a new metadata entry.
#' The entry for 'compare_groups_metadata' is now 'Condition' with comp_group1 = 'conditionA' and comp_group2 = 'conditionB'.
#'
#' For either approach, the metadata that defines the patients is going to be our entry for 'aggregate_groups_metadata'.
#' Since the same patients are present in both conditionA and conditionB, when can perform a paired comparison by setting paired = T. This corrects for any biological variability wihin the patients between the conditions.
#' When paired = TRUE, the function creates a table with two columns called 'com_group' and 'pair_group':
#'
#' \tabular{llrr}{
#'    \tab \tab \strong{com_group} \tab \strong{pair_group} \cr
#'   pat1-Tcell_conditionA \tab \tab Tcell_conditionA \tab pat1 \cr
#'   pat1-Tcell_conditionB \tab \tab Tcell_conditionB \tab pat1 \cr
#'   pat2-Tcell_conditionA \tab \tab Tcell_conditionA \tab pat2 \cr
#'   pat2-Tcell_conditionB \tab \tab Tcell_conditionB \tab pat2 \cr
#'   pat3-Tcell_conditionA \tab \tab Tcell_conditionA \tab pat3 \cr
#'   pat3-Tcell_conditionB \tab \tab Tcell_conditionB \tab pat3 \cr
#'   pat4-Tcell_conditionA \tab \tab Tcell_conditionA \tab pat4 \cr
#'   pat4-Tcell_conditionB \tab \tab Tcell_conditionB \tab pat4 \cr
#'   pat5-Tcell_conditionA \tab \tab Tcell_conditionA \tab pat5 \cr
#'   pat5-Tcell_conditionB \tab \tab Tcell_conditionB \tab pat5 \cr
#' }
#'
#' In case paired = FALSE, the 'pair_group' column is not created and the comparison is only made based on the 'com_group' column.
#'
#' @examples
#' \dontrun{
#'
#' library(pseuDE2)
#' library(Seurat)
#' library(EnhancedVolcano)
#'
#' pseude_res <- pseuDE2(sr,
#' aggregate_groups_metadata = 'Patients',
#' compare_group_metadata = 'cellType_condition',
#' outputDir = 'output/',
#' comp_group1 = 'Tcell_conditionA',
#' comp_group2 = 'Tcell_conditionB',
#' paired = T,
#' cores = 3)
#'
#' res_counts <- pseude_res[[1]] # Table with log2FC, p-values and normalized counts.
#' de_out <- pseude_res[[2]] # DESeq2 object
#' 
#' # CREATE VOLCANOPLOT ----
#' EnhancedVolcano(res_counts,
#'                 rownames(res_counts),
#'                 x ="Log2FoldChange",
#'                 y ="padj",
#'                 drawConnectors = TRUE,
#'                 cutoffLineWidth = 0.3,
#'                 widthConnectors = 0.2,
#'                 max.overlaps = 60,
#'                 xlab = bquote(~Log[2]~ "fold change"),
#'                 title = paste0("Pseudobulk DEG"),
#'                 pCutoff = 0.05,
#'                 FCcutoff = 0.499,
#'                 gridlines.major = FALSE,
#'                 gridlines.minor = FALSE
#' )
#'
#' DimPlot(sr,
#'         reduction = 'umap',
#'         group.by = 'Patients',
#'         raster = F) +
#'     DimPlot(sr,
#'             reduction = 'umap',
#'             group.by = 'cellType_condition',
#'             raster = F) +
#'     FeaturePlot(sr,
#'                 features = head(rownames(deg_table), 5),
#'                 raster = F,
#'                 order = T)
#'
#' }
#'
#' @export
pseuDE2 <- function(object,
                    aggregate_groups_metadata,
                    compare_groups_metadata,
                    outputDir,
                    comp_group1,
                    comp_group2 = NULL,
                    do_prefilter = T,
                    paired = F,
                    pct = 0.1,
                    order_results = 'padj',
                    saverds = F,
                    return_results = T,
                    saveName = NULL,
                    cores = 1){
    
    # INITIAL CHECKS
    if(! isS4(object) ) {
        stop('Object should be a Seurat object.')
    }
    
    if(! is.character(aggregate_groups_metadata) | ! aggregate_groups_metadata %in% colnames(object@meta.data) ) {
        stop('aggregate_groups_metadata should be the name of a metadata entry in object')
    }
    
    if(! is.character(compare_groups_metadata) | ! compare_groups_metadata %in% colnames(object@meta.data) ) {
        stop('compare_groups_metadata should be the name of a metadata entry in object')
    }
    
    if( ! dir.exists(outputDir) ) {
        writeLines(paste0('WARNING: outputDir does not exist.\nCreating following directory: ', outputDir))
        dir.create(outputDir)
    }
    
    if(! comp_group1 %in% object@meta.data[[compare_groups_metadata]] ) {
        stop('comp_group1 not found in compare_groups_metadata. Please enter a value that is present in compare_groups_metadata.')
    }
    
    if(! is.null(comp_group2) ) {
        if(! comp_group2 %in% object@meta.data[[compare_groups_metadata]] ) {
            stop('comp_group2 not found in compare_groups_metadata. Please enter a value that is present in compare_groups_metadata.')
        }
    }
    
    if(! is.logical(do_prefilter) ) {
        stop('do_prefilter should be TRUE (T) or FALSE (F).')
    }
    
    if(! is.logical(paired) ) {
        stop('paired should be TRUE (T) or FALSE (F).')
    }
    
    if ( pct > 1.0 | pct < 0.0 | ! is.numeric(pct)) {
        stop('pct requires a numeric value between 0.0 and 1.0 or NULL.')
    }
    
    if (! order_results %in% c("padj", "log2FoldChange", "pvalue", "baseMean", "lfcSE") ) {
        writeLines("WARNING: order_results should one of 'padj', 'log2FoldChange', 'pvalue', 'baseMean' or 'lfcSE'.\nResetting order_results = 'padj'")
        order_results <- 'padj'
    }
    
    if(! is.logical(saverds) ) {
        stop('saverds should be TRUE (T) or FALSE (F).')
    }
    
    if(! is.logical(return_results) ) {
        stop('return_results should be TRUE (T) or FALSE (F)')
    }
    
    if(! is.character(saveName) & ! is.null(saveName) ) {
        stop('saveName should be NULL or a character string.')
    }
    
    if(! is.numeric(cores) | cores < 1) {
        stop('cores should be a positive integer value')
    }
    
    if( length(grep('-', object@meta.data[[compare_groups_metadata]])) != 0 ){
        object@meta.data[[compare_groups_metadata]] <- stringr::str_replace(object@meta.data[[compare_groups_metadata]], '-', '_')
        comp_group1 <- stringr::str_replace(comp_group1, '-', '_')
        comp_group2 <- stringr::str_replace(comp_group2, '-', '_')
    }
    
    if( length(grep('-', object@meta.data[[aggregate_groups_metadata]])) != 0 ){
        object@meta.data[[aggregate_groups_metadata]] <- stringr::str_replace(object@meta.data[[aggregate_groups_metadata]], '-', '_')
    }
    
    
    # SETUP ENVIRONMENT
    writeLines('Loading Seurat and DESeq2 packages.')
    require(Seurat)
    require(DESeq2)
    
    
    if( cores != 1 ) {
        cores <- BiocParallel::register(BiocParallel::MulticoreParam(cores))
    }
    
    
    # GET COUNTS MATRIX
    writeLines('Joining layers RNA-assay.')
    DefaultAssay(object) <- 'RNA'
    object <- JoinLayers(object)
    
    
    # PREPARE OBJECT
    object@meta.data[[compare_groups_metadata]] <- as.character(object@meta.data[[compare_groups_metadata]])
    
    if ( length(unique(object@meta.data[[compare_groups_metadata]])) > 2 ) { # If there are more than two different values in compare_groups_metadata
        if (! is.null(comp_group2) | ! all(is.na(comp_group2)) ) {
            cells_use <- rownames(object@meta.data)[object@meta.data[[compare_groups_metadata]] == comp_group1 | object@meta.data[[compare_groups_metadata]]  == comp_group2]
            object <- subset(object, cells = cells_use) # Subset object to include only comp_group1 and comp_group2
        } else {
            object@meta.data[['pseudobulk']] <- object@meta.data[[compare_groups_metadata]]
            object@meta.data[['pseudobulk']][object@meta.data[['pseudobulk']] != comp_group1] <- 'x' # Set cell not in comp_group1 to 'x'
            compare_groups_metadata <- 'pseudobulk'
            comp_group2 <- 'x'
        }
    } else if ( length(unique(object@meta.data[[compare_groups_metadata]])) == 2 ) {
        if ( is.null(comp_group2)) # If value comp_group2 is not set and there are only two values in compare_groups_metadata, then set comp_group2 to the value other than comp_group1
            comp_group2 <- object@meta.data[[compare_groups_metadata]][object@meta.data[[compare_groups_metadata]] != comp_group1][1]
    } else {
        stop(paste0('There are too few values in ', compare_groups_metadata, '. Ensure there are at least two values to compare.'))
    }
    
    
    # GET COUNTS MATRIX
    sc_counts <- object@assays[['RNA']]@layers[['counts']]
    
    if (class(sc_counts)[1] != "dgCMatrix") {
        sc_counts <- Matrix::Matrix(as.matrix(sc_counts), sparse = TRUE)
    }
    
    colnames(sc_counts) <- colnames(object)
    rownames(sc_counts) <- rownames(object)
    
    
    # DEFINE GROUPS TO COMPARE
    compare_groups <- unique(object@meta.data[[compare_groups_metadata]]) # Should be exactly two names
    compare_groups <- compare_groups[! is.na(compare_groups)]
    
    
    # CREATE AGGREGATION DATAFRAME
    agg_df <- data.frame('agg_comp' = paste0(object@meta.data[[aggregate_groups_metadata]], '-', object@meta.data[[compare_groups_metadata]]),
                         'com_group' = object@meta.data[[compare_groups_metadata]],
                         'agg_group' = object@meta.data[[aggregate_groups_metadata]],
                         row.names = rownames(object@meta.data))
    
    
    # REMOVE GENE IT IT HAS EXPRESSION IN LESS THAN PCT PERCENTAGE OF CELLS IN EITHER COMPARE_GROUPS
    if(! is.null(pct) & pct != 0.0 ) {
        writeLines(paste0('Removing genes that are expressed in less than ', pct * 100, ' percent of cells in either group.'))
        pct_genes_group1 <- Matrix::rowSums(sc_counts[, rownames(agg_df)[agg_df$com_group == comp_group1]] > 0) / dim(sc_counts[, rownames(agg_df)[agg_df$com_group == comp_group1]])[2]
        pct_genes_group2 <- Matrix::rowSums(sc_counts[, rownames(agg_df)[agg_df$com_group == comp_group2]] > 0) / dim(sc_counts[, rownames(agg_df)[agg_df$com_group == comp_group2]])[2]
        pct_genes_keep <- unique(c(names(pct_genes_group1[pct_genes_group1 > pct]), names(pct_genes_group2[pct_genes_group2 >= pct]))) #Keep genes that are expressed in more than pct of the cells.
        
        writeLines(paste0(length(pct_genes_keep), ' of ',dim(sc_counts)[1] , ' genes are kept.'))
        
        sc_counts <- sc_counts[pct_genes_keep,]
    }
    
    
    # CREATE AGGREGATE MATRIX
    writeLines('Aggregating counts matrix.')
    agg_form <- stats::as.formula('~ 0 + agg_comp') # Intersect is not important for aggregation. Removed with '0'.
    
    agg_map <- Matrix::sparse.model.matrix(agg_form, agg_df, row.names = FALSE)
    colnames(agg_map) <- sub("^agg_comp", "", colnames(agg_map))
    
    agg_counts = sc_counts %*% agg_map
    rownames(agg_counts) <- rownames(sc_counts)
    
    writeLines('Rounding counts values to make it compatible with DESeq2.')
    agg_counts <- round(agg_counts)
    
    
    # CREATE COLDATA ENTRY FOR DESEQ2
    if(paired){
        pair_group <- unlist(lapply(stringr::str_split(colnames(agg_counts), '-'), '[[', 1))
        com_group <- unlist(lapply(stringr::str_split(colnames(agg_counts), '-'), '[[', 2))
        agg_table <- data.frame('com_group' = com_group,
                                'pair_group' = pair_group,
                                row.names = colnames(agg_counts))
        
        Ncells <- aggregate(rep(1, nrow(agg_df)), by = list(com_group = agg_df$com_group, agg_group = agg_df$agg_group), sum)
        
        design <- stats::as.formula('~ + pair_group + com_group')
    } else {
        com_group <- unlist(lapply(stringr::str_split(colnames(agg_counts), '-'), '[[', 2))
        agg_table <- data.frame('com_group' = com_group,
                                row.names = colnames(agg_counts))
        
        Ncells <- data.frame(table(agg_df$com_group))
        colnames(Ncells) <- c('com_group', 'x')
        
        design <- stats::as.formula('~ + com_group')
    }
    
    print(Ncells)    
    if( sum(Ncells$x < 100) > 0 ) {
        writeLines('WARNING: Some groups have < 100 cells. Check reliability of results!')
    }
    
    
    # PERFORM DIFFERENTIAL EXPRESSION ANALYSIS
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = agg_counts,
                                          colData = agg_table,
                                          design = design)
    
    if (do_prefilter) { # More than half of the samples should have at least two fpm counts in at least one of the two groups.
        writeLines('Performing DESeq2 analysis with prefiltering.')
        dds_frag <- DESeq2::fpm(dds)
        
        con_id1 <- grep(paste0('-', comp_group1), colnames(dds_frag))
        con_id2 <- grep(paste0('-', comp_group2), colnames(dds_frag))
        kp_gene <- rowSums(dds_frag[, con_id1] >= 2) >= round(length(con_id1) / 2) |
            rowSums(dds_frag[, con_id2] >= 2) >= round(length(con_id2) / 2)
        writeLines(paste0(table(kp_gene)[[2]], ' of ',dim(sc_counts)[1] , ' genes (', round(table(kp_gene)[[2]] / dim(sc_counts)[1] *100, digits = 2), ' percent) are kept after prefiltering.'))
        
        de_out <- DESeq2::DESeq(dds[kp_gene, ], parallel = TRUE)
    } else {
        writeLines('Performing DESeq2 analysis without prefiltering.')
        de_out <- DESeq2::DESeq(dds, parallel = TRUE)
    }
    
    
    # GENERATE RESULTS WITH ASHR-SHRUNK LOG2FC
    writeLines('Preparing output file.')
    res <- data.frame(DESeq2::lfcShrink(de_out,
                                        contrast = c('com_group', comp_group1, comp_group2), type = "ashr",
                                        quiet = TRUE))
    
    # COMPUTE NORMALIZED COUNTS
    norm_counts <- data.frame(DESeq2::counts(de_out, normalized = TRUE))
    
    res_counts <- cbind(res, norm_counts)
    res_counts <- res_counts[order(res_counts[[order_results]], decreasing = F),]
    
    
    # SAVE RESULTS
    if(! is.null(saveName)){
        write.csv(res_counts, file = file.path(outputDir, paste0('log2fc_', saveName,'.csv')), quote = F)
    } else {
        write.csv(res_counts, file = file.path(outputDir, 'log2fc.csv'), quote = F)
    }
    
    if(saverds){
        if(! is.null(saveName)){
            saveRDS(de_out, file = file.path(outputDir, paste0('DESeq2object_', saveName,'.RDS')))
        } else {
            saveRDS(de_out, file = file.path(outputDir, 'DESeq2object.RDS'))
        }
    }
    
    if(return_results){
        return(list(res_counts, de_out))
    }
}
