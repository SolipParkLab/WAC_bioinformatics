## Loading libraries
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(ggplot2)
library(pheatmap)


### ... Loading tumor datasets ----
# Recommended to create one folder per cancer type.
cancer_types_tumor <- c('BRCA', 'CCRCC', 'COAD', 'GBM', 'HNSCC',
                        'LSCC', 'LUAD', 'OV', 'PDAC', 'UCEC')
# Loading tumor datasets
datasets_tumor <- lapply(cancer_types_tumor, function(type){
  dir <- sprintf('./data/%s/%s_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt', type, type)
  return(read.delim(dir))
})
names(datasets_tumor) <- cancer_types_tumor


### ... Loading normal datasets ----
# Recommended to create one folder per cancer type.
cancer_types_normal <- c('CCRCC', 'HNSCC', 'LSCC', 'LUAD', 'PDAC')
# Loading normal datasets
datasets_normal <- lapply(cancer_types_normal, function(type){
  dir <- sprintf('./data/%s/%s_RNAseq_gene_RSEM_coding_UQ_1500_log2_Normal.txt', type, type)
  return(read.delim(dir))
})
names(datasets_normal) <- cancer_types_normal


### ... Creating gene list ----
# Genes we will analyse and their matching Ensembl IDs
genes <- data.frame('Gene_Symbol' = c('WAC',
                                      'MTOR', 'RPTOR',
                                      'MLST8',
                                      'TELO2', 'TTI1', 'TTI2',
                                      'RUVBL1', 'RUVBL2', 'RPAP3', 'PIH1D1'),
                    'Ensembl_ID' = c('ENSG00000095787.23',
                                     'ENSG00000198793.13', 'ENSG00000141564.15',
                                     'ENSG00000167965.18',
                                     'ENSG00000100726.15', 'ENSG00000101407.13', 'ENSG00000129696.13',
                                     'ENSG00000175792.12', 'ENSG00000183207.14', 'ENSG00000005175.10', 'ENSG00000104872.11'))


### ... Filtering datasets to select our genes ----
# Tumor
filtered_datasets_tumor <- lapply(datasets_tumor, function(dataset){
  df <- filter(dataset, idx %in% genes$Ensembl_ID)
  df <- merge(genes, df,
              by.x = 'Ensembl_ID',
              by.y = 'idx')
  return(df)
})
names(filtered_datasets_tumor) <- cancer_types_tumor
# Normal
filtered_datasets_normal <- lapply(datasets_normal, function(dataset){
  df <- filter(dataset, idx %in% genes$Ensembl_ID)
  df <- merge(genes, df,
              by.x = 'Ensembl_ID',
              by.y = 'idx')
  return(df)
})
names(filtered_datasets_normal) <- cancer_types_normal


### ... Preparing input for heatmap ----
# Tumor
input_tumor <- lapply(cancer_types_tumor, function(cancer_type){
  df <- filtered_datasets_tumor[[cancer_type]]
  rownames(df) <- df$Gene_Symbol
  df[,c('Gene_Symbol', 'Ensembl_ID')] <- NULL
  df[, cancer_type] <- rowMeans(df, na.rm = T)
  df$Gene_Symbol <- row.names(df)
  return(df[,tail(seq_along(df), 2)])
})
input_tumor <- input_tumor %>% reduce(full_join, by = 'Gene_Symbol')
# Normal
input_normal <- lapply(cancer_types_normal, function(cancer_type){
  df <- filtered_datasets_normal[[cancer_type]]
  rownames(df) <- df$Gene_Symbol
  df[,c('Gene_Symbol', 'Ensembl_ID')] <- NULL
  df[, cancer_type] <- rowMeans(df, na.rm = T)
  df$Gene_Symbol <- row.names(df)
  return(df[,tail(seq_along(df), 2)])
})
input_normal <- input_normal %>% reduce(full_join, by = 'Gene_Symbol')
input_normal[,c('BRCA', 'GBM', 'COAD', 'OV', 'UCEC')] <- NA # Adding them because the have no normal samples


### ... Heatmap tumor only ----
t_only <- input_tumor
rownames(t_only) <- t_only$Gene_Symbol
t_only$Gene_Symbol <- NULL
t_only <- t_only[,sort(colnames(t_only))]
t_only <- t_only[genes$Gene_Symbol,]
t_only <- t(t_only)
# Computing Z-score
Z <- apply(t_only, 2, function(col){
  media <- mean(col, na.rm = T)
  desviacion <- sd(col, na.rm = T)
  return((col - media) / desviacion)
})
# Z-score heatmap
pheatmap(t(Z)[c(2:11, 1),],
         breaks = seq(-2.5, 2, 0.1),
         legend_breaks = seq(-2, 2, 1),
         color = colorRampPalette(c('blue', 'white', 'red'))(n=46),
         border_color = 'black',
         scale = 'none',
         cluster_rows = F,
         cluster_cols = T,
         na_col = 'lightgrey',
         main = 'RNA expression in tumor (Z score)')
# RNA expression heatmap
pheatmap(t(t_only[,c(2:11, 1)]),
         breaks = seq(0, 13, 0.1),
         legend_breaks = seq(0, 12, 2),
         color = colorRampPalette(c('white', '#feefac', 'yellow', 'red'))(n=135),
         border_color = 'black',
         scale = 'none',
         cluster_rows = F,
         cluster_cols = T,
         na_col = 'lightgrey',
         main = 'RNA expression in tumor')


### ··· Correlations ----
# Tumor
corr_tumor <- lapply(cancer_types_tumor, function(cancer_type){
  df <- filtered_datasets_tumor[[cancer_type]]
  rownames(df) <- df$Gene_Symbol
  df[,c('Gene_Symbol', 'Ensembl_ID')] <- NULL
  
  correlations <- expand.grid('Gene1' = rownames(df), 'Gene2' = rownames(df))
  correlation_result <- apply(correlations, 1, function(row){
    pval <- cor.test(as.numeric(df[row['Gene1'],]), as.numeric(df[row['Gene2'],]), method = 'pearson')$p.value
    estimate <- cor.test(as.numeric(df[row['Gene1'],]), as.numeric(df[row['Gene2'],]), method = 'pearson')$estimate
    return(c('Pvalue' = pval, 'Estimate' = estimate))
  })
  correlation_result <- as.data.frame(t(correlation_result))
  correlations <- cbind(correlations, correlation_result)
  colnames(correlations) <- c('Gene_1', 'Gene_2', paste0(cancer_type,'_Pval'), paste0(cancer_type,'_Estimate'))
  return(correlations)
})
corr_tumor <- corr_tumor %>% reduce(full_join, by = c('Gene_1', 'Gene_2'))


### ... Heatmaps for significant correlations in tumor ----
nsig <- corr_tumor %>% 
  pivot_longer(cols = colnames(corr_tumor)[!colnames(corr_tumor) %in% c('Gene_1', 'Gene_2')],
               values_to = 'Value',
               names_to = c('Canc_Type', 'Value_Type'),
               names_sep = '_')
nsig <- pivot_wider(nsig,
                    names_from = 'Value_Type',
                    values_from = 'Value')
nsig <- filter(nsig, Gene_1 != Gene_2)
# Selecting positive correlations
POSITIVE <- nsig %>% 
  na.omit() %>% 
  group_by(Gene_1, Gene_2) %>% 
  summarise('N_Sig_Positive' = sum(Pval <= 0.01 & Estimate > 0.2))
POSITIVE <- pivot_wider(POSITIVE,
                        values_from = 'N_Sig_Positive',
                        names_from = 'Gene_2')
POSITIVE <- as.data.frame(POSITIVE)
rownames(POSITIVE) <- POSITIVE$Gene_1
POSITIVE$Gene_1 <- NULL
POSITIVE <- POSITIVE[genes$Gene_Symbol, genes$Gene_Symbol]
# Heatmap
pheatmap(POSITIVE,
         legend_breaks = seq(0, 10, 2),
         color = colorRampPalette(c('white', '#ECC30B', '#FF1010'))(n=110),
         border_color = 'black',
         scale = 'none',
         cluster_rows = T,
         cluster_cols = T,
         na_col = 'darkgrey')
