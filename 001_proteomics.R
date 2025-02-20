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
  dir <- sprintf('./data/%s/%s_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt', type, type)
  return(read.delim(dir))
})
names(datasets_tumor) <- cancer_types_tumor


### ... Loading normal datasets ----
# Recommended to create one folder per cancer type.
cancer_types_normal <- c('CCRCC', 'COAD', 'HNSCC',
                        'LSCC', 'LUAD', 'OV', 'PDAC', 'UCEC')
# Loading normal datasets
datasets_normal <- lapply(cancer_types_normal, function(type){
  dir <- sprintf('./data/%s/%s_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt', type, type)
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
input_normal[,c('BRCA', 'GBM')] <- NA # Adding them because the have no normal samples


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
  print(media)
  print(desviacion)
  return((col - media) / desviacion)
})
# Z-score heatmap
pheatmap(t(Z[,c(2:11, 1)]),
         breaks = seq(-2.5, 2, 0.1),
         legend_breaks = seq(-2, 2, 1),
         color = colorRampPalette(c('blue', 'white', 'red'))(n=46),
         border_color = 'black',
         scale = 'none',
         cluster_rows = F,
         cluster_cols = T,
         na_col = 'lightgrey',
         main = 'Protein expression in tumor (Z-score)')
# Protein expression heatmap
pheatmap(t(t_only[,c(2:11, 1)]),
         breaks = seq(0, 28, 0.1),
         legend_breaks = seq(0, 28, 3),
         color = colorRampPalette(c('white', '#FFFFDC', '#feefac', 'yellow', '#FF3A3A'))(n=281),
         border_color = 'black',
         scale = 'none',
         cluster_rows = F,
         cluster_cols = T,
         na_col = 'lightgrey',
         main = 'Protein expression in tumor')


### ... Heatmap difference tumor-normal ----
# Creating normal only
n_only <- input_normal
rownames(n_only) <- n_only$Gene_Symbol
n_only$Gene_Symbol <- NULL
n_only <- n_only[,sort(colnames(n_only))]
n_only <- n_only[genes$Gene_Symbol,]
n_only <- t(n_only)
# P-values difference
pvals_t_vs_n <- lapply(cancer_types_normal, function(cancer_type){
  df_t <- filtered_datasets_tumor[[cancer_type]]
  rownames(df_t) <- df_t$Gene_Symbol
  df_t[,c('Gene_Symbol', 'Ensembl_ID')] <- NULL
  df_n <- filtered_datasets_normal[[cancer_type]]
  rownames(df_n) <- df_n$Gene_Symbol
  df_n[,c('Gene_Symbol', 'Ensembl_ID')] <- NULL
  wilcox_pvals <- unlist(lapply(rownames(df_t), function(gene){
    return(wilcox.test(as.numeric(df_t[gene,]), as.numeric(df_n[gene,]), alternative = 'two.sided')$p.value)
  }))
  return_df <- data.frame('Gene' = rownames(df_t),
                          'Cancer.Type' = rep(cancer_type, length(rownames(df_t))),
                          'Wilcox.Pval' = wilcox_pvals)
  return(return_df)
})
pvals_t_vs_n <- bind_rows(pvals_t_vs_n)
pvals_t_vs_n <- pivot_wider(pvals_t_vs_n,
                            names_from = 'Cancer.Type',
                            values_from = 'Wilcox.Pval')
pvals_t_vs_n[,c('BRCA', 'GBM')] <- NA
pvals_t_vs_n <- as.data.frame(pvals_t_vs_n)
rownames(pvals_t_vs_n) <- pvals_t_vs_n$Gene
# Mapping P-values into significance levels
txt <- as.data.frame(pvals_t_vs_n[genes$Gene_Symbol, c('Gene', cancer_types_tumor)])
row.names(txt) <- txt$Gene
txt$Gene <- NULL
txt <- ifelse(txt <= 0.001, '***',
              ifelse(txt <= 0.01, '**',
                     ifelse(txt <= 0.05, '*', '')))
txt <- t(txt)
txt[is.na(txt)] <- ''
# Heatmap
diff <- t_only - n_only
pheatmap(t(diff[,c(5:9, 1)]),
         breaks = seq(-0.5, 0.5, 0.01),
         legend_breaks = seq(-0.5, 0.5, 0.1),
         color = colorRampPalette(c('blue', 'white', 'red'))(n=100),
         border_color = 'black',
         scale = 'none',
         cluster_rows = F,
         cluster_cols = T,
         na_col = 'darkgrey',
         display_numbers = t(txt[,c(5:9, 1)]),
         number_color = 'black',
         main = 'Protein expression difference')


### ... Computing average expression, difference and P-value ----
# Tumor
all_samples_tumor <- filtered_datasets_tumor[cancer_types_normal] %>% reduce(full_join, by = c('Ensembl_ID', 'Gene_Symbol'))
rownames(all_samples_tumor) <- all_samples_tumor$Gene_Symbol
all_samples_tumor[,c("Gene_Symbol", "Ensembl_ID")] <- NULL
# Normal
all_samples_normal <- filtered_datasets_normal[cancer_types_normal] %>% reduce(full_join, by = c('Ensembl_ID', 'Gene_Symbol'))
rownames(all_samples_normal) <- all_samples_normal$Gene_Symbol
all_samples_normal[,c("Gene_Symbol", "Ensembl_ID")] <- NULL
# P-values
all_pvals <- unlist(lapply(rownames(all_samples_tumor), function(gene){
  return(wilcox.test(as.numeric(all_samples_tumor[gene,]), as.numeric(all_samples_normal[gene,]), alternative = 'two.sided')$p.value)
}))
all_pvals_df <- data.frame('Gene' = row.names(all_samples_tumor),
                           'Wilcox.Pval' = all_pvals)
rownames(all_pvals_df) <- all_pvals_df$Gene
all_pvals_df <- all_pvals_df[c(3, 4, 6, 9, 10, 2),]
all_pvals_df$Gene <- NULL
txt <- all_pvals_df
txt <- ifelse(txt <= 0.001, '***',
              ifelse(txt <= 0.01, '**',
                     ifelse(txt <= 0.05, '*', '')))
# Averages tumor
all_avg_tumor <- rowMeans(all_samples_tumor, na.rm = T)
all_avg_tumor <- data.frame('Gene' = row.names(all_samples_tumor),
                            'Avg.Pancancer' = all_avg_tumor)
all_avg_tumor <- all_avg_tumor[c(3, 4, 6, 9, 10, 2),]
all_avg_tumor$Gene <- NULL
# Averages normal
all_avg_normal <- rowMeans(all_samples_normal, na.rm = T)
all_avg_normal <- data.frame('Gene' = row.names(all_samples_normal),
                             'Avg.Pancancer' = all_avg_normal)
all_avg_normal <- all_avg_normal[c(3, 4, 6, 9, 10, 2),]
all_avg_normal$Gene <- NULL
# Difference
all_diff <- all_avg_tumor - all_avg_normal
pheatmap(all_diff,
         breaks = seq(-0.5, 0.5, 0.01),
         legend_breaks = seq(-0.5, 0.5, 0.1),
         color = colorRampPalette(c('blue', 'white', 'red'))(n=100),
         border_color = 'black',
         scale = 'none',
         cluster_rows = F,
         cluster_cols = F,
         na_col = 'darkgrey',
         display_numbers = txt,
         number_color = 'black',
         main = 'Protein expression difference',
         cellwidth = 40)
         

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


