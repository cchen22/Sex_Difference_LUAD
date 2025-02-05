## ----------------------------------------------------------
## Purpose of script: Data preprocess and TIGER analysis
## Author: Chen Chen
## Date Created: 2023-06-13
## Copyright (c) Chen Chen, 2023
## Email: chenchen9945@gmail.com
## -----------------------------------------------------------


library(tidyverse)
library(recount3)


# function to convert rownames from ensembl to gene names
# if symbol duplicated, take average
# remove comptelyte NA rows
convert_exprmat_ensembl2symbol = function(expression_matrix,ensembl_to_symbol){
  converted_matrix = expression_matrix %>%
    as.data.frame() %>%
    rownames_to_column(var = "ensembl_id") %>%
    left_join(ensembl_to_symbol, by = "ensembl_id") %>%
    group_by(gene_symbol) %>%
    summarise(across(-ensembl_id, ~ mean(., na.rm = TRUE)))%>%
    dplyr::filter(gene_symbol!="") %>%
    arrange(gene_symbol) %>%
    column_to_rownames(var = "gene_symbol")
  
  converted_matrix = converted_matrix[rowSums(is.na(converted_matrix)) != ncol(converted_matrix), ]
  
  return(converted_matrix)
}

# function for PCA plot
make_pca_plot = function(expr,title,group){
  pca_result = prcomp(t(expr), scale. = TRUE)
  # Create a data frame with PCA scores and batch information
  pca_df = data.frame(
    Sample = colnames(expr),
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    group = group
  )
  p = ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
    geom_point() +
    labs(title = title, x = "PC1", y = "PC2")
  print(p)
  return(p)
}


## data download
cancer = "LUAD"
library(recount3)
library(recount)
TCGA_lung <- recount3::create_rse_manual(
  project = cancer,
  project_home = "data_sources/tcga",
  organism = "human",
  annotation = "gencode_v26",
  type = "gene"
)

assays(TCGA_lung)$counts = recount3::transform_counts(TCGA_lung)
assays(TCGA_lung)$TPM = recount::getTPM(TCGA_lung)


# I. read RSE objects: LUAD -----------------------------------
## RSE data
rse = readRDS("data/tcga_luad.rds")

# II. Preprocessing 
## 1. filter samples and map TCGA barcode; 
## 2. normalize/filter genes
## 3. batch effects and outliers: PCA batch; PCA check outlier samples;
## 4. map gene names to symbol

## 1. get counts 
assays(rse)$counts = recount3::transform_counts(rse)

## get log2 tpm values
assays(rse)$TPM = recount::getTPM(rse)
assays(rse)$expr = log2(assays(rse)$TPM+1)

## TCGA includes normal tissue, remove them from analysis
table(rse$tcga.gdc_cases.samples.sample_type)
tumor_s = rse$external_id[rse$tcga.gdc_cases.samples.sample_type == "Primary Tumor"]
rse = rse[,tumor_s]#63856*540
dim(rse)

## map UUID to TCGA barcodes; 
## barcode has duplicates keep the one with less zero
barcode = colData(rse)$tcga.cgc_sample_id
colnames(rse) = barcode
patientID = substr(barcode,1,12)
patient.dup = patientID[duplicated(patientID)]
remove = NULL
for (i in 1:length(patient.dup)){
  tmp = which(patientID==patient.dup[i])
  num.zero = colSums(assays(rse)$counts[,tmp]==0)
  remove = c(remove,tmp[-which.min(num.zero)])
}
rse = rse[,-remove]
dim(rse)

## 2. filter low expressed genes
## at least expressed 10 in %20 samples
keep = edgeR::filterByExpr(assays(rse)$counts, design = NULL,
                           group = NULL, lib.size = NULL,
                           min.count = 10, min.prop = 0.2)
sum(keep)
rse = rse[keep,]
dim(rse)


# Perform PCA on the batch effects
batch_info = as.factor(colData(rse)$tcga.cgc_case_batch_number)
pca_plot = make_pca_plot(assays(rse)$expr,
                         title = "LUAD Batch Information", 
                         group = batch_info)
ggsave("plot/pca/luad_pca_batch_info.png",pca_plot,width = 6,height = 6)

## remove outlier samples below line y = 0.8x-200
pca_plot = pca_plot + geom_abline(slope = -0.6, intercept = 200, color = "black", linetype = "dashed") 
ggsave("plot/pca/luad_pca_batch_outlier.png",pca_plot,width = 6,height = 6)

pca_result = prcomp(t(assays(rse)$expr), scale. = TRUE)
pc_scores = as.data.frame(pca_result$x)

curve_PC2 = -0.6 *pc_scores$PC1 + 200
below_curve = rownames(pc_scores[pc_scores$PC2<curve_PC2,])
rse = rse[,below_curve]
dim(rse)

## 3. PCA on Y chromosome genes
chrY_expr = assays(rse[rowRanges(rse)@seqnames=="chrY",])$expr
dim(chrY_expr)
gender = rse@colData$tcga.gdc_cases.demographic.gender
pca_plot = make_pca_plot(chrY_expr,
                         title = "LUAD chrY 34 genes", 
                         group = gender)
ggsave("plot/pca/luad_pca_gender_info.png",pca_plot,width = 6,height = 6)

## add a verical cutoff
pca_plot= pca_plot + geom_vline(xintercept = -2.5, linetype = "dashed", color = "black")
ggsave("plot/pca/luad_pca_gender_cut.png",pca_plot,width = 6,height = 6)

## filter unclear male and female samples
pca_result = prcomp(t(chrY_expr), scale. = TRUE)
pc_scores = as.data.frame(pca_result$x)
misclassified_f = rownames(pc_scores[pc_scores$PC1>=-2.5 & gender=="female",])
misclassified_m = rownames(pc_scores[pc_scores$PC1<(-2.5) & gender=="male",])
misclassified_samples = c(misclassified_f,misclassified_m)
rse = rse[,!colnames(rse)%in%misclassified_samples]
dim(rse)


## 4. gene name mapping 
ensembl_to_symbol = data.frame("ensembl_id" = rowRanges(rse)$gene_id,
                               "gene_symbol" = rowRanges(rse)$gene_name)
# map rownames from ensemble to gene symbol. 
expr = convert_exprmat_ensembl2symbol(assays(rse)$expr,ensembl_to_symbol)
dim(expr)
which(rownames(expr)=="SRY") # ENSG00000184895.7 #58 >0 (128 counts)
which(rownames(expr)=="HSFY2") # ENSG00000169953.11 #151 >0 (1007 counts)


# III. save expression mats: LUAD, LUSC, GTEx -----------------------------

## save expression matrix
saveRDS(expr,"data/pp/luad_expression.rds")
saveRDS(rse,"data/pp/luad_rse.rds")

