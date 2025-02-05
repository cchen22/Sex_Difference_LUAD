library(tidyverse)
library(ggplot2)
library(recount)
library(data.table)
library(sva)

#### part1. load data
# expression data
expression = data.frame(fread("data/GSE68465/GSE68465_expression.txt"))
expression[, -c(1,2)] = log2(1+expression[, -c(1,2)])

# phenotype data
phenotypes = data.frame(fread("data/GSE68465/GSE68465_phenotype.txt"))

#### part2. batch correction
# check batch effects
sample_IDs = colnames(expression)[-c(1,2)]
sample_source = phenotypes$source_name_ch1[match(sample_IDs, phenotypes$geo_accession)]
pca = prcomp(t(data.matrix(expression[,-c(1,2)])))
U = pca$x
plot(U[, 1], U[, 2], main = "PCA of GSE68465", col = factor(sample_source), pch = 19, xlab = "PC1", ylab = "PC2")

# Batch Correction
expression_corrected <- ComBat(expression[,-c(1,2)], batch = sample_source)
rownames(expression_corrected) = expression[,1]
pca_corrected = prcomp(t(data.matrix(expression_corrected)))
U <- pca_corrected$x
plot(U[, 1], U[, 2], main = "PCA of GSE68465 (batch-corrected)", col = factor(sample_source), pch = 19, xlab = "PC1", ylab = "PC2")

expression = cbind(expression[,1:2],expression_corrected)

#### part3. map probe to unique gene name
# remove empty gene names
# keep the gene names shown in dorothea targets
# For every gene symbol, compute the standard deviation (SD) for all the probes. 
# Keep the probe with highest SD and discard the rest.

expression = expression[expression$gene_name!="",]
gene_symbols = expression$gene_name
gene_symbols2 = strsplit(gene_symbols, split = " /// ", fixed= T)
target_genes = unique(dorothea::dorothea_hs_pancancer$target)

# Function to filter genes
filter_genes = function(genes, filter) {
  common_genes = intersect(genes, filter)
  if (length(common_genes) > 0) {
    return(common_genes[1])
  } else {
    return(genes[1])
  } 
}

# Apply the function to the list of genes
filtered_genes = lapply(gene_symbols2, filter_genes, filter = target_genes)
gene_symbols2 = unlist(filtered_genes)

probewise_sd = unlist(apply(expression[,-c(1:2)], MARGIN = 1, sd))
df = data.frame(cbind(gene_symbols2, probewise_sd, names(probewise_sd) ))
filtered = df %>% 
  group_by(gene_symbols2) %>% 
  dplyr::slice(which.max(probewise_sd))
genes_filtered = filtered$gene_symbols2
probes_filtered = filtered$V3
expression_filtered = expression[which(rownames(expression) %in% probes_filtered),]
rownames(expression_filtered) = genes_filtered[match(rownames(expression_filtered), probes_filtered)]

expr = as.matrix(expression_filtered[,-c(1:2)])


#### part4. check Y chrom genes
# Get Y genes
# Get chromosome and gene names
info.TCGA = readRDS("data/tcga_luad.rds")
gene_info.TCGA = rowRanges(info.TCGA)
gene_info.TCGA = data.frame(gene_info.TCGA)
Y_genes = gene_info.TCGA$gene_name[which(gene_info.TCGA$seqnames == "chrY")]

# Check sex-misannotation
sum(is.na(phenotypes$Sex.ch1))
phenotypes = phenotypes[!is.na(phenotypes$Sex.ch1),]
gender = phenotypes$Sex.ch1
expr = expr[,phenotypes$geo_accession]

# PCA plot to detect sex mis-specification
Y_expr = expr[rownames(expr) %in% Y_genes,]
gender = factor(as.character(gender))
pca = prcomp(t(Y_expr))
U <- pca$x
theme_set(bigstatsr::theme_bigstatsr(0.8))
pca_gender = qplot(U[, 1], U[, 2], main = "PCA of Y genes: GSE68465", colour = gender) +
  coord_equal() +
  xlab("PC1") +
  ylab("PC2")
ggsave("plot/pca/pca_GSE68465_gender.png",pca_gender,width = 6,height = 4)

# Samples to exclude because of sex-misannotation
female_wrong = colnames(Y_expr)[which(U[,1] >0 & gender == "Female")]
male_wrong = colnames(Y_expr)[which(U[,1] <0 & gender == "Male")]

expr = expr[,-which(colnames(expr) %in% c(female_wrong, male_wrong))]
phenotypes = phenotypes[phenotypes$geo_accession%in%colnames(expr),]

# save
saveRDS(expr,"data/GSE68465/expr.rds")
saveRDS(phenotypes,"data/GSE68465/pheno.rds")



