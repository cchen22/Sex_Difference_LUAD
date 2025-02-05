library(limma)
library(edgeR)
library(SummarizedExperiment)
TCGA_deg_function = function(rse){
  counts = assays(rse)$counts
  dge = DGEList(counts = counts)
  dge = calcNormFactors(dge,method = "TMM")
  
  ## load TCGA clinic variables
  metadata=rse
  id = metadata@colData$tcga.cgc_sample_id
  gender = metadata@colData$tcga.gdc_cases.demographic.gender[match(colnames(counts),id)]
  race = metadata@colData$tcga.gdc_cases.demographic.race[match(colnames(counts),id)]
  age = metadata@colData$tcga.xml_age_at_initial_pathologic_diagnosis[match(colnames(counts),id)]
  stage = metadata@colData$tcga.gdc_cases.diagnoses.tumor_stage[match(colnames(counts),id)]
  smoke = metadata@colData$tcga.xml_tobacco_smoking_history[match(colnames(counts),id)]
  
  gender = factor(gender,levels = c("female","male"))
  
  #race[race=="american indian or alaska native"] = "native"
  race[race=="black or african american"] = "black"
  race[race=="white"] = "white"
  race[race!="black" & race!="white"] = "other"
  race = factor(race)
  
  age[which(is.na(age))] = mean(age,na.rm=TRUE)
  
  stage[stage=="stage i" | stage=="stage ia" | stage == "stage ib"] = "I"
  stage[stage=="stage ii" | stage=="stage iia" | stage == "stage iib"] = "II"
  stage[stage=="stage iii" |stage=="stage iiia" | stage == "stage iiib"] = "III"
  stage[stage=="stage iv"] = "IV"
  stage[stage=="not reported"] = "Unknown"
  stage = factor(stage)
  
  smoke[is.na(smoke)] = "Unknown"
  smoke[smoke=="1"] = "No"
  smoke[smoke!="No" & smoke!="Unknown"] = "Yes"
  smoke = factor(smoke)
  
  
  design = model.matrix(~0+ gender + race + age + stage + smoke)
  v = voom(dge, design, plot=TRUE)
  DEG_fit = lmFit(v,design)
  contrast = makeContrasts(FvsM = genderfemale - gendermale, levels = colnames(design))
  DEG_fit = contrasts.fit(DEG_fit, contrast)
  DEG_fit = eBayes(DEG_fit)
  DEGs = topTable(DEG_fit,number=Inf,coef = "FvsM",sort.by = "p",p.value = 1,lfc = 0)
  
  gene_info = data.frame("ensembl_id" = rowRanges(rse)$gene_id,
                         "gene_symbol" = rowRanges(rse)$gene_name,
                         "chr_info" = rowRanges(rse)@seqnames)
  DEGs = DEGs %>%
    rownames_to_column(var = "ensembl_id") %>%
    left_join(gene_info, by = "ensembl_id") %>%
    dplyr::filter(gene_symbol!="")
  
  return(DEGs)
  
}
TCGA_limma_diff = function(TFA,metadata){
  require(limma)
  ## load TCGA clinic variables
  id = metadata@colData$tcga.cgc_sample_id
  gender = metadata@colData$tcga.gdc_cases.demographic.gender[match(colnames(TFA),id)]
  race = metadata@colData$tcga.gdc_cases.demographic.race[match(colnames(TFA),id)]
  age = metadata@colData$tcga.xml_age_at_initial_pathologic_diagnosis[match(colnames(TFA),id)]
  stage = metadata@colData$tcga.gdc_cases.diagnoses.tumor_stage[match(colnames(TFA),id)]
  smoke = metadata@colData$tcga.xml_tobacco_smoking_history[match(colnames(TFA),id)]
  
  gender = factor(gender,levels = c("female","male"))
  
  #race[race=="american indian or alaska native"] = "native"
  race[race=="black or african american"] = "black"
  race[race=="white"] = "white"
  race[race!="black" & race!="white"] = "other"
  race = factor(race)
  
  age[which(is.na(age))] = mean(age,na.rm=TRUE)
  
  stage[stage=="stage i" | stage=="stage ia" | stage == "stage ib"] = "I"
  stage[stage=="stage ii" | stage=="stage iia" | stage == "stage iib"] = "II"
  stage[stage=="stage iii" |stage=="stage iiia" | stage == "stage iiib"] = "III"
  stage[stage=="stage iv"] = "IV"
  stage[stage=="not reported"] = "Unknown"
  stage = factor(stage)
  
  smoke[is.na(smoke)] = "Unknown"
  smoke[smoke=="1"] = "No"
  smoke[smoke!="No" & smoke!="Unknown"] = "Yes"
  smoke = factor(smoke)
  
  ## limma analysis
  design = model.matrix(~0+ gender + race + age + stage + smoke)
  #TFA_fit = lmFit(TFA,design,method = "robust")
  TFA_fit = lmFit(TFA,design)
  contrast = makeContrasts(FvsM = genderfemale - gendermale, levels = colnames(design))
  TFA_fit = contrasts.fit(TFA_fit, contrast)
  TFA_fit = eBayes(TFA_fit)
  DATs = topTable(TFA_fit,number=Inf,coef = "FvsM",sort.by = "p")
  cat("number of p<0.05:", sum(DATs$adj.P.Val<0.05))
  return(DATs)
}

##### part1. luad DEG using limma-voom
rse = readRDS("data/pp/luad_rse.rds")
DEGs = TCGA_deg_function(rse)
saveRDS(DEGs,"result/pp/luad_degs.rds")

#### part2. luad TFA difference
tiger = readRDS("result/pp/luad_tiger.rds")
TFA = tiger$Z #734*502
rse = readRDS("data/pp/luad_rse.rds")
identical(colnames(TFA),colnames(rse))

DATs = TCGA_limma_diff(TFA,rse)
DATs$"gene symbol" = rownames(DATs)
saveRDS(DATs,"result/pp/luad_dats_limma.rds")




