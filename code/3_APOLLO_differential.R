library(limma)
library(edgeR)
apollo_limma_diff = function(pro,clinic){
  
  id = clinic$submitter_id
  pro = pro[,id] 
  
  gender = clinic$gender
  gender = factor(gender,levels = c("female","male"))
  
  smoke = clinic$tobacco_smoking_status
  smoke[smoke=="Lifelong Non-Smoker"] = "No"
  smoke[smoke!="No"] = "Yes"
  smoke = factor(smoke)
  
  ## limma analysis
  design = model.matrix(~0+ gender + smoke)
  pro_fit = lmFit(pro,design)
  contrast = makeContrasts(FvsM = genderfemale - gendermale, levels = colnames(design))
  pro_fit = contrasts.fit(pro_fit, contrast)
  pro_fit = eBayes(pro_fit)
  DE = topTable(pro_fit,number=Inf,coef = "FvsM",sort.by = "p")
  cat("number of p<0.05:", sum(DE$adj.P.Val<0.05))
  return(DE)
}

### 0. clinical data
clinic = readRDS("data/pp/apollo_clinic.rds")

### 1. differential protein expression
pro = readRDS("data/pp/apollo_pro.rds")
DEPs = apollo_limma_diff(pro,clinic)
DEPs$`gene symbol` =rownames(DEPs)
saveRDS(DEPs,"result/pp/apollo_deps.rds")

### 2. differential phosphorylation
phos = readRDS("data/pp/apollo_phos_imputed.rds")
DPPs = apollo_limma_diff(phos,clinic)
DPPs$"gene symbol" = sub("_.*", "", rownames(DPPs))
saveRDS(DPPs,"result/pp/apollo_dpps.rds")


#### 3. PTM-SEA results diff
PTMgsea = cmapR::parse.gctx("data/PTM-SEA/apollo_phos-scores.gct")
res = PTMgsea@mat
res_diff = apollo_limma_diff(res,clinic)
res_diff$"gene symbol" = sub("^[^_]*_", "", rownames(res_diff))
res_diff$"class" = sapply(strsplit(rownames(res_diff), "-"), "[[", 1)
saveRDS(res_diff,"result/pp/apollo_dPTMgsea_limma.rds")


