library(tidyverse)
library(limma)
gse_limma_diff = function(TFA,clinic){
  
  id = clinic$geo_accession
  TFA = TFA[,id] 
  
  gender = clinic$Sex.ch1
  gender[gender=="Male"] = "male"
  gender[gender=="Female"] = "female"
  gender = factor(gender,levels = c("female","male"))
  
  smoke = clinic$smoking_history.ch1
  smoke[smoke=="Never smoked"] = "No"
  smoke[smoke=="Currently smoking" | smoke=="Smoked in the past"] = "Yes"
  smoke[smoke!="Yes" & smoke!="No"] = "Unknown"
  smoke = factor(smoke)
  
  race = clinic$race.ch1
  race[race=="Black or African American"] = "black"
  race[race=="White"] = "white"
  race[race!="black" & race!="white"] = "other"
  race = factor(race)
  
  stage = ifelse(str_detect(clinic$disease_stage.ch1,"T1"),"I",
                 ifelse(str_detect(clinic$disease_stage.ch1,"T2"),"II",
                        ifelse(str_detect(clinic$disease_stage.ch1,"T3"),"III",
                               ifelse(str_detect(clinic$disease_stage.ch1,"T4"),"IV","Unknown"))))
  stage = factor(stage)
  
  
  age = clinic$age.ch1
  
  ## limma analysis
  design = model.matrix(~0+ gender + race + age + smoke + stage)
  TFA_fit = lmFit(TFA,design,method = "robust")
  #TFA_fit = lmFit(TFA,design)
  contrast = makeContrasts(FvsM = genderfemale - gendermale, levels = colnames(design))
  TFA_fit = contrasts.fit(TFA_fit, contrast)
  TFA_fit = eBayes(TFA_fit)
  DATs = topTable(TFA_fit,number=Inf,coef = "FvsM",sort.by = "p")
  cat("number of p<0.05:", sum(DATs$adj.P.Val<0.05))
  return(DATs)
}

#### 0. clinical
clinic = readRDS("data/GSE68465/pheno.rds")
gender = clinic$Sex.ch1
gender[gender=="Female"] = "female"
gender[gender=="Male"] = "male"
gender = factor(gender,levels = c("female","male"))

smoke = clinic$smoking_history.ch1
smoke[smoke=="Unknown" | smoke=="--"] = "Unknown"
smoke[smoke=="Never smoked"] = "No"
smoke[smoke!="No" & smoke!="Unknown"] = "Yes"
smoke = factor(smoke)

table(gender,smoke)

#### 1. DIFF TFA analysis
tiger = readRDS("result/pp/GSE68465_tiger.rds")
TFA = tiger$Z
TFA = TFA[,clinic$geo_accession]
DATs = gse_limma_diff(TFA,clinic)
DATs$"gene symbol" = rownames(DATs)
saveRDS(DATs,"result/pp/GSE68465_dats_limma.rds")

#### 2. DIFF expr analysis
DEGs = gse_limma_diff(expr,clinic)
DEGs$"gene symbol" = rownames(DEGs)
saveRDS(DEGs,"result/pp/GSE68465_degs_limma.rds")

