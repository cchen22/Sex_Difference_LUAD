Proteomics_limma_diff = function(pro,clinic){
  require(limma)
  id = rownames(clinic)
  pro = pro[,id] # remove normal samples ".N"
  
  ## load TCGA clinic variables
  gender = clinic$sex[match(colnames(pro),id)]
  race = clinic$race[match(colnames(pro),id)]
  age = clinic$age[match(colnames(pro),id)]
  stage = clinic$tumor_stage_pathological[match(colnames(pro),id)]
  smoke = clinic$tobacco_smoking_history[match(colnames(pro),id)]
  
  gender[gender=="Female"] = "female"
  gender[gender=="Male"] = "male"
  gender = factor(gender,levels = c("female","male"))
  
  #race[race=="american indian or alaska native"] = "native"
  race[race=="Black or African American"] = "black"
  race[race=="White"] = "white"
  race[race!="black" & race!="white"] = "other"
  race = factor(race)
  
  age = as.numeric(age)
  age[which(is.na(age))] = mean(age,na.rm=TRUE)
  
  stage[stage=="Stage I"] = "I"
  stage[stage=="Stage II"] = "II"
  stage[stage=="Stage III"] = "III"
  stage[stage=="Stage IV"] = "IV"
  stage = factor(stage)
  
  smoke[smoke=="Smoking history not available"] = "Unknown"
  smoke[smoke=="Lifelong non-smoker: Less than 100 cigarettes smoked in lifetime"] = "No"
  smoke[smoke!="No" & smoke!="Unknown"] = "Yes"
  smoke = factor(smoke)
  
  ## limma analysis
  pro = pro[rowSums(is.na(pro))<(min(table(gender))/2),]# at least half female not NA
  
  design = model.matrix(~0+ gender + race + age + stage + smoke)
  pro_fit = lmFit(pro,design)
  contrast = makeContrasts(FvsM = genderfemale - gendermale, levels = colnames(design))
  pro_fit = contrasts.fit(pro_fit, contrast)
  pro_fit = eBayes(pro_fit)
  DE = topTable(pro_fit,number=Inf,coef = "FvsM",sort.by = "p")
  cat("number of p<0.05:", sum(DE$adj.P.Val<0.05))
  return(DE)
}
PTM_limma_diff = function(ptm,pro,clinic){
  perform_limma_analysis <- function(phosphosite, ptm_data, pro_data, design) {
    # Extract the protein expression data for the corresponding protein
    protein_id <- gsub("_.*", "", phosphosite)  # Extract protein name from phosphosite
    protein_expr <- pro_data[protein_id, , drop = FALSE]  # Subset the protein expression data
    
    # Update the design matrix with the protein expression as a covariate
    updated_design <- cbind(design, t(protein_expr))
    colnames(updated_design) = c(colnames(design),"protein")
    
    # Fit the linear model using limma on PTM data for the specific phosphosite
    y <- ptm_data[phosphosite, , drop = FALSE]  # Subset the PTM data
    fit <- lmFit(y, updated_design)
    
    # Define contrasts (assuming you're interested in gender comparison)
    contrast <- makeContrasts(FvsM = genderfemale - gendermale, levels = colnames(updated_design))
    fit <- contrasts.fit(fit, contrast)
    fit <- eBayes(fit)
    
    # Get differentially expressed results
    DE <- topTable(fit, number = Inf, coef = "FvsM", sort.by = "p")
    return(DE)
  }
  
  require(limma)
  
  id = rownames(clinic)
  pro = pro[,id] # remove normal samples ".N"
  ptm = ptm[,id] # remove normal samples ".N"
  
  ## load TCGA clinic variables
  gender = clinic$sex[match(colnames(pro),id)]
  race = clinic$race[match(colnames(pro),id)]
  age = clinic$age[match(colnames(pro),id)]
  stage = clinic$tumor_stage_pathological[match(colnames(pro),id)]
  smoke = clinic$tobacco_smoking_history[match(colnames(pro),id)]
  
  gender[gender=="Female"] = "female"
  gender[gender=="Male"] = "male"
  gender = factor(gender,levels = c("female","male"))
  
  #race[race=="american indian or alaska native"] = "native"
  race[race=="Black or African American"] = "black"
  race[race=="White"] = "white"
  race[race!="black" & race!="white"] = "other"
  race = factor(race)
  
  age = as.numeric(age)
  age[which(is.na(age))] = mean(age,na.rm=TRUE)
  
  stage[stage=="Stage I"] = "I"
  stage[stage=="Stage II"] = "II"
  stage[stage=="Stage III"] = "III"
  stage[stage=="Stage IV"] = "IV"
  stage = factor(stage)
  
  smoke[smoke=="Smoking history not available"] = "Unknown"
  smoke[smoke=="Lifelong non-smoker: Less than 100 cigarettes smoked in lifetime"] = "No"
  smoke[smoke!="No" & smoke!="Unknown"] = "Yes"
  smoke = factor(smoke)
  
  # Basic design matrix without protein expression
  basic_design <- model.matrix(~ 0 + gender + race + age + stage + smoke)
  
  # Loop over each phosphosite and perform analysis
  ptm = ptm[gsub("_.*", "", rownames(ptm))%in%rownames(pro), ]
  results <- lapply(rownames(ptm), function(phosphosite) {
    perform_limma_analysis(phosphosite, ptm, pro, basic_design)
  })
  final_results <- do.call(rbind, results)
  rownames(final_results) = rownames(ptm)
  
  final_results = final_results[order(final_results$adj.P.Val),]
  final_results$adj.P.Val = p.adjust(final_results$adj.P.Val) #not adjusted yet
  return(final_results)
}

#### part1. clinical data: gender and smoke
clinic = readRDS("data/luad_cptac_clinic.rds")
gender = clinic$sex
gender[gender=="Female"] = "female"
gender[gender=="Male"] = "male"
gender = factor(gender,levels = c("female","male"))

smoke = clinic$tobacco_smoking_history
smoke[smoke=="Smoking history not available"] = "Unknown"
smoke[smoke=="Lifelong non-smoker: Less than 100 cigarettes smoked in lifetime"] = "No"
smoke[smoke!="No" & smoke!="Unknown"] = "Yes"
smoke = factor(smoke)

table(gender,smoke)

#### part 2. Proteomics
pro = readRDS("data/luad_proteomics_imputed.rds")
DEPs = Proteomics_limma_diff(pro,clinic)
DEPs$`gene symbol` = rownames(DEPs)
saveRDS(DEPs,"result/pp/luad_deps_limma.rds")

#### part 3. Phosphorylation
phos = readRDS("data/luad_phosphoproteomics_imputed.rds")
DPPs = Proteomics_limma_diff(phos,clinic)
DPPs$"gene symbol" = sub("_.*", "", rownames(DPPs))
saveRDS(DPPs,"result/pp/luad_dpps_limma.rds")

DPPs2 = PTM_limma_diff(phos,pro,clinic)
DPPs2$"gene symbol" = sub("_.*", "", rownames(DPPs2))
saveRDS(DPPs2,"result/pp/luad_dpps2_limma.rds")

#### part 4. acetylation
ace = readRDS("data/luad_acetylproteomics_imputed.rds")
DAPs = Proteomics_limma_diff(ace,clinic)
DAPs$"gene symbol" = sub("_.*", "", rownames(DAPs))
saveRDS(DAPs,"result/pp/luad_daps_limma.rds")

DAPs2 = PTM_limma_diff(ace,pro,clinic)
DAPs2$"gene symbol" = sub("_.*", "", rownames(DAPs2))
saveRDS(DAPs2,"result/pp/luad_daps2_limma.rds")

#### part 5. kinase activity difference
phos_activity = readRDS("result/pp/phos_activity_tiger.rds")
phos_activity = phos_activity[,rownames(clinic)]
DAKs = Proteomics_limma_diff(phos_activity,clinic)
DAKs$"gene symbol" =rownames(DAKs)
saveRDS(DAKs,"result/pp/luad_daks_limma.rds")


#### part 6. PTM-SEA results
PTMgsea = cmapR::parse.gctx("data/PTM-SEA/luad_phos-scores.gct")
res = PTMgsea@mat
res =res[,rownames(clinic)]
res_diff = Proteomics_limma_diff(res,clinic)
res_diff$"gene symbol" = sub("^[^_]*_", "", rownames(res_diff))
res_diff$"class" = sapply(strsplit(rownames(res_diff), "-"), "[[", 1)
saveRDS(res_diff,"result/pp/luad_dPTMgsea_limma.rds")




