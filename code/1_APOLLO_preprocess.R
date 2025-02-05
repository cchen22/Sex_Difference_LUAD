library(tidyverse)
library(GenomicDataCommons)
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

info = cases() %>% 
  GenomicDataCommons::filter(~ project.project_id == "APOLLO-LUAD")%>%
  GenomicDataCommons::select('submitter_id')%>%
  results_all()
info = data.frame("case_id"=info$case_id,"submitter_id"=info$submitter_id)

## read in proteomics data
pro = read_csv("data/apollo/APOLLO1_Level3_Gprot_v061418.n87.csv",
               col_select = -(1:3))
pro=pro%>%
  group_by(Gene) %>%
  summarise(across(everything(), ~ mean(., na.rm = TRUE)))%>% # duplicate gene take mean
  dplyr::filter(!is.na(Gene)) %>%
  arrange(Gene) %>%
  column_to_rownames(var = "Gene")
dim(pro) #7614*87

## convert patient data to correct format 
patient.id = colnames(pro)
patient.id = gsub("\\.", "-", patient.id)
patient.id = gsub("^(.*?-.*?)(-.*)?$", "\\1", patient.id)
colnames(pro) = patient.id
saveRDS(pro,"data/pp/apollo_pro.rds")

## get clinical data from GDC
UUID = plyr::mapvalues(patient.id,info$submitter_id,info$case_id)
clinic = gdc_clinical(UUID, include_list_cols = FALSE)
clinical_data = inner_join(clinic$demographic,clinic$exposures,by="case_id") %>%
  dplyr::select("case_id","gender","tobacco_smoking_status") %>%
  inner_join(info,by="case_id")
saveRDS(clinical_data,"data/pp/apollo_clinic.rds")


## non-smoker number is too small so merge them?
table(clinical_data$tobacco_smoking_status,clinical_data$gender)


### read in official impuation data
phos = read.table("data/apollo/PHOSPHO_MS_RPPA_combined_non_global_norm_KNNImpute_K_10_measFrac_0.5_skip_RPPA_multimap.txt",
                  header = T,
                  check.names = F) %>%
  group_by(Phosphosite) %>%
  summarise(across(everything(), ~ mean(., na.rm = TRUE)))%>% # duplicate gene take mean
  arrange(Phosphosite) %>%
  column_to_rownames(var = "Phosphosite")
saveRDS(phos,"data/pp/apollo_phos_imputed.rds")




