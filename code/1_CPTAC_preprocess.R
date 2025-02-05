
# Load the package
library(reticulate)

# filter duplicate rownames by less number of NAs
# if same number of NAs than by larger rowSums
select_and_filter_rows <- function(mat) {
  row_missing_counts <- rowSums(is.na(mat))
  unique_row_names <- unique(rownames(mat))
  select_rows <- logical(nrow(mat))
  for (name in unique_row_names) {
    rows_with_name <- rownames(mat) == name
    min_missing_count <- min(row_missing_counts[rows_with_name])
    if (sum(rows_with_name) > 1) {
      select_rows[which(rows_with_name)[which.min(row_missing_counts[rows_with_name])]] <- TRUE
    } else {
      select_rows[rows_with_name] <- TRUE
    }
  }
  filtered_matrix <- mat[select_rows, , drop = FALSE]
  return(filtered_matrix)
}

########## data download #########
cptac = import("cptac", convert = FALSE) 
utils = import("cptac.utils", convert = FALSE)
cc = cptac$Luad()

# proteomics data
pro_py = cc$get_proteomics() 
pro_py = utils$reduce_multiindex(pro_py, flatten = TRUE) 
pro = py_to_r(pro_py) 
pro = t(pro)
rownames(pro) = sub("_.*", "", rownames(pro))
pro =  select_and_filter_rows(pro)
saveRDS(pro,"data/luad_proteomics.rds")

# phosphoproteomics data
phos_py = cc$get_phosphoproteomics()
phos_py = utils$reduce_multiindex(phos_py, flatten = TRUE) 
phos = py_to_r(phos_py) 
phos = t(phos)
rownames(phos) = sub("^(.*?_.*?)_.*", "\\1", rownames(phos))
phos =  select_and_filter_rows(phos)
saveRDS(phos,"data/luad_phosphoproteomics.rds")

# acetylproteomics data
ace_py = cc$get_acetylproteomics()
ace_py = utils$reduce_multiindex(ace_py, flatten = TRUE) 
ace = py_to_r(ace_py)
ace = t(ace)
ace[ace=="NaN"]=NA
ace = apply(ace, c(1, 2), as.numeric)
rownames(ace) = sub("^(.*?_.*?)_.*", "\\1", rownames(ace))
ace =  select_and_filter_rows(ace)
saveRDS(ace,"data/luad_acetylproteomics.rds")

# clinical data
clinic_py = cc$get_clinical()
clinic = py_to_r(clinic_py) 
dim(clinic)
saveRDS(clinic,"data/luad_cptac_clinic.rds")


######### data preprocessing ################
library(impute)
clinic =readRDS("data/luad_cptac_clinic.rds")

### proteomics
pro = readRDS("data/luad_proteomics.rds")
pro = pro[,rownames(clinic)] # remove normal samples ".N"
pro = pro[rowSums(is.na(pro))<ncol(pro)/2,]

set.seed(123)
pro.imputed = impute.knn(pro,k=5, rowmax = 0.5, colmax = 0.8)
pro = pro.imputed$data
saveRDS(pro,"data/luad_proteomics_imputed.rds")


### phosphorylation
phos = readRDS("data/luad_phosphoproteomics.rds")
phos = phos[,rownames(clinic)] # remove normal samples ".N"
phos = phos[rowSums(is.na(phos))<ncol(phos)/2,] #50% not missing

set.seed(123)
phos.imputed = impute.knn(phos,k=5)
phos = phos.imputed$data
saveRDS(phos,"data/luad_phosphoproteomics_imputed.rds")

### acetylation
ace = readRDS("data/luad_acetylproteomics.rds")
ace = ace[,rownames(clinic)] # remove normal samples ".N"
ace = ace[rowSums(is.na(ace))<ncol(ace)/2,]

set.seed(123)
ace.imputed = impute.knn(ace,k=5, rowmax = 0.5, colmax = 0.5)
ace = ace.imputed$data
saveRDS(ace,"data/luad_acetylproteomics_imputed.rds")




