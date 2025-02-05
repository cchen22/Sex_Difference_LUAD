library(tidyverse)
library(glmnet)
library(ggplot2)
library(dplyr)
library(ggforce)

## lasso with bootsrap
# perform_lasso = function(Y, X) {
#   library(glmnet)
#   coefficients_matrix = matrix(NA, nrow = nrow(X), ncol = nrow(Y))
#   X = t(X)
#   for (i in 1:nrow(Y)) {
#     cv_fit = cv.glmnet(X, Y[i, ], alpha = 1,standardize=F,
#                        lower.limits = c(rep(0,15),rep(-Inf,13)),
#                        upper.limits = c(rep(Inf,15),rep(0,13)),)
#     beta = coef(cv_fit, s = "lambda.min")  
#     coefficients_matrix[,i] = as.vector(beta[-1])
#   }
#   return(coefficients_matrix)
# }
perform_lasso_with_bootstrap = function(Y, X, n_pos=15,n_neg=13,
                                        n_bootstrap = 100, 
                                        seed = 123) {
  library(glmnet)
  
  # standarize
  Y = t(scale(t(Y)))
  X = t(scale(t(X)))
  
  # seed
  set.seed(seed)
  
  # Transpose X for glmnet compatibility
  X = t(X)
  
  # Initialize matrices to store coefficients and their distributions
  coefficients_matrix = matrix(NA, nrow = ncol(X), ncol = nrow(Y))
  coefficients_distributions = array(NA, dim = c(ncol(X), nrow(Y), n_bootstrap))
  trinary_matrix = matrix(0, nrow = ncol(X), ncol = nrow(Y)) # Trinary matrix
  
  # Perform bootstrap for each row of Y
  for (i in 1:nrow(Y)) {
    for (j in 1:n_bootstrap) {
      # Bootstrap resampling
      set.seed(seed+j)
      sample_indices = sample(nrow(X), replace = TRUE)
      X_bootstrap = X[sample_indices, ]
      Y_bootstrap = Y[i, sample_indices]
      
      # Fit Lasso model on bootstrap sample
      #set.seed(seed+i)
      cv_fit = cv.glmnet(X, Y[i, ], alpha = 1,standardize=F,
                         lower.limits = c(rep(0,n_pos),rep(-Inf,n_neg)),
                         upper.limits = c(rep(Inf,n_pos),rep(0,n_neg)),)
      beta = coef(cv_fit, s = "lambda.min")
      coefficients_distributions[, i, j] = as.vector(beta[-1])
    }
    
    # Average coefficient across all bootstrap samples
    coefficients_matrix[, i] = apply(coefficients_distributions[, i, ], 1, mean)
  }
  
  return(list(coefficients = coefficients_matrix, 
              coeff_dist = coefficients_distributions))
}
get_trinary = function(res,prob1=0.05,prob2=0.95,point_cut=0.2){
  trinary_matrix = matrix(0, nrow = nrow(res$coefficients), ncol = ncol(res$coefficients))
  for (i in  1:ncol(res$coefficients)){
    ci_lower = apply(res$coeff_dist[, i, ], 1, quantile, probs = prob1)
    ci_upper = apply(res$coeff_dist[, i, ], 1, quantile, probs = prob2)
    #point = res$coefficients[,i]
    trinary_matrix[, i] = ifelse(ci_lower > point_cut, 1, 
                                 ifelse(ci_upper < (-point_cut), -1, 0))
  }
  return(trinary_matrix)
}


## read all gene names
DAPs = readRDS("result/pp/luad_daps_limma.rds") #aceytal

## all histones, which are differentially acetylated?
genes = unique(DAPs$`gene symbol`)
H1 = genes[str_starts(genes,"H1")] ##linker
H2A = c(genes[str_starts(genes,"H2A")],genes[str_starts(genes,"MACROH2A")])
H2B = genes[str_starts(genes,"H2B")]
H3 = genes[str_starts(genes,"H3")]
H4 = genes[str_starts(genes,"H4")] ## not found in this dataset


histones = c(H1,H2A,H2B,H3,H4)
sum(DAPs$`gene symbol`%in%histones) #67
DAPs_histone = DAPs[DAPs$`gene symbol`%in%histones,]


## correlate histone acetylation with HDAC/HAT abundance/phosphorylation
## differential co-expression between HDACs/HATs and histones
## network plot or heatmap plot, clusterprofiler

pro = readRDS("data/luad_proteomics_imputed.rds")
ace = readRDS("data/luad_acetylproteomics_imputed.rds")
clinic = readRDS("data/luad_cptac_clinic.rds")
gender = clinic$sex

identical(colnames(pro),rownames(clinic))
identical(colnames(ace),rownames(clinic))

ace_histone = ace[str_detect(rownames(ace),paste0(histones,collapse = "|")),]
H1_K = rownames(ace)[str_detect(rownames(ace),paste0(H1,collapse = "|"))]
H2A_K = rownames(ace)[str_detect(rownames(ace),paste0(H2A,collapse = "|"))]
H2B_K = rownames(ace)[str_detect(rownames(ace),paste0(H2B,collapse = "|"))]
H3_K = rownames(ace)[str_detect(rownames(ace),paste0(H3,collapse = "|"))]

ace_histone = ace[c(H1_K,H2A_K,H2B_K,H3_K),]


### this list is curate from 
### Jonas, 
### WikiPedia, 
### http://www3.iiserpune.ac.in/~coee/histome/ptm.php?mod_type=lysine_acetylation, 
### and PTM paper supplement Table 5: https://www.cell.com/cell/fulltext/S0092-8674(23)00781-X?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS009286742300781X%3Fshowall%3Dtrue#supplementaryMaterial

HATs = c("CREBBP", "EP300","NCOA1","NCOA2","NCOA3", "ATF2","TAF1","HAT1", "KAT2A","KAT2B", "KAT5","KAT6A", "KAT6B", "KAT7", "KAT8")
HDACs = c("HDAC1", "HDAC2", "HDAC3", "HDAC4", "HDAC5", "HDAC6", "HDAC7", "HDAC8", "HDAC9", "HDAC10", "HDAC11", "SIRT1", "SIRT2", "SIRT3","SIRT6")

HATs = intersect(HATs,rownames(pro))
HDACs = intersect(HDACs,rownames(pro))
histone_regs = unique(c(HATs,HDACs))
pro_reg = pro[histone_regs,]

## use lasso
res_F = perform_lasso_with_bootstrap(X = pro_reg[,gender=="Female"],
                                     Y = ace_histone[,gender=="Female"],
                                     n_bootstrap = 100,
                                     seed = 42)
res_M = perform_lasso_with_bootstrap(X = pro_reg[,gender=="Male"],
                                     Y = ace_histone[,gender=="Male"],
                                     n_bootstrap = 100,
                                     seed = 42)

## get trinary matrix
res_F$trinary_matrix = get_trinary(res_F,
                                   prob1=0.05,
                                   prob2=0.95,
                                   point_cut=0.2)

#saveRDS(res_F,"result/pp/lasso_acetyl_F.rds")
#res_F = readRDS("result/pp/lasso_acetyl_F.rds")

mat_F = res_F$trinary_matrix #signed binary mat
rownames(mat_F) = rownames(pro_reg)
colnames(mat_F) = rownames(ace_histone)

#saveRDS(mat_F,"result/pp/net_acetyl_F.rds")
#write.csv(mat_F,"result/pp/net_acetyl_F.csv")

res_M$trinary_matrix = get_trinary(res_M,
                                   prob1=0.05,
                                   prob2=0.95,
                                   point_cut=0.2)

#saveRDS(res_M,"result/pp/lasso_acetyl_M.rds")
#res_M = readRDS("result/pp/lasso_acetyl_M.rds")

mat_M = res_M$trinary_matrix #signed binary mat
rownames(mat_M) = rownames(pro_reg)
colnames(mat_M) = rownames(ace_histone)

#saveRDS(mat_M,"result/pp/net_acetyl_M.rds")
#write.csv(mat_M,"result/pp/net_acetyl_M.csv")

#mat_diff = mat_F - mat_M
mat_diff = abs(mat_F) - abs(mat_M)
DeDi = rowSums(mat_diff) # degree change
sort(DeDi)


########### Plot arrow map ######
# read data
#female_matrix = readRDS("result/pp/net_acetyl_F.rds")
#male_matrix = readRDS("result/pp/net_acetyl_M.rds")
female_matrix = mat_F
male_matrix = mat_M


# filter rows and column
row_keep = rowSums(female_matrix!=0)>0 | rowSums(male_matrix!=0)>0
col_keep = colSums(female_matrix!=0)>0 | colSums(male_matrix!=0)>0

female_matrix = female_matrix[row_keep,col_keep]
male_matrix = male_matrix[row_keep,col_keep]

saveRDS(female_matrix,"result/pp/lasso_final_F.rds")
saveRDS(male_matrix,"result/pp/lasso_final_M.rds")


HATs_filt = intersect(HATs,rownames(male_matrix))
HDACs_filt = intersect(HDACs,rownames(male_matrix))
H1_K_filt = intersect(H1_K,colnames(male_matrix))
H2A_K_filt = intersect(H2A_K,colnames(male_matrix))
H2B_K_filt = intersect(H2B_K,colnames(male_matrix))
H3_K_filt = intersect(H3_K,colnames(male_matrix))

# group
row_groups = c(rep("HATs",length(HATs_filt)),rep("HDACs",length(HDACs_filt)))
col_groups = c(rep("H1",length(H1_K_filt)),rep("H2A",length(H2A_K_filt)),
               rep("H2B",length(H2B_K_filt)),rep("H3",length(H3_K_filt)))

# clustering on the row to find the row order
heat = pheatmap::pheatmap(male_matrix[HATs_filt,],
                          cluster_rows = T,
                          cluster_cols = F)
HATs_order = heat$tree_row$order

heat = pheatmap::pheatmap(male_matrix[HDACs_filt,],
                          cluster_rows = T,
                          cluster_cols = F)
HDACs_order = heat$tree_row$order

male_matrix = male_matrix[c(HATs_order,HDACs_order+length(HATs_order)),]
female_matrix = female_matrix[rownames(male_matrix),]

# Prepare data for plotting and add group information
data_graph <- expand.grid(x = 1:length(col_groups), y = 1:length(row_groups))
data_graph$row_group <- rep(row_groups, each = length(col_groups))
data_graph$col_group <- rep(col_groups, times = length(row_groups))
data_graph$female_value <- as.vector(t(female_matrix))
data_graph$male_value <- as.vector(t(male_matrix))
data_graph$xname <- rep(colnames(female_matrix), times = length(row_groups))
data_graph$yname <- rep(rownames(female_matrix), each = length(col_groups))
data_graph$`lasso coefficients` <- factor(with(data_graph, ifelse(female_value == 1 | male_value == 1, "Positive", 
                                                                  ifelse(female_value == -1 | male_value == -1, "Negative", "None"))),
                                          levels = c("Positive", "Negative", "None"))


data = data_graph
data$X = data$x
data$Y = data$y
data$Value1 = data$female_value
data$Value2 = data$male_value
data$row_group = as.factor(data$row_group)
data$col_group = as.factor(data$col_group)


# Set the spacing between arrows
arrow_spacing = 0.1
arrow_size = 1.2
arrow_ratio = 0.1  # Adjust this value to control the ratio between arrowhead and line

# Create the grid plot with parallel arrows
gg=ggplot(data, aes(x = xname, y = yname)) +
  geom_point(data = subset(data, Value1==0 & Value2==0),
             size = 0, aes(x=as.numeric(X),y=as.numeric(Y)))+
  geom_segment(data = subset(data, Value1 != 0),
               aes(x = as.numeric(X) - arrow_spacing, xend = as.numeric(X) - arrow_spacing, 
                   y =  as.numeric(Y) - Value1*0, yend = as.numeric(Y) + Value1*0.3, color = "Female"),
               arrow = arrow(length = unit(arrow_size, "inches") * arrow_ratio, 
                             ends = "last", 
                             angle = 20),
               size = arrow_size
  ) +
  geom_segment(data = subset(data, Value2 != 0),
               aes(x = as.numeric(X) + arrow_spacing, xend = as.numeric(X) + arrow_spacing, 
                   y = as.numeric(Y) - Value2*0, yend = as.numeric(Y) + Value2*0.3, color = "Male"),
               arrow = arrow(length = unit(arrow_size, "inches") * arrow_ratio, 
                             ends = "last", 
                             angle = 20),
               size = arrow_size
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "#D3D3D3"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16,face = "bold"),
    strip.background = element_rect(fill = "#D3D3D3"),
    strip.text = element_text(size = 16, face = "bold")
  ) +
  labs(color = "") +
  facet_grid(rows = vars(row_group), scales = "free", space = "fixed", switch = "y",drop = T)+
  scale_y_continuous(breaks = unique(data$Y), labels = unique(data$yname))+
  scale_x_continuous(breaks = unique(data$X), labels = unique(data$xname))+
  scale_color_manual(values = c("Female" = scales::alpha("#F8766D", 1), 
                                "Male" = scales::alpha("#619CFF", 1)))

gg
ggsave("plot/heatmap/acetylation_arrow_map_clustered.png",gg,height = 4,width = 7)

# ###### plot legend ######
# library(scales)
# 
# # Create arrow legend plot
# legend_data <- expand.grid(
#   color = c("#F8766D","#619CFF"),
#   direction = c("up", "down"),
#   stringsAsFactors = FALSE
# )
# 
# # Adjust positions for 1x4 layout
# legend_data$x <- c(1,5,9,13)  # 4 positions along the x-axis
# legend_data$y <- rep(0, 4)           # Same y position for all
# 
# # Define text and arrow directions
# legend_data$text <- rep(c("Female", "Male"), times = 2)
# legend_data$dir_y <- -0.8
# 
# # Create the custom legend
# legend_plot = ggplot(legend_data, aes(x = x, y = y, color = color)) +
#   geom_segment(aes(xend = x, 
#                    y = ifelse(direction == "up", y, y - 0.3),
#                    yend = ifelse(direction == "up", y + 0.3, y - 0.6)),
#                arrow = arrow(length = unit(1, "inches") * 0.8, 
#                              #type = "closed", 
#                              ends = "last", 
#                              angle = 20),
#                size = 4) +
#   geom_text(aes(x = x, y = y + 0.8, label = text), vjust = 0,size=14) +
#   geom_text(aes(x = x, y = y - 1.4, label = direction), vjust = 0,size=14) +
#   scale_color_identity() +
#   theme_void() +
#   theme(legend.position = "none") +
#   xlim(-2, 16) + ylim(-2.2, 2.2) 
# # coord_fixed(ratio = 1)
# 
# legend_plot
# ggsave("plot/heatmap/histone_reg_legend.png",legend_plot,width = 8,height = 3.5)
# 
# 
# 
