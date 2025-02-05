
## read in gsea results
deg_gsea = readRDS("result/pp/DEG_kegg_gsea.rds")
dep_gsea = readRDS("result/pp/DEP_kegg_gsea.rds")

dpp_ora_F = readRDS("result/pp/DPP_kegg_ora_F.rds")
dpp_ora_M = readRDS("result/pp/DPP_kegg_ora_M.rds")

dap_ora_F = readRDS("result/pp/DAP_kegg_ora_F.rds")
dap_ora_M = readRDS("result/pp/DAP_kegg_ora_M.rds")


deg = deg_gsea@result
deg_c = deg[deg$p.adjust<0.05,]$Description

dep = dep_gsea@result
dep_c = dep[dep$p.adjust<0.05,]$Description

dpp = unique(rbind(dpp_ora_F@result,dpp_ora_M@result))
dpp_c = dpp[dpp$p.adjust<0.05,]$Description

dap = unique(rbind(dap_ora_F@result,dap_ora_M@result))
dap_c = dap[dap$p.adjust<0.05,]$Description


# pairwise common kegg pathway
vectors <- list(dep_c, deg_c, dap_c, dpp_c)
pairwise_intersections <- combn(vectors, 2, function(x) intersect(x[[1]], x[[2]]), simplify = FALSE)
names(pairwise_intersections) <- combn(length(vectors), 2, function(x) paste("Intersection of Vector", x[1], "and", x[2]), simplify = FALSE)
common_kegg = unique(unlist(pairwise_intersections)) 

# heatmap
deg_M = deg%>%
  dplyr::filter(Description%in%common_kegg) %>%
  dplyr::filter(NES<0) %>%
  mutate(data_type = "mRNA",score=ifelse(p.adjust<0.05,-1,0)) %>%
  dplyr::select(Description,data_type,score)

dep_M = dep%>%
  dplyr::filter(Description%in%common_kegg) %>%
  dplyr::filter(NES<0) %>%
  mutate(data_type = "protein",score=ifelse(p.adjust<0.05,-1,0)) %>%
  dplyr::select(Description,data_type,score)

dpp_M = dpp_ora_M@result%>%
  dplyr::filter(Description%in%common_kegg) %>%
  mutate(data_type = "phosphorylation",score=ifelse(p.adjust<0.05,-1,0)) %>%
  dplyr::select(Description,data_type,score)

dap_M = dap_ora_M@result%>%
  dplyr::filter(Description%in%common_kegg) %>%
  mutate(data_type = "acetylation",score=ifelse(p.adjust<0.05,-1,0)) %>%
  dplyr::select(Description,data_type,score)

df_M = bind_rows(deg_M, dep_M, dpp_M, dap_M) %>%
  spread(key = data_type, value = score) %>%
  column_to_rownames("Description")
df_M[is.na(df_M)] <- 0
df_M$"mRNA" = 0
df_M = as.matrix(df_M)
df_M = df_M[,c('mRNA',"protein","phosphorylation","acetylation")]


deg_F = deg%>%
  dplyr::filter(Description%in%common_kegg) %>%
  dplyr::filter(NES>0) %>%
  mutate(data_type = "mRNA",score=ifelse(p.adjust<0.05,1,0)) %>%
  dplyr::select(Description,data_type,score)

dep_F = dep%>%
  dplyr::filter(Description%in%common_kegg) %>%
  dplyr::filter(NES>0) %>%
  mutate(data_type = "protein",score=ifelse(p.adjust<0.05,1,0)) %>%
  dplyr::select(Description,data_type,score)

dpp_F = dpp_ora_F@result%>%
  dplyr::filter(Description%in%common_kegg) %>%
  mutate(data_type = "phosphorylation",score=ifelse(p.adjust<0.05,1,0)) %>%
  dplyr::select(Description,data_type,score)

dap_F = dap_ora_F@result%>%
  dplyr::filter(Description%in%common_kegg) %>%
  mutate(data_type = "acetylation",score=ifelse(p.adjust<0.05,1,0)) %>%
  dplyr::select(Description,data_type,score)

df_F = bind_rows(deg_F, dep_F, dpp_F, dap_F) %>%
  spread(key = data_type, value = score) %>%
  column_to_rownames("Description")
df_F[is.na(df_F)] <- 0
df_F = as.matrix(df_F)
df_F = df_F[,c('mRNA',"protein","phosphorylation","acetylation")]

## merge
df1 <- df_F %>% as.data.frame() %>% rownames_to_column(var = "rowname")
df2 <- df_M %>% as.data.frame() %>% rownames_to_column(var = "rowname")
df_final <- full_join(df1, df2, by = "rowname") %>% column_to_rownames("rowname")
df_final[is.na(df_final)] = 0
df_final = as.matrix(df_final)
df_final = df_final[common_kegg,]

# Plotting the heatmap
heat = pheatmap::pheatmap(
  df_final,
  scale = 'none',
  color = colorRampPalette(c("blue", "white", "red"))(255),
  cluster_rows = T,
  treeheight_row = 0,
  cluster_cols = F,
  show_colnames = T,
  show_rownames = T,
  labels_col = rep(c('mRNA',"protein","phosphorylation","acetylation"),2),
  border_color = 'grey50',
  angle_col = 45,
  legend = F,
  gaps_col = 4,
  gaps_row = c(3,5,10,33,39),
  fontsize = 12,
  fontsize_col = 12,
  cellwidth = 13,
  cellheight = 13
)

row_order = heat$tree_row$order
row_order = row_order[-c(1:3, 6:10,40:45)]


heat = pheatmap::pheatmap(
  df_final[row_order,],
  scale = 'none',
  color = colorRampPalette(c("#4472C4", "white", "#ec5e56"))(255),
  cluster_rows = F,
  treeheight_row = 0,
  cluster_cols = F,
  show_colnames = T,
  show_rownames = T,
  labels_col = rep(c('mRNA',"protein","phosphorylation","acetylation"),2),
  border_color = 'grey50',
  angle_col = "45",
  legend = F,
  gaps_col = 4,
  #gaps_row = c(3,5,10,33,39),
  gaps_row = c(2,25),
  fontsize = 12,
  fontsize_col = 12,
  cellwidth = 13,
  cellheight = 13
)

ggsave("plot/heatmap/luad_landscape_clean.png",heat,width = 7,height = 7)


