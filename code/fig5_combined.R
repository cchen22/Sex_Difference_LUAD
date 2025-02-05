library(tidyverse) 
library(purrr)
library(cowplot)  

## heatmap
DEGs = readRDS("result/pp/luad_degs.rds") %>%
  mutate("gene symbol" = gene_symbol,"score" = sign(logFC)*ifelse(adj.P.Val<0.1,1,0)) %>%
  select("gene symbol","score")

DEPs = readRDS("result/pp/luad_deps_limma.rds") %>%
  mutate("score" = sign(logFC)*ifelse(adj.P.Val<0.1,1,0)) %>%
  select("gene symbol","score")

DPPs = readRDS("result/pp/luad_dpps_limma.rds") %>%
  mutate("score" = sign(logFC)*ifelse(adj.P.Val<0.1,1,0)) %>%
  select("gene symbol","score")

DAPs = readRDS("result/pp/luad_daps_limma.rds") %>%
  mutate("score" = sign(logFC)*ifelse(adj.P.Val<0.1,1,0)) %>%
  select("gene symbol","score")

DATs = readRDS("result/pp/luad_dats_limma.rds") %>%
  mutate("score" = sign(logFC)*ifelse(adj.P.Val<0.1,1,0)) %>%
  select("gene symbol","score")

#DAKs = readRDS("result/pp/luad_daks_limma.rds") %>%
DAKs = readRDS("result/pp/luad_dPTMgsea_limma.rds") %>%
  mutate("score" = sign(logFC)*ifelse(adj.P.Val<0.1,1,0)) %>%
  select("gene symbol","score")
DAKs$`gene symbol`[DAKs$`gene symbol`=="P38A/MAPK14"] = "MAPK14"
DAKs$`gene symbol`[DAKs$`gene symbol`=="AurA/AURKA"] = "AURKA"


#list_of_dfs = list(DEGs, DEPs, DATs, DAKs)
list_of_dfs = list(DEGs, DEPs, DPPs, DAPs, DATs, DAKs)
df = purrr::reduce(list_of_dfs, full_join, by = "gene symbol") %>%
  mutate(across(everything(), ~replace_na(., 0)))
#colnames(df) = c("symbol","Gene","Protein","TF","Kinase")
colnames(df) = c("symbol","Gene","Protein","Phosphorylation","Acetylation","TF","Kinase")


## target genes
gene_name = c("AURKA","MAPK14","NR3C1","AR","ESR1",
              "HDAC1","HDAC6","HDAC8","HDAC10","MAPK1")
df = df[df$symbol %in% gene_name,]
df = df %>%
  group_by(symbol) %>%
  mutate(non_zero_count = sum(select(., -symbol) != 0)) %>%
  filter(non_zero_count == max(non_zero_count)) %>%
  slice(1) %>%
  ungroup() %>%
  select(-non_zero_count)

df$"Network"  = ifelse(df$symbol=="AR",-1,
                       ifelse(df$symbol=="ESR1",1,
                              ifelse(df$symbol=="MAPK1",-1,0)))
df$"HDAC" = ifelse(str_detect(df$symbol,"HDAC"),-1,0)

# merge HDAC rows
df = df %>%
  mutate(symbol = if_else(str_detect(symbol,"HDAC"), "HDAC", symbol)) %>%
  group_by(symbol) %>%
  summarise(across(where(is.numeric), ~if(all(is.na(.))) NA else max(., na.rm = TRUE)))

## pheatmap
df_new = column_to_rownames(df,var = "symbol")
df_new  =  as.matrix(df_new)
df_new  = df_new[c("AURKA","MAPK14","MAPK1","NR3C1","AR","ESR1","HDAC"),]
heat = pheatmap::pheatmap(
  df_new[,-(1:4)],
  scale = 'none',
  color = colorRampPalette(c("#4472C4", "white", "#ec5e56"))(255),
  cluster_rows = F,
  treeheight_row = 0,
  cluster_cols = F,
  show_colnames = T,
  show_rownames = F,
  border_color = 'grey50',
  angle_col = "45",
  legend = F,
  fontsize = 12,
  fontsize_col = 12,
  cellwidth = 16,
  cellheight = 16
)

ggsave("plot/heatmap/drug_gene_heat.png",heat,width = 1,height = 2.5)



## line plot
drug_gene_df = readRDS("result/pp/drug_gene_df.rds")
drug_gene_df = drug_gene_df[drug_gene_df$pvalue<0.05,]
drug_gene_df$gene = factor(drug_gene_df$gene,levels = c("AURKA","MAPK14","MAPK1","NR3C1","AR","ESR1","HDAC"))
drug_gene_df = drug_gene_df[order(drug_gene_df$gene),]
drug_gene_df$gene = as.character(drug_gene_df$gene)

drug_gene_df$"len.drug" = as.integer(factor(drug_gene_df$drug,levels = unique(drug_gene_df$drug)))
tmp = unique(sort(drug_gene_df$len.drug, decreasing = TRUE))
drug_gene_df$len.drug = match(drug_gene_df$len.drug, tmp)

drug_gene_df$"len.gene" = as.integer(factor(drug_gene_df$gene,levels = unique(drug_gene_df$gene)))
tmp = unique(sort(drug_gene_df$len.gene, decreasing = TRUE))
drug_gene_df$len.gene = match(drug_gene_df$len.gene, tmp)+9

drug_gene_df$"group" = 1:nrow(drug_gene_df)

df = data.frame("label" = c(drug_gene_df$drug,drug_gene_df$gene),
                "type" = c(rep("Drug",length(drug_gene_df$drug)),rep("Gene",length(drug_gene_df$gene))),
                "len" = c(drug_gene_df$len.drug,drug_gene_df$len.gene),
                "group" = c(drug_gene_df$group,drug_gene_df$group))



line = df %>% 
  ggplot(aes(x = type, y = len,group=group)) +  
  geom_line() +  
  geom_point(size = 2) +  
  geom_text(aes(label = label),
            size=4,
            nudge_x = ifelse(df$type == "Drug", -0.02, 0.02),
            hjust = ifelse(df$type == "Drug", 1, 0)) +
  scale_x_discrete(limits = c("Drug", "Gene"),
                   expand = expansion(add = c(0.5, 0.2),
                                      mult = c(1.2, 0.2)),
                   position = "top") +  
  scale_color_brewer(palette = "Paired") +  
  theme(
    plot.background = element_blank(),  
    panel.background = element_blank(),
    legend.position = "none",  
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  
    axis.text.y = element_blank(),  
    axis.text.x = element_blank(),
    axis.title = element_blank(), 
    axis.ticks = element_blank()
  )

ggsave("plot/heatmap/lineplot.png",line,width = 5.8,height = 4)



