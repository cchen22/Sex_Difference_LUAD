library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(BioEnricher)

## data
DEPs = readRDS("result/pp/luad_deps_removeY.rds") #pro
DPPs2 = readRDS("result/pp/luad_dpps2_removeY.rds") #phos
DAPs2 = readRDS("result/pp/luad_daps2_removeY.rds") #ace
DEGs = readRDS("result/pp/luad_degs_removeY.rds")



##### GSEA DEPs ######
res = DEPs
geneList = res$t
names(geneList) = res$`gene symbol`
geneList = sort(geneList,decreasing = T)

## gse kegg
gene_id = bitr(names(geneList), 
               fromType="SYMBOL", 
               toType = "ENTREZID",
               OrgDb = org.Hs.eg.db) 
names(geneList) = gene_id$ENTREZID
kk = gseKEGG(geneList     = geneList,
             keyType = "kegg",
             organism     = 'hsa',
             minGSSize    = 10,
             maxGSSize    = 500,
             pvalueCutoff = 0.05,
             verbose      = FALSE,
             seed = 123)
kk = setReadable(kk,OrgDb = org.Hs.eg.db, keyType="ENTREZID")
saveRDS(kk,"result/pp/DEP_kegg_gsea.rds")

kk = readRDS("result/pp/DEP_kegg_gsea.rds")
write.csv(kk@result,"result/pp/DEP_kegg_gsea.csv",row.names = F)

sum(kk@result$p.adjust<0.05)
sum(kk@result$p.adjust<0.025)

## dotplot
gg = dotplot(kk,showCategory=sum(kk@result$p.adjust<0.05),
             color = "p.adjust",
             x="NES",
             size = NULL,
             title="",
             font.size = 20,
             label_format = 60)
ggsave("plot/dotplot/GSEA_DEPs_all.png",gg,width = 10,height = 12)


#### ORA DPPs #####
DPPs = DPPs2 %>%
  dplyr::group_by(`gene symbol`) %>%
  dplyr::filter(P.Value == min(P.Value, na.rm = TRUE)) %>%
  dplyr::ungroup()

phos_luad = DPPs[DPPs$P.Value<0.05, "gene symbol"]%>%pull(`gene symbol`)
phos_luad_F = DPPs[DPPs$P.Value<0.05 & DPPs$logFC>0,"gene symbol"]%>%pull(`gene symbol`)
phos_luad_M = DPPs[DPPs$P.Value<0.05 & DPPs$logFC<0,"gene symbol"]%>%pull(`gene symbol`)

## male 
gene_id = bitr(phos_luad_M,
               fromType="SYMBOL",
               toType = "ENTREZID",
               OrgDb = org.Hs.eg.db)

kk <- enrichKEGG(gene         = gene_id$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
saveRDS(kk,"result/pp/DPP_kegg_ora_M.rds")

kk = readRDS("result/pp/DPP_kegg_ora_M.rds")
write.csv(kk@result[kk@result$p.adjust<0.05,],
          "result/pp/DPP_kegg_ora_M.csv",row.names = F)


sum(kk@result$p.adjust<0.05)
sum(kk@result$p.adjust<0.25)

gg = dotplot(kk,showCategory=sum(kk@result$p.adjust<0.05),size = "Count",
             title="Male",font.size = 20,label_format = 60)
ggsave("plot/dotplot/ORA_DPPs2_M.png",gg,width = 10,height = 6)


## female
gene_id = bitr(phos_luad_F,
               fromType="SYMBOL",
               toType = "ENTREZID",
               OrgDb = org.Hs.eg.db)

kk <- enrichKEGG(gene         = gene_id$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
saveRDS(kk,"result/pp/DPP_kegg_ora_F.rds")

kk = readRDS("result/pp/DPP_kegg_ora_F.rds")
write.csv(kk@result[kk@result$p.adjust<0.05,],
          "result/pp/DPP_kegg_ora_F.csv",row.names = F)

sum(kk@result$p.adjust<0.05)
sum(kk@result$p.adjust<0.25)

gg = dotplot(kk,showCategory=sum(kk@result$p.adjust<0.05),size = "Count",
             title="Female",font.size = 20,label_format = 60)
ggsave("plot/dotplot/ORA_DPPs2_F.png",gg,width = 10,height = 6)

#### ORA DAPs #####
DAPs = DAPs2 %>%
  dplyr::group_by(`gene symbol`) %>%
  dplyr::filter(P.Value == min(P.Value, na.rm = TRUE)) %>%
  dplyr::ungroup()

acey_luad = DAPs[DAPs$P.Value<0.05,"gene symbol"]%>%pull(`gene symbol`)
acey_luad_F = DAPs[DAPs$P.Value<0.05 & DAPs$logFC>0,"gene symbol"]%>%pull(`gene symbol`)
acey_luad_M = DAPs[DAPs$P.Value<0.05 & DAPs$logFC<0,"gene symbol"]%>%pull(`gene symbol`)

## male results good
gene_id = bitr(acey_luad_M,
               fromType="SYMBOL",
               toType = "ENTREZID",
               OrgDb = org.Hs.eg.db)

kk <- enrichKEGG(gene         = gene_id$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
saveRDS(kk,"result/pp/DAP_kegg_ora_M.rds")

kk = readRDS("result/pp/DAP_kegg_ora_M.rds")
write.csv(kk@result[kk@result$p.adjust<0.05,],
          "result/pp/DAP_kegg_ora_M.csv",row.names = F)


sum(kk@result$p.adjust<0.05)
sum(kk@result$p.adjust<0.25)

gg = dotplot(kk,showCategory=sum(kk@result$p.adjust<0.05),size = "Count",
             title="Male",font.size = 20,label_format = 60)
ggsave("plot/dotplot/ORA_DAPs2_M.png",gg,width = 10,height = 8)

## female results need to cut at 0.25
gene_id = bitr(acey_luad_F,
               fromType="SYMBOL",
               toType = "ENTREZID",
               OrgDb = org.Hs.eg.db)

kk <- enrichKEGG(gene         = gene_id$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.25)
saveRDS(kk,"result/pp/DAP_kegg_ora_F.rds")

kk = readRDS("result/pp/DAP_kegg_ora_F.rds")
write.csv(kk@result[kk@result$p.adjust<0.05,],
          "result/pp/DAP_kegg_ora_F.csv",row.names = F)

sum(kk@result$p.adjust<0.05)
sum(kk@result$p.adjust<0.25)

gg = dotplot(kk,showCategory=sum(kk@result$p.adjust<0.05),size = "Count",
             title="Female",font.size = 20,label_format = 60)
ggsave("plot/dotplot/ORA_DAPs2_F.png",gg,width = 10,height = 8)


##### GSEA DEGs ######
res = DEGs
geneList = res$t
names(geneList) = res$`gene_symbol`
geneList = sort(geneList,decreasing = T)

## gse kegg
gene_id = bitr(names(geneList), 
               fromType="SYMBOL", 
               toType = "ENTREZID",
               OrgDb = org.Hs.eg.db) 
names(geneList) = gene_id$ENTREZID
kk = gseKEGG(geneList     = geneList,
             keyType = "kegg",
             organism     = 'hsa',
             minGSSize    = 10,
             maxGSSize    = 500,
             pvalueCutoff = 0.05,
             verbose      = FALSE,
             seed = 123)
kk = setReadable(kk,OrgDb = org.Hs.eg.db, keyType="ENTREZID")
saveRDS(kk,"result/pp/DEG_kegg_gsea.rds")

kk = readRDS("result/pp/DEG_kegg_gsea.rds")
write.csv(kk@result[kk@result$p.adjust<0.05,],
          "result/pp/DEG_kegg_gsea.csv",row.names = F)

sum(kk@result$p.adjust<0.05)
sum(kk@result$p.adjust<0.0001)

## dotplot
gg = dotplot(kk,showCategory=20,
             color = "p.adjust",
             x="NES",
             size = NULL,
             title="",
             font.size = 20,
             label_format = 60)
ggsave("plot/dotplot/GSEA_DEGs_top20.png",gg,width = 10,height = 12)

