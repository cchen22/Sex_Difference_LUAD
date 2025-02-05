library(pheatmap)
library(tidyverse)
library(viridis)
library(scales)
heatmap_plot = function(df,res,pcutoff,gender,title){
  MOI_F = rownames(res[res$adj.P.Val<pcutoff & res$logFC>0,])%>%unique()
  MOI_M = rownames(res[res$adj.P.Val<pcutoff & res$logFC<0,])%>%unique()
  df = df[c(MOI_F,MOI_M),]
  
  df_M = rowMeans(df[,gender=="male"])
  df_F = rowMeans(df[,gender=="female"])
  df_final = cbind(df_F,df_M)
  rownames(df_final)=rownames(df)
  colnames(df_final)=c("F","M")
  
  pheatmap::pheatmap(
    df_final,
    scale = 'none',
    #color = magma(25, direction = 1),
    cluster_rows = F,
    cluster_cols = F,
    show_colnames = T,
    show_rownames = T,
    border_color = 'grey50',
    angle_col = 0,
    legend = T,
    main = title,
    fontsize = 12,
    cellwidth = 12,
    cellheight = 12
  )
  
}

pcutoff = 0.05

### data: diff test results
DATs = readRDS("result/pp/luad_dats_removeY.rds") #TF
DEPs = readRDS("result/pp/luad_deps_removeY.rds") #pro
DPPs = readRDS("result/pp/luad_dpps2_removeY.rds") #phos
DAKs = readRDS("result/pp/luad_daks_removeY.rds") #kinase
DEGs = readRDS("result/pp/luad_degs_removeY.rds")
DAPs = readRDS("result/pp/luad_daps2_removeY.rds")
DPTMgsea = readRDS("result/pp/luad_dPTMgsea_removeY.rds")

## data2: TFA, proteome, phosphorylation, KA
tiger = readRDS("result/pp/luad_tiger.rds")
TFA = tiger$Z #721*528
expr = readRDS("data/pp/luad_expression.rds")
pro = readRDS("data/luad_proteomics_imputed.rds")
phos = readRDS("data/luad_phosphoproteomics_imputed.rds")
ace = readRDS("data/luad_acetylproteomics_imputed.rds")
KA = readRDS("result/pp/phos_activity_tiger.rds")

PTMgsea = cmapR::parse.gctx("data/PTM-SEA/luad_phos-scores.gct")
PTMgsea = PTMgsea@mat

## clinical data of cptac
clinic = readRDS("data/luad_cptac_clinic.rds")
gender = clinic$sex
gender[gender=="Female"] = "female"
gender[gender=="Male"] = "male"
gender_CPTAC = factor(gender,levels = c("female","male"))

## clincal data of tcga
rse = readRDS("data/pp/luad_rse.rds")
identical(colnames(TFA),colnames(rse))
gender = rse@colData$tcga.gdc_cases.demographic.gender
gender_TCGA = factor(gender,levels = c("female","male"))


#################
#### heatmap####

## TFA
df = t(scale(t(TFA)))
res = DATs
heat = heatmap_plot(df,res,pcutoff=0.05,gender = gender_TCGA,title = "TFA")
ggsave("plot/heatmap/luad_TFA.png",heat,width = 2,height = 4.5)


## kinase activity from TIGER
df = t(scale(t(KA)))
res = DAKs
heat = heatmap_plot(df,res,pcutoff=0.25,gender = gender_CPTAC,title = "KA1")
ggsave("plot/heatmap/luad_KA.png",heat,width = 2,height = 1.5)

## pathway activity of PTM GSEA
pathway_act = PTMgsea[c(grep("PATH",rownames(PTMgsea)),grep("KINASE",rownames(PTMgsea))),]
df = t(scale(t(pathway_act)))
res = DPTMgsea[DPTMgsea$class=="PATH"|DPTMgsea$class=="KINASE",]
heat = heatmap_plot(df,res,pcutoff=0.1,gender = gender_CPTAC,title = "KA2")
ggsave("plot/heatmap/luad_PTMpathway.png",heat,width = 3.5,height = 3)


