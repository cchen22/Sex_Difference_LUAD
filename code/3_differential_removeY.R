library(biomaRt)

## read in DIFF results

## TCGA + CPTAC
DATs = readRDS("result/pp/luad_dats_limma.rds") #TF
DEPs = readRDS("result/pp/luad_deps_limma.rds") #pro
DPPs = readRDS("result/pp/luad_dpps_limma.rds") #phos
DPPs2 = readRDS("result/pp/luad_dpps2_limma.rds") #phos
DAPs = readRDS("result/pp/luad_daps_limma.rds") #ace
DAPs2 = readRDS("result/pp/luad_daps2_limma.rds") #ace
DAKs = readRDS("result/pp/luad_daks_limma.rds") #kinase
DEGs = readRDS("result/pp/luad_degs.rds")
DPTMgsea = readRDS("result/pp/luad_dPTMgsea_limma.rds")


## post-process: remove Y chromosome genes
gene_symbol = unique(c(DATs$`gene symbol`,
                       DEPs$`gene symbol`,
                       DPPs$`gene symbol`,
                       DAPs$`gene symbol`,
                       DAKs$`gene symbol`,
                       DEGs$gene_symbol,
                       DPTMgsea$`gene symbol`))
ensembldb = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_info = getBM(attributes = c("external_gene_name","chromosome_name"), 
                  filters = "external_gene_name", 
                  values = gene_symbol, 
                  mart = ensembldb)
Y_gene = gene_info[gene_info$chromosome_name%in%c("Y"),1]

## remove Y_genes 
DATs = DATs[!DATs$`gene symbol`%in%Y_gene,]
DEPs = DEPs[!DEPs$`gene symbol`%in%Y_gene,]
DPPs = DPPs[!DPPs$`gene symbol`%in%Y_gene,]
DPPs2 = DPPs2[!DPPs2$`gene symbol`%in%Y_gene,]
DAPs = DAPs[!DAPs$`gene symbol`%in%Y_gene,]
DAPs2 = DAPs2[!DAPs2$`gene symbol`%in%Y_gene,]
DAKs = DAKs[!DAKs$`gene symbol`%in%Y_gene,]
DPTMgsea = DPTMgsea[!DPTMgsea$`gene symbol`%in%Y_gene,]

DEGs = DEGs[!DEGs$`gene_symbol`%in%Y_gene,]


## save
saveRDS(DATs,"result/pp/luad_dats_removeY.rds")
saveRDS(DEPs,"result/pp/luad_deps_removeY.rds")
saveRDS(DPPs,"result/pp/luad_dpps_removeY.rds")
saveRDS(DPPs2,"result/pp/luad_dpps2_removeY.rds")
saveRDS(DAPs,"result/pp/luad_daps_removeY.rds")
saveRDS(DAPs2,"result/pp/luad_daps2_removeY.rds")
saveRDS(DAKs,"result/pp/luad_daks_removeY.rds")
saveRDS(DPTMgsea,"result/pp/luad_dPTMgsea_removeY.rds")

saveRDS(DEGs,"result/pp/luad_degs_removeY.rds")




