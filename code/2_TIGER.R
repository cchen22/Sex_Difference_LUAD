#library(TIGER)
library(netZooR)
library(dorothea)
library(OmnipathR)
library(tidyverse)

#### 0. dorothea prior remove E
db=dorothea::entire_database
db = db %>%  ## ABCD + doubleE
  dplyr::filter((confidence!="E") | (is_evidence_tfbs==T & is_evidence_inferred == T))
TFs = unique(db$tf) #753 TFs
df1 = dorothea::dorothea_hs_pancancer # TCGA
TCGA_prior = TIGER::el2adj(df1[df1$tf%in%TFs,-2]) #735 * 20419
saveRDS(TCGA_prior,"data/pp/TCGA_prior.rds")


#### 1. TCGA TFA TIGER
prior = readRDS("data/pp/TCGA_prior.rds")
expr = readRDS("data/pp/luad_expression.rds")
res = TIGER(expr,prior,TFexpressed = F)
saveRDS(res,"result/pp/luad_tiger.rds")

#### 2. GSE68465 TFA TIGER
prior = readRDS("data/pp/TCGA_prior.rds")
expr = readRDS("data/GSE68465/expr.rds")
res = TIGER(expr,prior)
saveRDS(res,"result/pp/GSE68465_tiger.rds")


#### 3. CPTAC KA TIGER
enz_sub = import_omnipath_enzsub()
enz_sub
phos_net = enz_sub %>% 
  filter(modification=="phosphorylation"|modification=="dephosphorylation") %>%
  mutate(target = paste0(substrate_genesymbol, "_", residue_type, residue_offset)) %>%
  mutate(source = enzyme_genesymbol)%>%
  mutate(mor = ifelse(modification == "phosphorylation", 1, -1))%>%
  group_by(source, target) %>% # some conflict edge sign
  mutate(mor = ifelse(any(mor == 1 & mor == -1), 1, mor))%>%
  select(source,target,mor)%>%
  unique()
phos_adj = TIGER::el2adj(phos_net)

clinic = readRDS("data/luad_cptac_clinic.rds")
phos = readRDS("data/luad_phosphoproteomics_imputed.rds")
phos = phos[,rownames(clinic)] # remove normal samples ".N"

phosite = intersect(colnames(phos_adj),rownames(phos))
keep = rowSums(phos_adj[,phosite]!=0)>=3 #at least 3 edges
kinase = rownames(phos_adj[keep,])

phos_tiger = TIGER(phos[phosite,],phos_adj[kinase,phosite],
                   TFexpressed = F,signed = T)
phos_activity = phos_tiger$Z
saveRDS(phos_activity,"result/pp/phos_activity_tiger.rds")


#### 4. APOLLO KA TIGER
enz_sub = import_omnipath_enzsub()
enz_sub
phos_net = enz_sub %>% 
  filter(modification=="phosphorylation"|modification=="dephosphorylation") %>%
  mutate(target = paste0(substrate_genesymbol, "_", residue_type, residue_offset)) %>%
  mutate(source = enzyme_genesymbol)%>%
  mutate(mor = ifelse(modification == "phosphorylation", 1, -1))%>%
  group_by(source, target) %>% # some conflict edge sign
  mutate(mor = ifelse(any(mor == 1 & mor == -1), 1, mor))%>%
  select(source,target,mor)%>%
  unique()
phos_adj = TIGER::el2adj(phos_net)

phos = readRDS("data/pp/apollo_phos_imputed_V2.rds")
phosite = intersect(colnames(phos_adj),rownames(phos))
keep = rowSums(phos_adj[,phosite]!=0)>=3 #at least 3 edges
kinase = rownames(phos_adj[keep,])

phos_tiger = TIGER(phos[phosite,],phos_adj[kinase,phosite],
                   TFexpressed = F,signed = T)
phos_activity = phos_tiger$Z
saveRDS(phos_activity,"result/pp/apollo_KA_tiger.rds")














