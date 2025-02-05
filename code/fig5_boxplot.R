library(tidyverse)
library(ggpubr)
library(scales)

##1. drug response dataset, filter by LUAD and MAPK
cdr = read_csv("data/CellLineDrug/Repurposing_Public_23Q2_Extended_Primary_Data_Matrix.csv")
cdr = rename_with(cdr,~ "drug ID", 1)

##2. drug informatio
drug = read_csv("data/CellLineDrug/Repurposing_Public_23Q2_Extended_Primary_Compound_List.csv")

##3. cancer cell lines, LUAD
ccl = read_csv("data/CellLineDrug/Model.csv")
keep = ccl$DepmapModelType%in%"LUAD"
ccl = ccl[keep,]

##4. add smoking information from model_list_20230902.csv data
model_list = read_csv("data/CellLineDrug/model_list_20230902.csv")
ccl = left_join(ccl,model_list,by=c("SangerModelID"="model_id"))
table(ccl$smoking_status,ccl$gender)

# Female Male
# Smoker      11   10
# Unknown     18   29

##4. double loops through target genes and drugs
gene_list = c("AURKA","MAPK14","NR3C1","AR","ESR1","HDAC",
              "MAPK1",
              "PPARG",
              "IKBKB",
              "RELA",
              "GSK3B",
              "NFKB1")
drug_MOA = c("AURORA KINASE",
             "MAPK|MAP KINASE",
             "GLUCOCORTICOID",
             "ANDROGEN",
             "ESTROGEN|PROGESTERONE",
             "HDAC",
             "FGFR INHIBITOR, KIT INHIBITOR, PDGFR TYROSINE KINASE RECEPTOR INHIBITOR, RAF INHIBITOR, RET",
             "PPAR",
             "IKK|NFKB",
             "NFKB|HIV",
             "GLYCOGEN SYNTHASE KINASE",
             ".+")

drug_gene_df = NULL
for (j in 1:length(gene_list)){
  gene_name = gene_list[j]
  MOA_name = drug_MOA[j]
  drug_ex = drug %>% 
    dplyr::filter(!is.na(repurposing_target)) %>%
    dplyr::filter(str_detect(repurposing_target,gene_name)) %>%
    dplyr::filter(str_detect(MOA,MOA_name)) 
  
  cdr_ex = cdr %>% 
    pivot_longer(!`drug ID`) %>%
    inner_join(drug_ex,by=c("drug ID"="IDs")) %>%
    inner_join(ccl,by=c('name'='ModelID'))
  
  
  ## inner loop on drugs
  drug_list = unique(cdr_ex$`drug ID`)
  drug_name = plyr::mapvalues(drug_list,cdr_ex$`drug ID`,cdr_ex$Drug.Name,warn_missing = F)
  
  for (i in 1:length(drug_list)){
    d = drug_list[i]
    y = cdr_ex[cdr_ex$`drug ID`==d,]
    y = y[y$Sex!="Unknown",]
    y_F = y[y$Sex=="Female",]
    y_M = y[y$Sex=="Male",]
    
    moa = drug_ex[drug_ex$IDs==d,"MOA"]
    pattern = paste(c("AGONIST", 
                      "ANTAGONIST", 
                      "MODULATOR",
                      "INHIBITOR",
                      "SERM",
                      "DESTABILIZER"),
                    collapse = "|")
    moa = stringr::str_extract(moa, pattern)
    
    ## filter if not enough cell lines
    if (sum(!is.na(y_F$value))<5 | sum(!is.na(y_M$value))<5){
      next
    }
    res = wilcox.test(y_F$value,y_M$value)
    pvalue = round(res$p.value,digits = 4)
    
    ## plot
    y$Sex = as.factor(y$Sex)
    tt = paste(str_wrap(drug_name[i],width = 12,whitespace_only = F), 
               paste0(gene_name," ",moa),
               paste0("(p = ",pvalue,")"),
               sep = "\n")
    
    g = ggboxplot(y,"Sex","value",
                  color = "Sex",
                  palette = c("#F8766D", "#619CFF"),
                  add = "jitter",
                  xlab = F,
                  ylab = paste0("Drug Response Value"),
                  title = tt
    )+theme(legend.position = "none",
            axis.text.x = element_text(size = 10),
            axis.text.y = element_text(size = 10),
            axis.title.y = element_text(size = 10),
            plot.title = element_text(size = 10,face = "bold"))+
      scale_y_continuous(breaks = seq(-3, 1, by = 1), limits = c(-3.2, 1.2))
    
    drug_gene_df = rbind(drug_gene_df,c(drug_name[i],gene_name,pvalue))
    
    if (pvalue<0.05){
      cat("P<0.05: ",drug_name[i],gene_name,moa,"\n")
      ggsave(paste0("plot/boxplot/p005/",gene_name,"_",drug_name[i],"_drug.png"),
             g,height = 2.5,width = 2.5)
    }else if(pvalue<0.1){
      cat("P<0.25: ",drug_name[i],gene_name,moa,"\n")
      ggsave(paste0("plot/boxplot/p01/",gene_name,"_",drug_name[i],"_drug.png"),
             g,height = 2.5,width = 2.5)
    }
    
  }
  
  
}

drug_gene_df = as.data.frame(drug_gene_df)
colnames(drug_gene_df) = c("drug","gene","pvalue")
saveRDS(drug_gene_df,"result/pp/drug_gene_df.rds")

