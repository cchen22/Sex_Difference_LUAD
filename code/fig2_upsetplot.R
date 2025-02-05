library(tidyverse)
library(ggupset)

## data
expr = readRDS("data/pp/luad_expression.rds")
pro = readRDS("data/luad_proteomics_imputed.rds")
phos = readRDS("data/luad_phosphoproteomics_imputed.rds")
ace = readRDS("data/luad_acetylproteomics_imputed.rds")

listInput = list("Genes"=rownames(expr),
                 "Proteins"=rownames(pro),
                 "Phosphorylation"=sub("_.*", "", rownames(phos)),
                 "Acetylation"=sub("_.*", "", rownames(ace)))
listInput = lapply(listInput,unique)
lapply(listInput,length)
names(listInput) = c("RNA [25693]","Protein [10340]",
                     "Phospho-protein [5352]","Acetyl-protein [2129]")


# ## multi-symbol checker: https://www.genenames.org/tools/multi-symbol-checker/
# genes = unique(unlist(listInput))
# write.table(genes,"data/test/all_genes",
#             sep = "\n",quote = F,row.names = F,col.names = F)

library(HGNChelper)
genes = unique(unlist(listInput))
new.hgnc.table = getCurrentHumanMap()
gene_check = checkGeneSymbols(x = genes, 
                              unmapped.as.na = F, 
                              map = new.hgnc.table,
                              species = "human")


listInput = lapply(listInput,function(x){
  plyr::mapvalues(x,
                  gene_check$x,
                  gene_check$Suggested.Symbol,
                  warn_missing = F)
})

lapply(listInput,length)

## convert to ggupset structure
invertedList = lapply(unique(unlist(listInput)), function(x) {
  names(which(sapply(listInput, function(y) x %in% y)))
})
names(invertedList) = unique(unlist(listInput))
upset_data =  tibble(
  Element = names(invertedList),
  Categories = I(invertedList)
)

## plot
gg = ggplot(upset_data,aes(x=Categories))+
  geom_bar() +
  #geom_text(stat='count', aes(label=after_stat(count)), vjust=-1) +
  scale_x_upset(n_intersections = Inf) +
  scale_y_continuous(name = "Gene Size")+
  theme_combmatrix(combmatrix.label.text = element_text(size=24),
                   combmatrix.label.extra_spacing = 10)+
  xlab("") +
  theme(axis.title=element_text(size=24,face="bold"),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(angle = 90, vjust = -20, hjust = 0.5) 
  )

ggsave("plot/upset/luad_landscape_upset.png",gg,width = 8,height = 6)
