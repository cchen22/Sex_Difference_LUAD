library(ssGSEA2)
library(biomaRt)
library(cmapR)

# Function to map gene symbols to UniProt IDs
get_uniprot_ids <- function(gene_symbols) {
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  uniprot_mapping <- getBM(attributes = c('external_gene_name', 'uniprotswissprot'),
                           filters = 'external_gene_name',
                           values = gene_symbols,
                           mart = ensembl)
  # Remove duplicates
  uniprot_mapping <- uniprot_mapping[!duplicated(uniprot_mapping$external_gene_name), ]
  return(uniprot_mapping)
}
# Function to convert phosphosite IDs to the required format
convert_phosphosites <- function(phospho_data) {
  gene_symbols <- sapply(strsplit(rownames(phospho_data), "_"), `[`, 1)
  uniprot_mapping <- get_uniprot_ids(unique(gene_symbols))
  processed_identifiers <- character()
  
  # Replace phosphosite IDs with UniProt format
  new_row_names <- vector("character", nrow(phospho_data))
  for (i in seq_len(nrow(phospho_data))) {
    site <- rownames(phospho_data)[i]
    parts <- strsplit(site, "_")[[1]]
    gene_symbol <- parts[1]
    site_position <- parts[2]
    uniprot_id <- uniprot_mapping$uniprotswissprot[uniprot_mapping$external_gene_name == gene_symbol]
    
    if (length(uniprot_id) > 0 && !is.na(uniprot_id[1]) && uniprot_id!="") {
      new_id <- paste0(uniprot_id[1], ";", site_position, "-p")
      if (!(new_id %in% processed_identifiers)) {
        new_row_names[i] <- new_id
        processed_identifiers <- c(processed_identifiers, new_id)
      } else {
        new_row_names[i] <- NA  # Assign NA to duplicates
      }
    } else {
      new_row_names[i] <- NA  # Assign NA if no mapping found
    }
  }
  
  # Update row names, removing NAs
  phospho_data <- phospho_data[!is.na(new_row_names), ]
  rownames(phospho_data) <- new_row_names[!is.na(new_row_names)]
  phospho_data = as.matrix(phospho_data)
  return(phospho_data)
}


## CPTAC LUAD
phos = readRDS("data/luad_phosphoproteomics_imputed.rds")
## convert gene symbol to uniprot id then convert to gct format
phospho_data = convert_phosphosites(phos)
test = GCT(mat = phospho_data)
write_gct(test, "data/PTM-SEA/luad_phos", precision = 4, appenddim = TRUE, ver = 3)

## SSGSEA2, no permutaion
res = run_ssGSEA2("data/PTM-SEA/luad_phos_n111x13338.gct",
                  output.prefix = "luad_phos",
                  gene.set.databases = "data/PTM-SEA/db_ptm.sig.db.all.v2.0.0/ptm.sig.db.all.uniprot.human.v2.0.0.gmt",
                  output.directory = "data/PTM-SEA/",
                  sample.norm.type = "none", 
                  weight = 0.75, 
                  correl.type = "rank", 
                  statistic = "area.under.RES",
                  output.score.type = "NES", 
                  nperm = 0, 
                  min.overlap = 5, 
                  extended.output = TRUE, 
                  global.fdr = FALSE,
                  log.file = "data/PTM-SEA/run.log")



## APOLLO LUAD
phos = readRDS("data/pp/apollo_phos_imputed.rds")
phos = as.matrix(phos)
phospho_data = convert_phosphosites(phos)
test = GCT(mat = phospho_data)
write_gct(test, "data/PTM-SEA/apollo_phos", precision = 4, appenddim = TRUE, ver = 3)

## SSGSEA2, no permutaion
res = run_ssGSEA2("data/PTM-SEA/apollo_phos_n87x1303.gct",
                  output.prefix = "apollo_phos",
                  gene.set.databases = "data/PTM-SEA/db_ptm.sig.db.all.v2.0.0/ptm.sig.db.all.uniprot.human.v2.0.0.gmt",
                  output.directory = "data/PTM-SEA/",
                  sample.norm.type = "none", 
                  weight = 0.75, 
                  correl.type = "rank", 
                  statistic = "area.under.RES",
                  output.score.type = "NES", 
                  nperm = 0, 
                  min.overlap = 1, 
                  extended.output = TRUE, 
                  global.fdr = FALSE,
                  log.file = "data/PTM-SEA/run.log")




