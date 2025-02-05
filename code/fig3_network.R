library(igraph)
library(ggraph)
library(biomaRt)
library(tidygraph)
library(tidyverse)
library(graphlayouts)
library(ggforce)
library(scatterpie)


###### Full network ###########
## data
paths = readRDS("result/pp/pathway_paths.rds")
POI_sub = unique(unlist(paths))
TF_luad = readRDS("result/pp/TF_luad.rds")
kin_luad = readRDS("result/pp/kin_luad.rds")

## define TFs, proteins, and kinases
TFs = POI_sub[POI_sub%in%TF_luad] 
Kinases = POI_sub[POI_sub%in%kin_luad] 
Proteins = setdiff(POI_sub,c(Kinases,TFs))

## remove paths that have TFs or kinases in the middle
paths =
  paths %>%
  # removing single node paths
  discard(
    function(p){
      length(p) == 1
    }
  ) %>%
  # middle nodes are Kinases or TFs
  discard(
    function(p){
      any(p[-c(1,length(p))] %in% c(TFs,Kinases))
    }
  )

## network edges
edges <- lapply(paths, function(path) {
  if (length(path) > 1) {
    # Create edges for each consecutive pair of genes in the path
    return(data.frame(from = head(path, -1), to = tail(path, -1)))
  } else {
    return(NULL)  # Ignore paths with a single gene
  }
})
edges <- do.call(rbind, edges)
edges <- unique(edges)
#saveRDS(edges,"result/pp/pathway_network.rds")

## igraph
g = graph_from_data_frame(edges, directed = TRUE)
g = as_tbl_graph(g)
V(g)$node_type <- ifelse(V(g)$name%in%TFs, "TFs",
                         ifelse(V(g)$name%in%Proteins, "Proteins", "Kinases"))


# ## tree layout
# tree_layout = layout_as_tree(g)
# tree_layout = data.frame('name'=V(g)$name,'x1' = tree_layout[,1],'y1'= tree_layout[,2])

## mannual layout
level_K = 16
level_P = 8
level_T = 0
manual_layout = tibble::tribble(
  ~name, ~x, ~y,
  "AURKA", 2, level_K,
  "MAPK14", 5, level_K,
  "GSK3B", 0.25, level_P-0.25,
  "NFKBIA", -4, level_P+0.5,
  "IKBKB", -2.75, level_P,
  "PRKACA", 1, level_P-2,
  "RELA", -1.75, level_P-2,
  "RPS6KA5", -1, level_P+0.75,
  "MAPK3", 2.5, level_P-0.5,
  "AR", -1.2, level_P-0.8,
  "CREB1", 8.5, level_P-1,
  "CRTC2", 8, level_P-4,
  "FKBP5", 3.5, level_P+0.6,
  "MAPK1", 4.5, level_P-0.5,
  "NCOA2", 5.25, level_P+0.6,
  "PPARG", 6.5, level_P,
  "SGK1", 7.8, level_P+0.9,
  "CSNK2A1", 1.9, level_P+1.5,
  "NFKB1", -7, level_T,
  "GLIS2", -5.5, level_T,
  "NFATC1", 10, level_T,
  "NR1D1", -4, level_T,
  "NR3C1", 1, level_T,
  "HNF4A", 7.75, level_T,
  "CREB3", -1, level_T,
  "SP1", -2.75, level_T,
  "ETV6", 6, level_T,
)

#layout = inner_join(manual_layout,tree_layout,by="name")


#saveRDS(g,"result/pp/igraph_pathway.rds")

g = g%>%
  tidygraph::inner_join(manual_layout,by="name")

gg = ggraph(g,x=x,y=y) +
  geom_edge_arc2(
    # geom_edge_link(
    #aes(colour = signed_weight),
    arrow = arrow(
      angle = 15,
      length = unit(0.13, "inches"),
      # ends = "last",
      type = "closed"
    ),
    alpha = .4,
    strength = 0.1,
    start_cap = circle(2.6, 'mm'),
    end_cap = circle(2.6, 'mm')
  ) +
  geom_node_point(alpha=0) +
  geom_node_text(
    aes(label = name,
        #color = node_type
    ),
    repel = F,   # Enable repelling
    fontface = "bold",  # Optional: for better visibility
    size = 6,
    check_overlap = T,
  ) +
  #scale_color_discrete(guide = guide_legend(title = 'Node type')) +
  theme_bw() +
  xlab("") +
  ylab("") +
  ggtitle("")+
  theme_minimal()+
  theme(
    legend.text = element_text(size = 12), # Increase the font size of legend text
    legend.title = element_text(size = 14), # Optionally, increase the font size of the legend title
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
  )
gg


#### add bar chart to each node ####


## node bar info: KA, TFA, and protein abundance
DATs = readRDS("result/pp/luad_dats_limma.rds") #TF
DEPs = readRDS("result/pp/luad_deps_limma.rds") #pro
DAKs = readRDS("result/pp/luad_daks_limma.rds") #kinase

DATs = DATs[DATs$`gene symbol`%in%V(g)$name,]
TFA = (-log(DATs$P.Value))*sign(DATs$logFC)
names(TFA) = DATs$`gene symbol`
df1 = data.frame(ID = names(TFA), TFA = TFA)

DAKs = DAKs[DAKs$`gene symbol`%in%V(g)$name,]
KA = (-log(DAKs$P.Value))*sign(DAKs$logFC)
names(KA) = DAKs$`gene symbol`
df2 = data.frame(ID = names(KA), KA = KA)

DEPs = DEPs[DEPs$`gene symbol`%in%V(g)$name,]
PE = (-log(DEPs$P.Value))*sign(DEPs$logFC)
names(PE) = DEPs$`gene symbol`
df3 = data.frame(ID = names(PE), PE = PE)

df = full_join(df1, df2, by = "ID")
df = full_join(df, df3, by = "ID")
df[is.na(df)] = 0

node_wide = data.frame(ID=V(g)$name, x=V(g)$x,y=V(g)$y)
node_wide = full_join(df,node_wide,by = "ID")
node_long = node_wide %>% 
  dplyr::select(ID:PE) %>%
  rename(node_name = ID) %>%
  mutate(id = 1:nrow(node_wide)) %>%
  gather("attr", "value", TFA:PE)


# create the heatmap charts
heatmap_list <- lapply(1:vcount(g), function(i) {
  tmp = node_long[node_long$id == i, ]
  tmp$value = ifelse(tmp$value > 0, "F", ifelse(tmp$value<0,"M","N"))
  tmp$attr = factor(tmp$attr,levels = c("PE","TFA","KA"))
  gt_plot <- ggplotGrob(ggplot(tmp, aes(x = attr, 
                                        y = 1, 
                                        fill = value)) +
                          geom_tile(color = "grey",size = 0.5) + # Set the border color
                          scale_fill_manual(values = c("F" = scales::alpha("#F8766D", 1), 
                                                       "M" = scales::alpha("#619CFF", 1), 
                                                       "N" = scales::alpha("lightgrey", 1)))+
                        
                          theme_minimal() +
                          theme(
                            legend.position = "none",
                            axis.ticks.x = element_blank(),
                            axis.ticks.y = element_blank(),
                            axis.text.x = element_blank(),  # Remove x-axis text
                            axis.title.x = element_blank(), # Remove x-axis title
                            axis.title.y = element_blank(),
                            panel.grid = element_blank(),
                            axis.text.y = element_blank()
                          ) 
  )
  panel_coords <- gt_plot$layout[gt_plot$layout$name == "panel", ]
  gt_plot[panel_coords$t:panel_coords$b, panel_coords$l:panel_coords$r]
})

# convert to custom annotation
annot_list <- lapply(1:vcount(g), function(i) {
  xmin <- node_wide$x[i] - 0.45
  xmax <- node_wide$x[i] + 0.45
  ymin <- node_wide$y[i] - 0.8
  ymax <- node_wide$y[i] -0.3
  annotation_custom(
    heatmap_list[[i]],
    xmin = xmin,
    xmax = xmax,
    ymin = ymin,
    ymax = ymax
  )
})


# put everything together
final_g = Reduce("+", annot_list, gg) + theme(legend.position = "bottom")
final_g
ggsave("plot/network/network_heat.png",final_g,width = 12,height = 8)




