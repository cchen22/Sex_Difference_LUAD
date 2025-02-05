# library
library(tidyverse)
library(viridis)

# TCGA
tcga = readRDS("data/pp/luad_rse.rds")
sex_t = tcga@colData %>% 
  as_tibble() %>%
  mutate(value = ifelse(tcga.gdc_cases.demographic.gender=="female","Female",
                        ifelse(tcga.gdc_cases.demographic.gender=="male","Male","Unknown sex"))) %>% 
  dplyr::count(value) %>%
  mutate(variable = "sex") 

smoke_t = tcga@colData %>% 
  as_tibble() %>%
  mutate(value = ifelse(is.na(tcga.xml_tobacco_smoking_history),"Unknown smoking",
                        ifelse(tcga.xml_tobacco_smoking_history==1,"Never smokers","Smokers"))) %>%
  dplyr::count(value) %>% 
  mutate(variable = "smoke") 

race_t = tcga@colData %>% 
  as_tibble() %>%
  mutate(value = ifelse(tcga.gdc_cases.demographic.race=="black or african american","Black or African American",
                        ifelse(tcga.gdc_cases.demographic.race=="white","White",
                               ifelse(tcga.gdc_cases.demographic.race=="not reported","Unknown race","Others")
                        ))) %>%
  dplyr::count(value) %>% 
  mutate(variable = "race")
stage_t = tcga@colData %>% 
  as_tibble() %>%
  mutate(value = ifelse(tcga.gdc_cases.diagnoses.tumor_stage%in%c("stage i","stage ia","stage ib"),"Stage I",
                        ifelse(tcga.gdc_cases.diagnoses.tumor_stage%in%c("stage ii","stage iia","stage iib"),"Stage II",
                               ifelse(tcga.gdc_cases.diagnoses.tumor_stage%in%c("stage iii","stage iiia","stage iiib"),"Stage III",
                                      ifelse(tcga.gdc_cases.diagnoses.tumor_stage%in%c("stage iv"),"Stage IV","Unknown stage"))))) %>%
  dplyr::count(value) %>% 
  mutate(variable = "stage")
df_t = rbind(sex_t,smoke_t,race_t,stage_t)

### CPTAC
cptac = readRDS("data/luad_cptac_clinic.rds")

sex_c = cptac %>% 
  as_tibble(.name_repair = "unique") %>%
  dplyr::count(sex) %>% 
  dplyr::mutate(variable = "sex") %>%
  dplyr::rename(value = sex)

smoke_c = cptac %>% 
  as_tibble(.name_repair = "unique") %>%
  mutate(value = ifelse(tobacco_smoking_history=="Smoking history not available","Unknown smoking",
                        ifelse(tobacco_smoking_history=="Lifelong non-smoker: Less than 100 cigarettes smoked in lifetime","Never smokers","Smokers"))) %>% 
  dplyr::count(value) %>%
  mutate(variable = "smoke")

race_c = cptac %>% 
  as_tibble(.name_repair = "unique") %>%
  mutate(value = ifelse(race=="Unknown","Unknown race",
                        ifelse(race=="White","White",
                               ifelse(race=="Black or African American","Black or African American","Others")
                        ))) %>%
  dplyr::count(value) %>% 
  mutate(variable = "race")

stage_c = cptac %>% 
  as_tibble(.name_repair = "unique") %>%
  mutate(value = tumor_stage_pathological) %>%
  dplyr::count(value) %>% 
  mutate(variable = "stage") %>%
  add_row(value = "Unknown stage",n = 0,variable = "stage")

df_c = rbind(sex_c,smoke_c,race_c,stage_c)

## merge two tibble
df_t = df_t %>% mutate(group = "TCGA")
df_c = df_c %>% mutate(group = "CPTAC")
df = bind_rows(df_t,df_c)

## plot
df$value = factor(df$value,levels = unique(df$value))
df$variable = factor(df$variable,levels = unique(df$variable))
df$group = factor(df$group)
colnames(df) = c("category","count","variable","group")

df$category = factor(df$category,levels = c("Female", "Male", "Never smokers", "Smokers", "Unknown smoking", 
                                            "Black or African American", "White", "Others", "Unknown race",  
                                            "Stage I", "Stage II", "Stage III", "Stage IV", "Unknown stage")
)

my_colors <- c(
  # Subgroup of 2
  "#4F81BD", "#1F3864", 
  # Subgroup of 3
  "#F4B084", "#E46C0A", "#984807", 
  # Subgroup of 4
  "#FFF2CC", "#FFD966", "#FFC000","#BF9000", 
  # Subgroup of 5
  "#C6E0B4", "#A9D08E", "#70AD47", "#548235", "#375623"
)

gg=ggplot(df) +      
  geom_bar(aes(x=variable, y=count, fill=category), 
           stat="identity", 
           colour = "grey50",
           position = "stack",
           alpha=0.8) +
  scale_fill_manual(values = my_colors) +
  facet_grid(~group)+
  xlab("") +
  ylab("Sample Size")+
  theme(axis.title=element_text(size=24,face="bold"),
        axis.text.y = element_text(size=20),
        axis.text.x = element_text(size=20,angle = 45,hjust = 1),
        legend.text = element_text(size = 20),
        legend.title = element_text(size=24,face = "bold"),
        strip.text = element_text(size = 24,face = "bold"))

ggsave("plot/barplot/luad_landscape_stackedbar.png",gg,width = 8,height = 6)


## chang to percentage plot
df = df %>%
  group_by(group, variable) %>%
  mutate(total = sum(count), 
         percentage = count / total * 100) %>%
  ungroup()

gg = ggplot(df) +      
  geom_bar(aes(x = variable, y = percentage, fill = category), 
           stat = "identity", 
           colour = "grey50",
           position = "stack",
           alpha = 0.8) +
  scale_fill_manual(values = my_colors) +
  facet_grid(~group,labeller = labeller(group=c("CPTAC"="CPTAC\n(n=111)",
                                                "TCGA"="TCGA\n(n=502)"))) +
  xlab("") +
  ylab("Percentage") +
  theme(axis.title = element_text(size = 24, face = "bold"),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 24, face = "bold"),
        strip.text = element_text(size = 24, face = "bold"))
ggsave("plot/barplot/luad_landscape_stackedbar_percent.png",gg,width = 8,height = 6)


