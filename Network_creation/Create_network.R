library('tidyverse')
library('igraph')
library('gtools')

rm(list=ls())
gc()

#read in metadata file and select relevant columns only 
metadata <- read.delim("Final-metadata-networks.txt",sep = '\t',header = T) %>% 
  mutate(across(c(Country, Continent, Disease_Type, Disease_Category, Disease_State), ~ gsub(" ", "_", .x))) %>% 
  as.data.frame() %>% 
  `rownames<-`(.[,1]) %>% 
  select(c(Sample_ID, Country, Continent, Disease_Type, Disease_Category, Disease_State)) #%>% 

#Select depth values by filtered metadata values 
depth <- read.csv("rel_abd_filtered_at_1_percent.txt", sep='\t', header = T, row.names = 1) %>%  
  select(any_of(rownames(metadata))) %>% 
  rownames_to_column(var = 'plasmids') %>% 
  pivot_longer(!plasmids, names_to = 'samples', values_to = 'value') %>%
  filter(value>0)


#import plasmid segments 
segment <- read.delim("unique_plasmid_clusters.txt",sep = '\t',header = F) %>% 
    `colnames<-` (c('plasmids','Segment'))

#merge plasmid depth in humans and plasmid segments to 1 data frame
Merged = as_tibble(merge.data.frame(x = segment,y = depth, by = 'plasmids') )

#add human metadata to merge data of plasmids, humans and segments to one data frame
Merged_segments = as_tibble(inner_join(x = Merged,y =metadata ,by= c('samples'='Sample_ID') ,all.x = T))

#Deduplicate the entries to create binary edges (not weighted)
df <- select(Merged_segments,  3, 2) %>% 
  distinct()

#creating a graph object from dataframe
g <- graph_from_data_frame(df, directed = FALSE, vertices = NULL)

V(g)$type <- bipartite_mapping(g)$type

saveRDS(g, file="original_bipartite_network.RDS")

#create the projections
g_projected <- bipartite_projection(g, multiplicity = TRUE)

#human projection 
humans_projection <- get.data.frame(g_projected$proj1) %>% 
  rename(source=from, target=to) %>% 
  filter(humans_projection$weight>100)

write.table(humans_projection, "projected_network.csv",  sep = "\t", col.names = T, row.names = F, quote=F) 
