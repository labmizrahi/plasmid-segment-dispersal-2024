library('tidyverse')
library('igraph')
library('gtools')

rm(list=ls())
gc()

#Read in  depth file and filter values for plasmids present in each sample 
depth <- read.csv("rel_abd_filtered_at_1_percent.txt", sep='\t', header = T) %>%
  rename(plasmids=X) %>% 
  pivot_longer(!plasmids, names_to = 'samples', values_to = 'value') %>%
  filter(value>0)

#Read in plasmid cluster file
segment <- read.delim("plasmid_segments.txt",sep = '\t',header = F) %>% 
    `colnames<-` (c('plasmids','segments'))

#Merge plasmid depth in humans and plasmid segments 
Merged_segments = as_tibble(inner_join(x = segment,y = depth, by = 'plasmids') )

#Get df of segments per human
df <- select(Merged_segments,  3, 2) %>% 
  distinct()

#Create a graph object from df
g <- graph_from_data_frame(df, directed = FALSE, vertices = NULL)
V(g)$type <- bipartite_mapping(g)$type

#Save network 
saveRDS(g, file="original_bipartite_network.RDS")

