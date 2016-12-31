library(ggplot2)
library(RColorBrewer)
source("./01_function_parse.R")


###Correlation Data Frames
cor_pathways_sum <- data.frame(cor(pathway_data_sum[,-1])) %>%
  mutate(V2 = rownames(.)) %>%
  gather(V1, correlation, -V2)

cor_subpathways_sum <- data.frame(cor(subpathway_data_sum[,-1])) %>%
  mutate(V2 = rownames(.)) %>%
  gather(V1, correlation, -V2)


func_data_transpose <- select(func_data, -Pathway, -Subpathway) %>%
  gather(SampleID, Abundance, -Specific_Pathway) %>%
  spread(Specific_Pathway, Abundance)

cor_raw_pathways <- data.frame(cor(func_data_transpose[,-1])) %>%
  mutate(V2 = rownames(.)) %>%
  gather(V1, correlation, -V2)

###PLOTS
pairs(pathway_data_sum[,-1])

ggplot(data=cor_pathways_sum, mapping=aes(x=V1, y=V2, fill=correlation)) +
  geom_tile() +
  scale_fill_gradientn(colours=brewer.pal(11, "BrBG"), limits=c(-1,1))+ 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))



ggplot(data=cor_subpathways_sum, mapping=aes(x=V1, y=V2, fill=correlation)) +
  geom_tile() +
  scale_fill_gradientn(colours=brewer.pal(11, "BrBG"), limits=c(-1,1))+ 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))


ggplot(data=cor_raw_pathways, mapping=aes(x=V1, y=V2, fill=correlation)) +
  geom_tile() +
  scale_fill_gradientn(colours=brewer.pal(11, "BrBG"), limits=c(-1,1)) #+ 
 # theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))