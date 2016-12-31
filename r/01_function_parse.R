library(dplyr)
library(readxl)
library(tidyr)

#was functions for Jarrett originally, but the data was too correlated
func_data <- read_excel("../data/Final picrust.xlsx")

#split up the pathway types for future grouping
func_data <- separate(func_data, KEGG_Pathways, "; ", into=c("Pathway", "Subpathway", "Specific_Pathway"))
func_data <- func_data[,-1] #get rid of redundant column

#reshape to long format
func_data_long <- gather(func_data, "SampleID", "Abundance", -Pathway, -Subpathway, -Specific_Pathway)
func_data_long$SampleID <- as.character(func_data_long$SampleID)

#some summary data frames
pathway_data <- group_by(func_data_long, SampleID, Pathway) %>%
  dplyr::summarise(sum_abund = sum(Abundance, na.rm=T), mean_abund = mean(Abundance, na.rm=T)) 

#Make it wide again
pathway_data_sum <- dplyr::select(pathway_data, -mean_abund) %>%
  spread(Pathway, sum_abund)


pathway_data_mean <- dplyr::select(pathway_data, -sum_abund) %>%
  spread(Pathway, mean_abund)


subpathway_data <- group_by(func_data_long, SampleID, Subpathway) %>%
  dplyr::summarise(sum_abund = sum(Abundance, na.rm=T), mean_abund = mean(Abundance, na.rm=T))

subpathway_data_sum <- dplyr::select(subpathway_data, -mean_abund) %>%
  spread(Subpathway, sum_abund)
