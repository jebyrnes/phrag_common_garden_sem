source("01_function_parse.R")

exp_data <- read_excel("../data/common_garden_data.xlsx", sheet=2) %>%
  dplyr::rename(SampleID = `#SampleID`)

exp_data_func <- inner_join(exp_data, pathway_data_sum)
