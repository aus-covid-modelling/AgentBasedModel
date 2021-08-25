#library(tidyverse)
library(tidyr)
library(readr)
library(dplyr)
library(rjson)
library(data.table)

args <- commandArgs(trailingOnly = T)
max_t  = 600
i <- 1
parfile <- fromJSON(file=paste0(args[1], ".json"))
folder <- paste0(parfile$output_directory, parfile$folder_name)

num_sims <- length(list.files(folder))
sim_numbers <- gsub(pattern = "sim_number_|.csv", "", list.files(folder)) %>% as.numeric()

completed_df <- lapply(sim_numbers, function(i) {
  file_name <- paste(folder,"/sim_number_",i,".csv",sep = "")
  individuals <- fread(file_name)
  infected_individuals <- individuals %>% dplyr::filter(!is.na(`Time of exposure`))

  infected_individuals <- mutate(infected_individuals, `Age bracket` = cut(Age, breaks = c(seq(0,80,by = 5),Inf),include.lowest = TRUE, right= FALSE)) 
  
  infected_individuals <- mutate(infected_individuals, `Date symptoms` = cut(`Time of symptom onset`, breaks = seq(0,max_t,by = 1),include.lowest = FALSE, right= TRUE,labels = FALSE)) #Currently ignores symptoms at 0. 
  
  summary_df <- dplyr::filter(infected_individuals) %>% group_by(`Date symptoms`, `Age bracket`,`Vaccine at infection`, `Doses at infection`,Symptomatic) %>% summarise(incidence = n()) #Filter it to be only symptomatic infections.
  
  completed_df <-summary_df %>% ungroup() %>%
  complete(expand(infected_individuals,`Date symptoms` = seq(0,max_t,by=1), `Age bracket`,`Vaccine at infection` = c("None","Pfizer","Moderna","AstraZeneca"), Symptomatic,`Doses at infection` = c(0,1,2)),fill = list(incidence = 0))%>% mutate(Sim = i)
}) %>% rbindlist()


# ggplot(dplyr::filter(individuals,`Time of symptom onset` < 600),aes(x = `Time of symptom onset`,fill = `Symptomatic`)) + geom_histogram(binwidth = 1) + scale_y_continuous("Incidence")

# ggplot(dplyr::filter(individuals,`Time of symptom onset` < 600),aes(x = `Time of infection`,y = `Secondary Infections`)) + geom_point() + geom_smooth()

completed_df <- dplyr::filter(completed_df, !(`Vaccine at infection`=="None"&(`Doses at infection`!=0)))
completed_df <- dplyr::filter(completed_df, !(`Vaccine at infection`!="None"&(`Doses at infection`==0)))
completed_df <- dplyr::filter(completed_df, !is.na(`Date symptoms`))
wide_complete <- pivot_wider(completed_df,names_from = Sim,values_from = incidence)
fwrite(wide_complete,paste(folder,"_summary.csv",sep = ""))
fwrite(completed_df,paste(folder,"_summary_long.csv",sep = ""))


