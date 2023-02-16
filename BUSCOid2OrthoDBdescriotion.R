library(tidyverse)
library(jsonlite)

setwd("~/Documents/GiantClamBUSCO")
###read data
Tmaxima_missing <- read_lines("Tmaxima_mollusca_missing.tsv")
Tmaxima_full<- read_lines("Tmaxima_mollusca_full.tsv")
inter_metazoa_busco_ids <- read_lines("inter_metazoa_busco_ids.txt")
inter_mollusca_busco_ids <-read_lines("inter_mollusca_busco_ids.txt")

busco_ids <- inter_metazoa_busco_ids
busco_ids <- inter_mollusca_busco_ids
busco_ids <- Tmaxima_missing
busco_ids <- Tmaxima_full

busco_ids <- c("157138at33208", "179485at33208")

names(busco_ids) <- busco_ids # so map_df gets .id right
# read_json doesn't work here
# a tibble of interpro profile associated with each busco_id
busco_ipr <- map_df(busco_ids, .id="busco_id",  function(busco_id){
  write(busco_id, stderr()) # just so we can monitor progress
  # get all info on the orthogroup
  odb_info <- read_json(paste0("https://v100.orthodb.org/group?id=", busco_id),
                        simplifyVector = TRUE)
  # return the interpro table
  #odb_info$data$interpro_domains
  # return the molecularfunction
  odb_info$data$molecular_function
  # return the functional categroy 
  #odb_info$data$functional_category
})

busco_ipr %>% select(busco_id, description)
write.table(busco_ipr, file = "Tmaxima_functional_category_full.txt", row.names=FALSE, sep="\t", quote = FALSE)
write.table(busco_ipr, file = "Tmaxima_functional_category_missing.txt", row.names=FALSE, sep="\t", quote = FALSE)
write.table(busco_ipr, file = "Tmaxima_interpro_domains_full.txt", row.names=FALSE, sep="\t", quote = FALSE)
write.table(busco_ipr, file = "Tmaxima_molecular_function_full.txt", row.names=FALSE, sep="\t", quote = FALSE)


library(ggplot2)
library(dplyr)
# qucik barplot
busco_ipr %>% count(id)
ggplot(busco_ipr %>% count(description), aes(x=description, y=n))+
  geom_bar(stat = "identity")+
  coord_flip()
# by count
busco_ipr %>% 
  count(description) %>%
  ggplot(aes(x=description, y=n, label=n))+
  geom_bar(stat = "identity")+
  geom_text(hjust = -0.2)+
  coord_flip()
# by percentage
busco_ipr %>% 
  count(description) %>%
  mutate(per = prop.table(n)) %>% 
  ggplot(aes(x=description, y=per, label = scales::percent(per)))+
  geom_bar(stat = "identity")+
  geom_text(hjust = 0)+
  coord_flip()
### two dataset on the same plot
full <- read.table("Tmaxima_functional_category_full.txt", sep="\t", header = TRUE)
missing <- read.table("Tmaxima_functional_category_missing.txt", sep="\t", header = TRUE)

full_plot <- full%>% 
  count(description) %>%
  mutate(per = prop.table(n))
full_plot$group <- "full"

missing_plot <- missing%>% 
  count(description) %>%
  mutate(per = prop.table(n))
missing_plot$group <- "missing"


data<-rbind(full_plot, missing_plot)
ggplot(data, aes(x=description, y=per, label = scales::percent(per), fill=group))+
  geom_bar(width = 0.8, stat = "identity", position=position_dodge(width = 0.8))+
  geom_text(hjust = 0, size=3, position=position_dodge(width = 1))+
  coord_flip()


### subset the full dataset to get the missing dataset
