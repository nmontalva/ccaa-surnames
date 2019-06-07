# setwd("~/research/nico/ccaa-surnames/bienes-raices_2017/")
library(tidyverse)
library(reldist)

commoners_csv <- "commoners.csv"
col_names <- c("community", "right_id",
               "firstname", "surname",
               "rights", "shares",
               "commune", "province", "region",
               "firstname1", "firstname2",
               "surname_father", "surname_mother", "sex")
read_commoners_csv <- function(filename=commoners_csv) {
  read_csv(commoners_csv,
           skip=1,
           col_names=col_names,
           col_types="ciccddcccccccc")
}

traits <- function(commoners, by="community") {
  args <- c(list(commoners), by)
  do.call(group_by_, args) %>%
    summarise(N=n(),
              S=n_distinct(surname_father)/N,
              R=mean(rights),
              G=gini(shares), # gini(rights),
              A=mean(sex == "M"))
}
