# gc()
rm(list=ls())

library(tidyverse)
library(piraid)
library(xpose)

data_dir <- "../data"
nonmem_dir <- "../models"
output_dir <- "../outputs"

%>%
  group_by(ID) %>%
  mutate(
    NewColumn = ifelse(Value > lag(Value), 1, 0)
  )