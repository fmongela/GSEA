rm(list = ls()) # Clear datalist.

library(here)
library(readr)
library(dplyr)

# load microarray results
ResultsTable_DE <- read_excel(here("data","ResultsTable_DE.xlsx"))

load agilent chip decription 
chip_list <- readr::read_delim(
  "D:/smyl/Labo/Mouse/KO_NCL/Transcriptomics/Agilent/028005/028005_mm1_1747914873788/GeneList/028005_D_GeneList_20240520.txt",
  delim = "\t")

list_expe <- select(ResultsTable_DE, ProbeName)
list_ref <- select(chip_list, ProbeID) %>% rename(ProbeName = ProbeID)
dplyr::setequal(list_expe,list_ref)