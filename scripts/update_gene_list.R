rm(list = ls()) # Clear datalist.

# Lib

library(here)
library(readr)
library(stringr) 
library(readxl)
library(dplyr)
library(HGNChelper)

# load microarray results
ResultsTable_DE <- read_excel(here("data","ResultsTable_DE.xlsx"))
old_gene_name <- ResultsTable_DE %>% pull(GeneName)

# Check and update the gene symbols
updated_symbols <- checkGeneSymbols(old_gene_name, species = "mouse")

new_tibble <- ResultsTable_DE %>%
  mutate(HGNC_name = updated_symbols$Suggested.Symbol) %>% 
  filter(complete.cases(.)) %>%
  dplyr::select(HGNC_name, logFC_expe_controle, adj.P.Val_expe_controle)

ranked_gene_list <- new_tibble %>%
  mutate(rank = sign(logFC_expe_controle ) * -log10(adj.P.Val_expe_controle)) %>%
  dplyr::select(HGNC_name, rank)

final_list <- ranked_gene_list %>% 
  group_by(HGNC_name) %>%
  summarize(Average_rank = mean(rank, na.rm = TRUE)) %>%
  filter(!str_detect(HGNC_name, "Rik$"))

write_tsv(final_list, here("results","ranked_list.rnk"))
  
# The 'updated_symbols' dataframe will contain the original input,
# a boolean indicating if the original symbol was approved,
# and the suggested updated symbol.
# print(updated_symbols)