# Load biomaRt
library(tidyverse)
library(readxl)
library(here)
library(biomaRt)

# Load probe data (columns: ProbeID, ExpressionValue)
probe_data <- read_excel(here("data","ResultsTable_DE.xlsx"))

# 3. Connect to Ensembl BioMart
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

annotation <- getBM(
  attributes = c("agilent_sureprint_g3_ge_8x60k", "external_gene_name"),
  filters = "agilent_sureprint_g3_ge_8x60k",
  values = probe_data$ProbeName,  # REPLACE with your column name
  mart = mart
)

# 5. Merge with expression data
updated_data <- left_join(probe_data, annotation, 
                          by = c("ProbeName" = "agilent_sureprint_g3_ge_8x60k")) %>%  # CHANGE THIS
  filter(
    !is.na(external_gene_name),          # Remove NA values
    external_gene_name != "",              # Remove empty strings
    !str_detect(external_gene_name, "^[-/]")  # Remove placeholders like "-" or "---"
  )

# 6. Save updated data
# write.csv(updated_data, "updated_expression_data.csv", row.names = FALSE)