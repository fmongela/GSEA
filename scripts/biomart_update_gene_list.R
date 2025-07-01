# Try to update info about sequences on the microarray
# using Ensemble fresh data
# Hopefully will help with subsequent GSEA!

# Load biomaRt
library(tidyverse)
library(readxl)
library(here)
library(biomaRt)
library(rstudioapi)

# Load probe data (columns: ProbeID, ExpressionValue)
probe_data <- read_excel(here("data","ResultsTable_DE.xlsx"))

# 3. Connect to Ensembl BioMart

choice <- showQuestion(
  title = "Data Loading",
  message = "Do you want to download data or load a local file?",
  ok = "Download",
  cancel = "Load Local"
)

if (choice) {
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

annotation <- getBM(
  attributes = c(
    "agilent_sureprint_g3_ge_8x60k",
    "external_gene_name",
    "transcript_biotype"
  ),
  filters = c("agilent_sureprint_g3_ge_8x60k", "transcript_biotype"),
  values = list(
    probe_data$ProbeName,
    c("protein_coding", "lncRNA", "miRNA", "snRNA", "snoRNA", "misc_RNA", "rRNA")
  ),
  mart = mart
)
saveRDS(annotation, here("results","annotation_data.rds"))
} else {
  file_path <- selectFile(caption = "Select Saved Data File",    # Interactive file picker
                          path = here("results"))  
  annotation <- readRDS(file_path)
}



# 5. Merge with expression data
updated_data <- left_join(probe_data, annotation, 
                          by = c("ProbeName" = "agilent_sureprint_g3_ge_8x60k")) %>%  # CHANGE THIS
  filter(
    !is.na(external_gene_name),          # Remove NA values
    external_gene_name != "",              # Remove empty strings
    !str_detect(external_gene_name, "^[-/]"),  # Remove placeholders like "-" or "---"
    Description != "Unknown", # Remove unknown description
    transcript_biotype == "protein_coding"
  )

mouse_cols <- c(paste0("controle_", 1:4), paste0("expe_", 1:4))

tibble_for_GSEA <- updated_data %>%
  dplyr::select(
    "external_gene_name", "Description",
    "controle_1", "controle_2", "controle_3", "controle_4",
    "expe_1", "expe_2", "expe_3", "expe_4" 
  ) %>%
  dplyr::rename(Name = external_gene_name) %>%
  group_by(Name) %>%
  summarise(
    Description = first(na.omit(Description)),
    across(all_of(mouse_cols), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

# 6. Save updated data
write_tsv(tibble_for_GSEA, here("results", "updated_expression_data.txt"))