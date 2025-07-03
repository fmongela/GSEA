# Try to update info about sequences on the microarray
# using Ensemble fresh data
# Hopefully will help with subsequent GSEA

# Load biomaRt !
library(clipr)
library(tidyverse)
library(readxl)
library(here)
library(biomaRt)
library(rstudioapi)

# Load probe data (columns: ProbeID, ExpressionValue)
probe_data <- read_excel(here("data","ResultsTable_DE.xlsx"))

# 3. Choice to connect to Ensembl BioMart or load available mart
# if downloaded, the mart in saved locally for further use 

choice <- showQuestion(
  title = "Data Loading",
  message = "Download mart from Ensembl OR load a local file?",
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

# 5. Merge mar with expression data

updated_data <- left_join(probe_data,
                          annotation,
                          by = c("ProbeName" = "agilent_sureprint_g3_ge_8x60k")) %>%
  filter(
    !is.na(external_gene_name),                # Remove NA values
    external_gene_name != "",                  # Remove empty strings
    !str_detect(external_gene_name, "^[-/]"),  # Remove placeholders like "-" or "---"
    Description != "Unknown" # Remove unknown description
  )

# Define a lncRNA and a protein subsets
    
lncRNA_GSEA_set <- updated_data %>%
  filter(transcript_biotype == "lncRNA",
    adj.P.Val_expe_controle <= 1
  )

Protein_coding_GSEA_set <- updated_data %>%
  filter(transcript_biotype == "protein_coding",
         adj.P.Val_expe_controle <= 1
  )

# Take care of duplicates

mouse_cols <- c(paste0("controle_", 1:4), paste0("expe_", 1:4), "adj.P.Val_expe_controle")

lncRNA_GSEA_set <- lncRNA_GSEA_set %>%
  dplyr::select(
    "external_gene_name", "Description",
    "controle_1", "controle_2", "controle_3", "controle_4",
    "expe_1", "expe_2", "expe_3", "expe_4" ,"adj.P.Val_expe_controle"
  ) %>%
  dplyr::rename(Name = external_gene_name) %>%
  group_by(Name) %>%
  summarise(
    Description = first(na.omit(Description)),
    across(all_of(mouse_cols), ~ median(.x, na.rm = TRUE)),
    .groups = "drop"
  )%>%
  arrange(adj.P.Val_expe_controle) %>%
  dplyr::select(-adj.P.Val_expe_controle)

Protein_coding_GSEA_set <- Protein_coding_GSEA_set %>%
  dplyr::select(
    "external_gene_name", "Description",
    "controle_1", "controle_2", "controle_3", "controle_4",
    "expe_1", "expe_2", "expe_3", "expe_4" ,"adj.P.Val_expe_controle"
  ) %>%
  dplyr::rename(Name = external_gene_name) %>%
  group_by(Name) %>%
  summarise(
    Description = first(na.omit(Description)),
    across(all_of(mouse_cols), ~ median(.x, na.rm = TRUE)),
    .groups = "drop"
  )%>%
  arrange(adj.P.Val_expe_controle) %>%
  dplyr::select(-adj.P.Val_expe_controle)

# 6. Save subset (lncRNA and proteins) data as comma separated bare list of gene names
# or use in stuff like https://singlecell.broadinstitute.org

write.table (
  t(as.matrix(Protein_coding_GSEA_set$Name)),
  file = here("results", "Protein_coding_GSEA_set.csv"),
  sep       = ",",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
write.table (
  t(as.matrix(lncRNA_GSEA_set$Name)),
  file = here("results", "lncRNA_GSEA_set.csv"),
  sep       = ",",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

# 6. Save full protein coding dataset as tab separated table
# extension is .txt for GSEA (maybe not important...)
write_tsv(Protein_coding_GSEA_set,
            here("results","Protein_coding_GSEA_full_set.txt"))

