rm(list = ls()) # Clear datalist.

library(here)
library(readr)
library(dplyr)


Summ_Main_table <- read_csv("computed_data/Summ_Main_table.csv")

###########
# ✅ 1. Overrepresentation Analysis (ORA)

sig_genes <- Summ_Main_table$Gene_name[Summ_Main_table$p_value < 0.05 & abs(Summ_Main_table$log2FC) > log2(1.5)]
write(sig_genes, file = here("computed_data", "sig_genes.txt"), sep= "\t" )

###########
# . Gene Set Enrichment Analysis (GSEA)

library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)

res_df <- Summ_Main_table %>%
  dplyr::select(-c(Change, Significant, Color))

res_df$score <- -log10(res_df$p_value) * sign(res_df$log2FC)
gene_list <- res_df$score
names(gene_list) <- res_df$Gene_name
gene_list <- gene_list[!is.na(gene_list)]
gene_list <- gene_list[!duplicated(names(gene_list))]
gene_list <- sort(gene_list, decreasing = TRUE)

gsea_res <- gseGO(geneList     = gene_list,
                  OrgDb        = org.Mm.eg.db,
                  keyType      = "SYMBOL",   # Or "ENTREZID" if needed
                  ont          = "BP",       # Can be "MF" or "CC" too
                  minGSSize    = 10,
                  maxGSSize    = 500,
                  pvalueCutoff = 0.05,
                  verbose      = TRUE)

# View top results
head(gsea_res)

# Plot enrichment of top term
gseaplot2(gsea_res, geneSetID = gsea_res@result$ID[1])

# Dotplot of top GO terms
dotplot(gsea_res, showCategory = 20, split = ".sign") + facet_grid(.~.sign)

write.csv(as.data.frame(gsea_res), "GSEA_GO_results.csv")
