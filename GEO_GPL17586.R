install.packages(c("tidyverse", "utils", "dplyr", "stringr"))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("affy", "GEOquery", "oligo", "biomaRt", "limma"))



library(affy)
library(GEOquery)
library(tidyverse)
library(oligo)
library(utils)
library(dplyr)
library(biomaRt)
library(limma)
library(stringr)
options(timeout = 6000)


Filename_with_Subtype <- read.csv("Filename_with_Subtype_GPL17586.csv")

#Normalize .CEL files from data1
celFiles <- list.celfiles(full.names = TRUE, path="./data1")
affyRaw <- read.celfiles(celFiles)
eset <- oligo::rma(affyRaw)
normalized_expr1 <- as.data.frame(exprs(eset))
normalized_expr1 <- rownames_to_column(normalized_expr1, var = "ID")

#Normalize .CEL files from data2
celFiles <- list.celfiles(full.names = TRUE, path="./data2")
affyRaw <- read.celfiles(celFiles)
eset <- oligo::rma(affyRaw)
normalized_expr2 <- as.data.frame(exprs(eset))
normalized_expr2 <- rownames_to_column(normalized_expr2, var = "ID")

#Normalize .CEL files from data3
celFiles <- list.celfiles(full.names = TRUE, path="./data3")
affyRaw <- read.celfiles(celFiles)
eset <- oligo::rma(affyRaw)
normalized_expr3 <- as.data.frame(exprs(eset))
normalized_expr3 <- rownames_to_column(normalized_expr3, var = "ID")

#Merge the normalized data frames
normalized_expr <- normalized_expr1 %>%
  inner_join(normalized_expr2, by = "ID") %>%
  inner_join(normalized_expr3, by = "ID")








#limma
group <- factor(Filename_with_Subtype$Subtype)
# Create the design matrix for differential expression
design <- model.matrix(~0 + group)  # Using ~0 removes the intercept, allowing direct contrasts
colnames(design) <- levels(group)

fit <- lmFit(normalized_expr, design)
# Define contrasts for the comparisons you want (e.g., Treatment1 vs Control, Treatment2 vs Control)
# You can also add Treatment1 vs Treatment2 if needed
contrast.matrix <- makeContrasts(
  BL1vsNormal = BL1 - Normal,
  BL2vsNormal = BL2 - Normal,
  IMvsNormal = IM - Normal,
  LARvsNormal = LAR - Normal,
  MvsNormal = M - Normal,
  MSLvsNormal = MSL - Normal,
  UNSvsNormal = UNS - Normal,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
# Extract top differentially expressed genes for each comparison
top_genes_BL1vsNormal <- topTable(fit2, coef = "BL1vsNormal", adjust.method = "BH", number = Inf)
top_genes_BL2vsNormal <- topTable(fit2, coef = "BL2vsNormal", adjust.method = "BH", number = Inf)
top_genes_IMvsNormal <- topTable(fit2, coef = "IMvsNormal", adjust.method = "BH", number = Inf)
top_genes_LARvsNormal <- topTable(fit2, coef = "LARvsNormal", adjust.method = "BH", number = Inf)
top_genes_MvsNormal <- topTable(fit2, coef = "MvsNormal", adjust.method = "BH", number = Inf)
top_genes_MSLvsNormal <- topTable(fit2, coef = "MSLvsNormal", adjust.method = "BH", number = Inf)
top_genes_UNSvsNormal <- topTable(fit2, coef = "UNSvsNormal", adjust.method = "BH", number = Inf)







# Fetching Gene Symbols
gse <- getGEO("GPL17586", GSEMatrix = TRUE)
# Extracting specific columns "ID" and "gene_assignment" 
feature_data <- gse@dataTable@table[, c("ID", "gene_assignment")]

# Extract and replace the text between the first and second '//' in the same column
feature_data$gene_assignment <- trimws(str_extract(feature_data$gene_assignment, "(?<=//)[^/]*(?=//)"))

# perform inner join
feature_data <- inner_join(feature_data, normalized_expr[,"ID", drop = FALSE], by = "ID")
colnames(feature_data) <- c("ID", "Gene_Symbol")










#Merging the Gene Symbols and wrtitng to CSV files
top_genes_BL1vsNormal <- inner_join(feature_data, top_genes_BL1vsNormal, by = "ID")
top_genes_BL2vsNormal <- inner_join(feature_data, top_genes_BL2vsNormal, by = "ID")
top_genes_IMvsNormal <- inner_join(feature_data, top_genes_IMvsNormal, by = "ID")
top_genes_LARvsNormal <- inner_join(feature_data, top_genes_LARvsNormal, by = "ID")
top_genes_MvsNormal <- inner_join(feature_data, top_genes_MvsNormal, by = "ID")
top_genes_MSLvsNormal <- inner_join(feature_data, top_genes_MSLvsNormal, by = "ID")
top_genes_UNSvsNormal <- inner_join(feature_data, top_genes_UNSvsNormal, by = "ID")

normalized_expr <- inner_join(feature_data, normalized_expr, by = "ID")


write.csv(top_genes_BL1vsNormal, "top_genes_BL1vsNormal_GPL17586.csv", row.names = FALSE, quote = FALSE)
write.csv(top_genes_BL2vsNormal, "top_genes_BL2vsNormal_GPL17586.csv", row.names = FALSE, quote = FALSE)
write.csv(top_genes_IMvsNormal, "top_genes_IMvsNormal_GPL17586.csv", row.names = FALSE, quote = FALSE)
write.csv(top_genes_LARvsNormal, "top_genes_LARvsNormal_GPL17586.csv", row.names = FALSE, quote = FALSE)
write.csv(top_genes_MvsNormal, "top_genes_MvsNormal_GPL17586.csv", row.names = FALSE, quote = FALSE)
write.csv(top_genes_MSLvsNormal, "top_genes_MSLvsNormal_GPL17586.csv", row.names = FALSE, quote = FALSE)
write.csv(top_genes_UNSvsNormal, "top_genes_UNSvsNormal_GPL17586.csv", row.names = FALSE, quote = FALSE)

write.csv(normalized_expr, "normalized_GeneSymbol_GPL17586.csv", row.names = FALSE, quote = FALSE)

