library(affy)
library(GEOquery)
library(tidyverse)
library(oligo)
library(utils)
library(dplyr)
library(biomaRt)
library(limma)
options(timeout = 6000)


#Normalization
Filename_with_Subtype <- read.csv("Filename_with_Subtype_GPL570.csv")
# List all CEL files
celFiles <- list.celfiles(full.names = TRUE, path = "./data")
#Read .CEL files
affyRaw <- read.celfiles(celFiles)
# Perform RMA normalization
eset <- oligo::rma(affyRaw)

normalized_expr <- as.data.frame(exprs(eset))
normalized_expr <- rownames_to_column(normalized_expr, var = "ID")







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


top_genes_BL1vsNormal <- topTable(fit2, coef = "BL1vsNormal", adjust.method = "BH", number = Inf)
top_genes_BL2vsNormal <- topTable(fit2, coef = "BL2vsNormal", adjust.method = "BH", number = Inf)
top_genes_IMvsNormal <- topTable(fit2, coef = "IMvsNormal", adjust.method = "BH", number = Inf)
top_genes_LARvsNormal <- topTable(fit2, coef = "LARvsNormal", adjust.method = "BH", number = Inf)
top_genes_MvsNormal <- topTable(fit2, coef = "MvsNormal", adjust.method = "BH", number = Inf)
top_genes_MSLvsNormal <- topTable(fit2, coef = "MSLvsNormal", adjust.method = "BH", number = Inf)
top_genes_UNSvsNormal <- topTable(fit2, coef = "UNSvsNormal", adjust.method = "BH", number = Inf)






#Extracting Gene Symbols
gse <- getGEO("GPL570", GSEMatrix = TRUE)
# Extracting specific columns "ID" and "Gene Symbol" 
feature_data <- gse@dataTable@table[, c("ID", "Gene Symbol")]

feature_data <- inner_join(feature_data, normalized_expr[,"ID", drop = FALSE], by = "ID")





#Merging the Gene Symbols and wrtitng to CSV files
top_genes_BL1vsNormal <- inner_join(feature_data, top_genes_BL1vsNormal, by = "ID")
top_genes_BL2vsNormal <- inner_join(feature_data, top_genes_BL2vsNormal, by = "ID")
top_genes_IMvsNormal <- inner_join(feature_data, top_genes_IMvsNormal, by = "ID")
top_genes_LARvsNormal <- inner_join(feature_data, top_genes_LARvsNormal, by = "ID")
top_genes_MvsNormal <- inner_join(feature_data, top_genes_MvsNormal, by = "ID")
top_genes_MSLvsNormal <- inner_join(feature_data, top_genes_MSLvsNormal, by = "ID")
top_genes_UNSvsNormal <- inner_join(feature_data, top_genes_UNSvsNormal, by = "ID")

normalized_expr <- inner_join(feature_data, normalized_expr, by = "ID")



write.csv(top_genes_BL1vsNormal, "top_genes_BL1vsNormal_GPL570.csv", row.names = FALSE, quote = FALSE)
write.csv(top_genes_BL2vsNormal, "top_genes_BL2vsNormal_GPL570.csv", row.names = FALSE, quote = FALSE)
write.csv(top_genes_IMvsNormal, "top_genes_IMvsNormal_GPL570.csv", row.names = FALSE, quote = FALSE)
write.csv(top_genes_LARvsNormal, "top_genes_LARvsNormal_GPL570.csv", row.names = FALSE, quote = FALSE)
write.csv(top_genes_MvsNormal, "top_genes_MvsNormal_GPL570.csv", row.names = FALSE, quote = FALSE)
write.csv(top_genes_MSLvsNormal, "top_genes_MSLvsNormal_GPL570.csv", row.names = FALSE, quote = FALSE)
write.csv(top_genes_UNSvsNormal, "top_genes_UNSvsNormal_GPL570.csv", row.names = FALSE, quote = FALSE)

write.csv(normalized_expr, "normalized_GeneSymbol_GPL570.csv", row.names = FALSE, quote = FALSE)



