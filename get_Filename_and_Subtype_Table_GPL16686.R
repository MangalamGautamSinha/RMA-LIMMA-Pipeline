library(GEOquery)
library(tidyverse)



gse1 <- getGEO("GSE80180", GSEMatrix = TRUE)
df1 <- gse1$GSE80180_series_matrix.txt.gz@phenoData@data[, c("supplementary_file", "genotype/variation:ch1")]
df1 <- df1 %>%
  mutate(supplementary_file = str_extract(supplementary_file, "GSM\\w+\\.CEL"))
colnames(df1) <- c("File_Name", "Subtype")



gse2 <- getGEO("GSE176128", GSEMatrix = TRUE)
df2 <- gse2$GSE176128_series_matrix.txt.gz@phenoData@data[, c("supplementary_file", "genotype:ch1")]
df2 <- df2 %>%
  mutate(supplementary_file = str_extract(supplementary_file, "GSM\\w+\\.CEL"))
colnames(df2) <- c("File_Name", "Subtype")



gse3 <- getGEO("GSE226289", GSEMatrix = TRUE)
df3 <- gse3$GSE226289_series_matrix.txt.gz@phenoData@data[, c("supplementary_file", "subclassification:ch1")]
df3 <- df3 %>%
  mutate(supplementary_file = str_extract(supplementary_file, "GSM\\d+_[\\w-]+\\.CEL"))
colnames(df3) <- c("File_Name", "Subtype")


Filename_with_Subtype <- rbind(df1, df2, df3)
write.csv(Filename_with_Subtype, "Filename_with_Subtype_GPL16686.csv", row.names = FALSE, quote = FALSE)

