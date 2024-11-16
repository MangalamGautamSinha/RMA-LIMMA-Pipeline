library(GEOquery)
library(tidyverse)



gse1 <- getGEO("GSE95700", GSEMatrix = TRUE)
df1 <- gse1$GSE95700_series_matrix.txt.gz@phenoData@data[, c("supplementary_file", "bd lehmann subtype:ch1")]
df1 <- df1 %>%
  mutate(supplementary_file = str_extract(supplementary_file, "GSM\\w+\\.CEL"))
colnames(df1) <- c("File_Name", "Subtype")



gse2 <- getGEO("GSE167213", GSEMatrix = TRUE)
df2 <- gse2$GSE167213_series_matrix.txt.gz@phenoData@data[, c("supplementary_file", "tissue:ch1")]
df2 <- df2 %>%
  mutate(supplementary_file = str_extract(supplementary_file, "GSM\\w+\\.CEL"))
colnames(df2) <- c("File_Name", "Subtype")



gse3 <- getGEO("GSE7904", GSEMatrix = TRUE)
df3 <- gse3$GSE7904_series_matrix.txt.gz@phenoData@data[, c("supplementary_file", "characteristics_ch1")]
df3 <- df3 %>%
  mutate(supplementary_file = str_extract(supplementary_file, "GSM\\w+\\.CEL"))
colnames(df3) <- c("File_Name", "Subtype")


Filename_with_Subtype <- rbind(df1, df2, df3)
write.csv(Filename_with_Subtype, "Filename_with_Subtype_GPL570.csv", row.names = FALSE, quote = FALSE)

