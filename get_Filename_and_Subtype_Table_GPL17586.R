library(GEOquery)
library(tidyverse)



gse1 <- getGEO("GSE86945", GSEMatrix = TRUE)
df1 <- gse1$GSE86945_series_matrix.txt.gz@phenoData@data[, c("supplementary_file.1", "lehmann subtype:ch1")]
df1 <- df1 %>%
  mutate(supplementary_file.1 = str_extract(supplementary_file.1, "GSM\\w+\\.CEL"))
colnames(df1) <- c("File_Name", "Subtype")



gse2 <- getGEO("GSE86946", GSEMatrix = TRUE)
df2 <- gse2$GSE86946_series_matrix.txt.gz@phenoData@data[, c("supplementary_file.1", "lehmann subtype:ch1")]
df2 <- df2 %>%
  mutate(supplementary_file.1 = str_extract(supplementary_file.1, "GSM\\w+\\.CEL"))
colnames(df2) <- c("File_Name", "Subtype")



gse3 <- getGEO("GSE73540", GSEMatrix = TRUE)
df3 <- gse3$GSE73540_series_matrix.txt.gz@phenoData@data[, c("supplementary_file", "source_name_ch1")]
df3 <- df3 %>%
  mutate(supplementary_file = str_extract(supplementary_file, "GSM\\w+\\.CEL"))
colnames(df3) <- c("File_Name", "Subtype")


Filename_with_Subtype <- rbind(df1, df2, df3)
write.csv(Filename_with_Subtype, "Filename_with_Subtype_GPL17586.csv", row.names = FALSE, quote = FALSE)

