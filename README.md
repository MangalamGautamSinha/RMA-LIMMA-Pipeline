# RMA-LIMMA-Pipeline

### Overview  

This repository contains R scripts for analyzing gene expression data from Affymetrix microarray platforms. It supports three different platforms: **GPL16686**, **GPL570**, and **GPL17586**. The pipeline includes:  

1. **Data Normalization**: Raw CEL files are normalized using RMA.  
2. **Differential Expression Analysis**: Differential expression is calculated using the `limma` package.  
3. **Gene Annotation**: Probe IDs are mapped to gene symbols using platform-specific methods.  
4. **Result Export**: Outputs are saved as CSV files for downstream analysis.  

### Platforms Supported  

- **GPL16686**  
- **GPL570**  
- **GPL17586**  

### Features  

- Preprocess and normalize raw CEL files from multiple datasets.  
- Perform customizable differential expression analysis for various subtype comparisons.  
- Annotate results with gene symbols for biological interpretation.  

### Prerequisites  

- **R** (version ≥ 4.0)  
- RStudio or any R-compatible IDE  

### Required R Packages  

Install required packages using the following commands:  

```r  
install.packages(c("tidyverse", "dplyr", "utils", "stringr"))  

if (!requireNamespace("BiocManager", quietly = TRUE))  
    install.packages("BiocManager")  

BiocManager::install(c("affy", "GEOquery", "oligo", "biomaRt", "limma", "org.Hs.eg.db"))  
```  

### Data Preparation  

- **CEL Files**: Organize CEL files for each platform in their respective `./data` directories.  
- **Subtype Metadata**: Prepare CSV files named as using respective get_Filename_and_Subtype_Table_GPLnnnn scripts:  
  - `get_Filename_and_Subtype_Table_GPL17586`  
  - `Filename_with_Subtype_GPL570.csv`  
  - `Filename_with_Subtype_GPL17586.csv`
- Make sure that the filename of CEL files in the ./data folder matches exactly with `get_Filename_and_Subtype_Table_GPLnnnn` CSV file.
 
### Running the Pipeline  

1. Clone this repository:  
   ```bash  
   git clone https://github.com/your_username/Microarray-Analysis-Pipeline.git  
   cd Microarray-Analysis-Pipeline  
   ```  

2. Run the platform-specific R scripts:  
   - `GPL16686_analysis.R`  
   - `GPL570_analysis.R`  
   - `GPL17586_analysis.R`  

3. Review the output files in the working directory.  

### Output  

Each script generates the following outputs:  

#### Normalized Expression Data  
- `normalized_GeneSymbol_<Platform>.csv`  

#### Differential Expression Results  
- `top_genes_BL1vsNormal_<Platform>.csv`  
- `top_genes_BL2vsNormal_<Platform>.csv`  
- `top_genes_IMvsNormal_<Platform>.csv`  
- `top_genes_LARvsNormal_<Platform>.csv`  
- `top_genes_MvsNormal_<Platform>.csv`  
- `top_genes_MSLvsNormal_<Platform>.csv`  
- `top_genes_UNSvsNormal_<Platform>.csv`  

### Customization  

You can modify the `makeContrasts` section in each script to fit specific research comparisons or adjust parameters in the `topTable` function to filter top genes differently.  

### Contributing  

Contributions are welcome! Fork the repository and submit pull requests to enhance functionality or support additional platforms.  
