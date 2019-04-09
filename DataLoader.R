# TITLE: DataLoader.R
# AUTHOR: Kevin O'Connor
# DATE CREATED: 3/7/19
# DATE MODIFIED: 4/9/19

DataLoader <- function(file_dir = getwd(), 
                       normalize_quantiles = TRUE, 
                       select_primary = TRUE){
  require(preprocessCore)
  
  stringsAsFactors <- FALSE
  
  target_rnaseq <- read.csv(file.path(file_dir, "target_laml_rnaseq_040119.csv"), header = TRUE)
  target_rnaseq <- target_rnaseq[-which(duplicated(target_rnaseq[, 1])), ]
  rownames(target_rnaseq) <- target_rnaseq[, 1]
  target_rnaseq <- target_rnaseq[, -1]
  
  target_clinical <- read.csv(file.path(file_dir, "target_laml_clinical_040119.csv"), header = TRUE)
  rownames(target_clinical) <- target_clinical[, 1]
  target_clinical <- target_clinical[, -1]
  
  tcga_rnaseq <- read.csv(file.path(file_dir, "tcga_laml_rnaseq_040119.csv"), header = TRUE)
  rownames(tcga_rnaseq) <- tcga_rnaseq[, 1]
  tcga_rnaseq <- tcga_rnaseq[, -1]
  
  tcga_clinical <- read.csv(file.path(file_dir, "tcga_laml_clinical_040119.csv"), header = TRUE)
  rownames(tcga_clinical) <- tcga_clinical[, 1]
  tcga_clinical <- tcga_clinical[, -1]
  
  # Quantile normalization.
  if(normalize_quantiles){
    tcga_rnaseq_unnorm <- tcga_rnaseq
    target_rnaseq_unnorm <- target_rnaseq
    tcga_rnaseq <- as.matrix(tcga_rnaseq) %>% normalize.quantiles() %>% as.data.frame()
    target_rnaseq <- as.matrix(target_rnaseq) %>% normalize.quantiles() %>% as.data.frame()
    rownames(tcga_rnaseq) <- rownames(tcga_rnaseq_unnorm)
    rownames(target_rnaseq) <- rownames(target_rnaseq_unnorm)
    colnames(tcga_rnaseq) <- colnames(tcga_rnaseq_unnorm) %>%
      gsub(pattern="[.]", replacement="-") %>% 
      substr(start=1, stop=12)
    colnames(target_rnaseq) <- colnames(target_rnaseq_unnorm) %>% 
      gsub(pattern="[.]", replacement="-") %>% 
      substr(start=1, stop=16)
  }

  # Select primary samples in target.
  if(select_primary){
    target_rnaseq <- target_rnaseq[, which(substr(colnames(target_rnaseq_unnorm), start=18, stop=20) == "09A")]  
  }
  
  return(list("target_rnaseq"   = target_rnaseq, 
              "target_clinical" = target_clinical, 
              "tcga_rnaseq"     = tcga_rnaseq, 
              "tcga_clinical"   = tcga_clinical))
}