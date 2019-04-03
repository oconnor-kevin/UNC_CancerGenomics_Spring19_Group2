# TITLE: DataLoader.R
# AUTHOR: Kevin O'Connor
# DATE CREATED: 3/7/19
# DATE MODIFIED: 4/1/19

DataLoader <- function(file_dir = getwd()){
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
  
  return(list("target_rnaseq"   = target_rnaseq, 
              "target_clinical" = target_clinical, 
              "tcga_rnaseq"     = tcga_rnaseq, 
              "tcga_clinical"   = tcga_clinical))
}