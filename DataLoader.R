# TITLE: DataLoader.R
# AUTHOR: Kevin O'Connor
# DATE CREATED: 3/7/19
# DATE MODIFIED: 4/9/19

DataLoader <- function(file_dir = getwd(), 
                       normalize_quantiles = TRUE, 
                       select_primary = TRUE,
                       filt = TRUE,
                       min_med = 5,
                       min_sd = 1){
  require(preprocessCore)
  
  stringsAsFactors <- FALSE
  
  # Read all the data and set row and column names.
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
  
  tcga_as <- read.csv(file.path(file_dir, "tcga_laml_aneuploidy_040119.csv"), header=TRUE)
  tcga_as_samples <- sapply(tcga_as[, 1], function(x) substr(x, start=1, stop=12))
    
  # Merge rnaseq data.
  rnaseq <- merge(tcga_rnaseq, target_rnaseq, by="row.names")
  rownames(rnaseq) <- rnaseq$Row.names
  rnaseq <- rnaseq[, -1]
  sample_flag <- c(rep("tcga", ncol(tcga_rnaseq)), rep("target", ncol(target_rnaseq)))
  
  # Quantile normalization.
  if(normalize_quantiles){
    rnaseq_unnorm <- rnaseq
    rnaseq <- as.matrix(rnaseq) %>% normalize.quantiles() %>% as.data.frame()
    rownames(rnaseq) <- rownames(rnaseq_unnorm)
    tcga_colnames <- colnames(tcga_rnaseq) %>%
      gsub(pattern="[.]", replacement="-") %>% 
      substr(start=1, stop=12)
    target_colnames <- colnames(target_rnaseq) %>% 
      gsub(pattern="[.]", replacement="-") %>% 
      substr(start=1, stop=16)
    colnames(rnaseq) <- c(tcga_colnames, target_colnames)
  }

  # Select primary samples in target.
  if(select_primary){
    good_samples <- c(1:ncol(tcga_rnaseq), ncol(tcga_rnaseq) + which(substr(colnames(target_rnaseq), start=18, stop=20) == "09A"))
    rnaseq <- rnaseq[, good_samples]  
    sample_flag <- sample_flag[good_samples]
  } else {
    good_samples <- 1:ncol(rnaseq)
  }
 
  # Filter.
  rnaseq_unfilt <- rnaseq
  if(filt){
    good_genes <- which(matrixStats::rowMedians(as.matrix(rnaseq)) > min_med)
    good_genes <- intersect(good_genes, which(genefilter::rowSds(rnaseq) > min_sd))
    rnaseq <- rnaseq[good_genes, ]
  } else {
    good_genes <- 1:nrow(rnaseq)
  }
  
  # Log2 transform.
  rnaseq <- log(rnaseq + 1, base = 2)
  
  # Standardize and median center.
  rnaseq <- rnaseq / genefilter::rowSds(rnaseq)
  rnaseq <- rnaseq - matrixStats::rowMedians(as.matrix(rnaseq))
  
  # Gather clinical data.
  tcga_clin_inds <- match(colnames(rnaseq)[which(sample_flag == "tcga")], rownames(tcga_clinical))
  target_clin_inds <- match(colnames(rnaseq)[which(sample_flag == "target")], rownames(target_clinical))
  # Get survival times for each patient.
  tcga_survtime <- tcga_clinical[tcga_clin_inds, 10]
  target_survtime <- target_clinical[target_clin_inds, 8]
  survtime <- c(tcga_survtime, target_survtime)
  # Genders
  tcga_gender <- tcga_clinical$"gender"[tcga_clin_inds] %>% as.character() %>% sapply(tolower)
  target_gender <- target_clinical$"Gender"[target_clin_inds] %>% as.character() %>% sapply(tolower)
  gender <- c(tcga_gender, target_gender) %>% as.character()
  # Race
  tcga_race <- tcga_clinical$"race_list"[tcga_clin_inds] %>% as.character() %>% sapply(tolower)
  target_race <- target_clinical$"Race"[target_clin_inds] %>% as.character() %>% sapply(tolower)
  race <- c(tcga_race, target_race) %>% as.character()
  race[which(race == "")] <- "unknown"
  # Age
  tcga_age <- tcga_clinical$"days_to_birth"[tcga_clin_inds]
  target_age <- target_clinical$"Age.at.Diagnosis.in.Days"[target_clin_inds]
  age <- c(-1*tcga_age, target_age)
  # Aneuploidy
  tcga_as <- dplyr::select(tcga_as, Aneuploidy.Score)
  rownames(tcga_as) <- tcga_as_samples
  
  # Combine into a single dataframe.
  clinical <- data.frame(gender = gender,
                         race = race,
                         age = age,
                         survtime = survtime)
  rownames(clinical) <- colnames(rnaseq)
  clinical <- merge(clinical, tcga_as, by="row.names", all.x = TRUE)
   
  return(list("rnaseq" = rnaseq, "clinical" = clinical))
}
