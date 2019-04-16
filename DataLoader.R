# TITLE: DataLoader.R
# AUTHOR: Kevin O'Connor
# DATE CREATED: 3/7/19
# DATE MODIFIED: 4/14/19

DataLoader <- function(file_dir = getwd(), 
                       normalize_quantiles = TRUE, 
                       filt = TRUE,
                       min_med = 5,
                       min_sd = 1.5){
  stringsAsFactors <- FALSE
  
  # Read all the data and set row and column names.
  rnaseq <- read.csv(file.path(file_dir, "tcga_laml_rnaseq_040119.csv"), header = TRUE)
  rownames(rnaseq) <- rnaseq[, 1]
  rnaseq <- rnaseq[, -1]
  
  clinical <- read.csv(file.path(file_dir, "tcga_laml_clinical_040119.csv"), header = TRUE)
  rownames(clinical) <- clinical[, 1]
  clinical <- clinical[, -1]
  
  tcga_as <- read.csv(file.path(file_dir, "tcga_laml_aneuploidy_040119.csv"), header=TRUE)
  tcga_as_samples <- sapply(tcga_as[, 1], function(x) substr(x, start=1, stop=12))
    
  # Quantile normalization.
  ## Quartile normalization function courtesy of Joel Parker.
  quartileNorm<- function(x, y = NA){
    uqs <- apply(x, 2, function(x){
      quantile(x[x>0 & !is.na(x)], 0.75)
    })
    if(is.na(y)){
      y <- median(uqs)
    }
    x.norm <- t(apply(x, 1, function(x,y){x * y}, y / uqs))
    dimnames(x.norm) <- dimnames(x)
    return(x.norm)
  }
  if(normalize_quantiles){
    rnaseq_unnorm <- rnaseq
    rnaseq <- quartileNorm(rnaseq)
    rownames(rnaseq) <- rownames(rnaseq_unnorm)
    colnames(rnaseq) <- colnames(rnaseq_unnorm) %>%
      gsub(pattern="[.]", replacement="-") %>% 
      substr(start=1, stop=12)
  }

  # Log2 transform.
  rnaseq <- log(rnaseq + 1, base = 2)
   
  # Filter.
  rnaseq_unfilt <- rnaseq
  if(filt){
    good_genes <- which(matrixStats::rowMedians(as.matrix(rnaseq)) > min_med)
    good_genes <- intersect(good_genes, which(genefilter::rowSds(rnaseq) > min_sd))
    rnaseq <- rnaseq[good_genes, ]
  } else {
    good_genes <- 1:nrow(rnaseq)
  }
  
  # Standardize and median center.
  rnaseq <- rnaseq / genefilter::rowSds(rnaseq)
  rnaseq <- rnaseq - matrixStats::rowMedians(as.matrix(rnaseq))
  
  # Gather clinical data.
  clin_inds <- match(colnames(rnaseq), rownames(clinical))
  # Get survival times for each patient.
  survtime <- clinical[clin_inds, 10]
  # Genders
  gender <- clinical$"gender"[clin_inds] %>% as.character() %>% sapply(tolower)
  # Race
  race <- clinical$"race_list"[clin_inds] %>% as.character() %>% sapply(tolower)
  race[which(race == "")] <- "unknown"
  # Age
  age <- -1*clinical$"days_to_birth"[clin_inds]
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
