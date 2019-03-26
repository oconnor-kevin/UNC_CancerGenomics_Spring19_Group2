# TITLE: DataLoader.R
# AUTHOR: Kevin O'Connor
# DATE CREATED: 3/7/19
# DATE MODIFIED: 3/12/19

DataLoader <- function(write.rnaseq = FALSE){
  # Download data.
  ## RNAseq data.
  rna.query <- TCGAbiolinks::GDCquery(project = "TARGET-AML",
                                      data.category = "Transcriptome Profiling",
                                      data.type = "Gene Expression Quantification",
                                      workflow.type = "HTSeq - Counts")
  TCGAbiolinks::GDCdownload(rna.query, method = "api", files.per.chunk = 10)
  ## Clinical data.
  clinical.query <- TCGAbiolinks::GDCquery(project = "TARGET-AML", 
                                           data.category = "Clinical")
  TCGAbiolinks::GDCdownload(clinical.query)
  
  
  # Prepare RNA Data for Analysis.
  ## Prepares a 'Summarized Experiment' object.
  rna.seqdat <- TCGAbiolinks::GDCprepare(query = rna.query, 
                                         save = TRUE, 
                                         save.filename = "rnaseq.rda")
  ## Create a matrix of the RNA-seq data with genes as rows and samples as columns
  full.rna.seq.mat <- SummarizedExperiment::assay(rna.seqdat)
  rownames(full.rna.seq.mat) <- SummarizedExperiment::rowData(rna.seqdat)$external_gene_name
  if (write.rnaseq){
    ## Write the RNA-Seq data to a matrix
    write.csv(full.rna.seq.mat, file = "full_rna_seq_mat.csv",
              row.names = TRUE) 
  }
  
  return(full.rna.seq.mat)
}