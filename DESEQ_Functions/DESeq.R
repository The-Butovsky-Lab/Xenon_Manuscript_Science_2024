###########################################
### Purpose of the code is to perform DESEQ
### Author: Madison Carpenter
###########################################

##################### Setup the Environment ####################################
library(tximport)
library(DESeq2)
library(biomaRt)
library(magrittr)
library(tidyverse)
library(WriteXLS)


################ Tx2gene function #####################################

tx2genefun <- function(reference) {
  # Determine Reference
  if(toupper(reference) == "HUMAN") {
    ref_type <- "hsapiens_gene_ensembl"
  } else{
    ref_type <- "mmusculus_gene_ensembl"}
  
  # useMart is a function used to connect to the selected BioMart database and dataset (ListEnsemblArchives() to list host sites)
  ensembl_hs_mart <- useMart(biomart="ensembl", dataset= ref_type)
  print("step1")
  
  ## getBM is a function used to retrieve information from the BioMart database
  ensembl_df <- getBM(attributes = c("ensembl_transcript_id_version", "ensembl_gene_id_version", "external_gene_name"), 
                      mart = ensembl_hs_mart) %>%
    as_tibble() %>%
    mutate(noVersion = sub('\\..*', '', ensembl_transcript_id_version)) %>% # Remove Version off transcript ID
    mutate(noVersion_ENSG = sub('\\..*', '', ensembl_gene_id_version))
  
  # Build tx2gene table from ensembl_df
  tx2gene <- ensembl_df %>%
    dplyr::select(noVersion, external_gene_name) %>%
    dplyr::rename(c('Transcript' = "noVersion" , "Gene" = "external_gene_name")) %>%
    data.frame()
  
  # Return the Tx2gene database
  return(tx2gene)
}


############################  TXIMPORT #########################################

# Setup Tximport Function
Tximportfun <- function(samplefiles, reference, metadata, design, tx2gene) {
  # Set status
  print("Running Tximport")

  # Run Tximport
  txi<-tximport(files=samplefiles,type="salmon", tx2gene = tx2gene, ignoreTxVersion = T)

  print("Creating DESeq object")
  
  # Setting up DESEQ Object
  design = noquote(paste0("~ ", noquote(design))) %>% as.formula()
  ddsTxi <- DESeqDataSetFromTximport(txi, colData = metadata, design = design)
  
  return(ddsTxi)
}

############################ DESEQ  ############################################

# Construct a function to run DESeq
DESeqFun <- function(treatcolumn, control, DESeqObject, var1 = NA, var2 = NA, test_type) {
  # Setup Variables
  DESeqObjects <- list()
  
  # Remove Quotes & Assign the control
  treatcolumns <- noquote(treatcolumn)
  DESeqObject@colData@listData[[treatcolumns]] <- relevel(DESeqObject@colData@listData[[treatcolumns]], ref = control) 
  
  if(toupper(test_type) == "WALD") {
    print("Running DESEQ with WALD")
    dds <- DESeq(DESeqObject, betaPrior = F)
    dds_result <- results(dds, contrast = c(treatcolumn,var1, var2), independentFiltering = T, pAdjustMethod = 'BH')
    dds_result <- lfcShrink(dds, contrast = c(treatcolumn, var1, var2), res = dds_result, type = 'normal')
  } else{
    print("Running DESEQ with LRT")
    dds <- DESeq(DESeqObject, test = 'LRT', reduced = ~ 1)
    dds_result <- results(dds, pAdjustMethod = 'BH', independentFiltering = T, cooksCutoff = T)
  }

  # Return Deseq object
  DeseqObjects <- list(dds, dds_result)
  return(DeseqObjects)
}







