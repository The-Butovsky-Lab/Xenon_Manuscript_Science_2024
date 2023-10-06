#################################
## Export Stats on DESEQ Analysis
## Updated: February 11, 2023
#################################



##################### Setup the Environment ####################################
library(tximport)
library(DESeq2)
library(WriteXLS)

##################### Extract Counts & Top Genes ###############################

# Construct function to extract stats & counts
extractstats <- function(DeseqObj,ContrastObj,graphlabel) {
  
  # Setup up Variables
  graphlabel <- c("gene", graphlabel)
  
  # Remove Unwatend Genes
  unwanted_genes = paste(c('^Gm', '^mt-', '^Vmn', '^Rpl', '^Rps', '^Olfr','Rik'), collapse = '|')
  
  ## Save DESEq counts & results
  DESeq_res <- ContrastObj %>% as.data.frame() %>% rownames_to_column(., var = "gene") %>% filter(!str_detect(gene, unwanted_genes))
   DESeqcounts <- DeseqObj %>% counts(normalized = TRUE) %>% as.data.frame() %>% rownames_to_column(., var = "gene") %>% filter(!str_detect(gene, unwanted_genes))
  names(DESeqcounts) <- graphlabel
  
  # Save significant results 
  sigDESeq_res <- DESeq_res[DESeq_res$pvalue<0.05 & !is.na(DESeq_res$pvalue),] %>% as.data.frame() %>% arrange(pvalue) # Generate Sig Results
  sig_counts <- DESeqcounts %>% .[.$gene %in% sigDESeq_res$gene,] %>% format(.,scientific=FALSE)
  adj_order <- match(sigDESeq_res$gene, sig_counts$gene)
  ordercounts  <- sig_counts[adj_order, ]
  rownames(ordercounts) <- NULL
  
  # Convert to numeric
  ordercounts <- ordercounts %>% column_to_rownames(., var = "gene") %>% 
  mutate_all(., function(x) as.numeric(as.character(x))) %>% rownames_to_column(., var = "gene")
  
  # Return Objects
  extractedobjects <- list(DESeq_res,DESeqcounts,sigDESeq_res,ordercounts)
  return(extractedobjects)
}

# Construct function to export stats & counts
exportstats <- function(DEResults, sigDEResults, DEcounts,sigDEcounts,directory, fileprefix) {
  # Export sig counts and genes as files
  WriteXLS::WriteXLS(DEResults, ExcelFileName = paste0(directory,fileprefix, '_statistics.xlsx')) # Stats
  WriteXLS::WriteXLS(sigDEResults, ExcelFileName = paste0(directory,fileprefix, '_statistics_pval05.xlsx'))
  
  # Export sig counts and genes as files
  WriteXLS::WriteXLS(DEcounts, ExcelFileName = paste0(directory,fileprefix, '_counts.xlsx')) # counts
  WriteXLS::WriteXLS(sigDEcounts, ExcelFileName = paste0(directory,fileprefix, '_counts_pval05.xlsx'))
}

