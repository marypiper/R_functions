# List of functions for RNA-Seq analyses

## GSEA analysis: conversion of Entrez IDs to gene symbols (using output from clusterProfiler)

```r

## Function to parse the associated Entrez IDs for each pathway and split at the "/", then extract these genes from the database of IDs and output the gene symbol
gsea_entrez_to_symbol <- function(list_idx){
path <- gseaKEGGSummary$core_enrichment %>% str_split(pattern="/") %>% .[[list_idx]]


idx <- which(annotations$entrez %in% path)

annotations[idx,  "symbol"]
}

# Intialize list to output symbols
list_gsea_symbols <- list()

# Output symbols for every pathway
list_gsea_symbols <- map(1:nrow(gseaKEGGSummary), gsea_entrez_to_symbol)

# Name pathways based on the summary IDs
names(list_gsea_symbols) <- gseaKEGGSummary$ID

# Function to convert the associated gene symbols to a single value separated by "/"
paste_symbols <- function(list_idx){
  paste(list_gsea_symbols[[list_idx]], collapse ="/")
}

# Collapse to a single value for every pathway
symbols <- map(1:length(list_gsea_symbols), paste_symbols)
names(symbols) <- names(list_gsea_symbols)

# Turn the single values into a vector
symbols <- unlist(symbols)

# Transpose the vector to be a column of symbols
gene_symbols <- t(data.frame(symbols))

# Add column to kegg summary
gseaKEGGSummaryGenes <- data.frame(gseaKEGGSummary, gene_symbols)
```

## Changing *C. elegans* sequence IDs in Wormcat annotations to symbols

```r
# Read in wormcat annotations downloaded from Wormcat (http://wormcat.com)
wormcat <- read_csv("wormcat_whole_genome_nov-16-2019.csv")

# Connect to AnnotationHub
ah <- AnnotationHub()

# Query AnnotationHub
cel_ens <- query(ah, c("Caenorhabditis elegans", "EnsDb"))


# Get most recent one
cel_ens <- cel_ens[["AH78732"]]

# Extract gene-level information
genedb <- genes(cel_ens, return.type = "data.frame") %>% 
  dplyr::select(gene_id, gene_name, symbol, description)  

annotations <- genedb

annotations$gene_name <-  str_replace(string=genedb$gene_name, pattern = "\\.", replacement = "_")

wormcat_merged <- left_join(wormcat, annotations, by = c("Wormbase ID" = "gene_id"))

## Function to parse the associated IDs for each pathway and split at the "/", then extract these genes from the database of IDs  and output the gene symbol
ora_seqID_to_symbol <- function(list_idx, cat_enrich_obj) {
  path <- cat_enrich_obj$geneID %>%
    str_split(pattern="/") %>% .[[list_idx]]
  idx <- which(wormcat_merged$`Sequence ID` %in% path)
  wormcat_merged[idx,  "symbol"]
}

# Intialize list to output symbols
get_symbols <- function(cat_enrich_obj){
  list_ora_symbols <- list()
  # Output symbols for every pathway
  list_ora_symbols <-
    map(1:nrow(cat_enrich_obj),
        ora_seqID_to_symbol, cat_enrich_obj)
  # Name pathways based on the summary IDs
  names(list_ora_symbols) <-
    cat_enrich_obj$ID
  # Function to convert the associated gene symbols to a single value separated by "/"
  paste_symbols <- function(list_idx){
    paste(list_ora_symbols[[list_idx]]$symbol, collapse ="/")
  }
  
  # Collapse to a single value for every pathway
  symbols <- map(1:length(list_ora_symbols), paste_symbols)
  names(symbols) <- names(list_ora_symbols)
  # Turn the single values into a vector
  symbols <- unlist(symbols)
  
  # Transpose the vector to be a column of symbols
  gene_symbols <- data.frame(symbols)
  # Add column to kegg summary
  cat_enrich <- data.frame(cat_enrich_obj, gene_symbols)
  return(cat_enrich)
}

cat1_sig_enrich <- get_symbols(cat1_sig_enrich)
```
