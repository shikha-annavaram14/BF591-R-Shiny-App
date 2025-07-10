if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Dr.eg.db")
library(org.Dr.eg.db)
library(tidyverse)
library(fgsea)
library(msigdbr)

counts <- read.csv("deseq_normalized_counts.csv") %>% as.data.frame()

id2gene <- mapIds(org.Dr.eg.db, keys = counts$gene,
       column = c('SYMBOL'), keytype = 'ENSEMBL') %>% as.data.frame() %>% rownames_to_column(var = 'gene') %>% dplyr::rename(symbol = 2) %>% drop_na()


de_results <- read.csv("differential_expr_results.csv") %>% as.data.frame() %>% dplyr::select(gene = 1, everything())


joined <- left_join(de_results, id2gene, by = 'gene') %>% dplyr::filter(!is.na(symbol)) %>% 
   dplyr::filter(!is.na(log2FoldChange)) %>% arrange(desc(log2FoldChange))

ranked_vector <- setNames(joined$log2FoldChange, joined$symbol)

gene_sets <- msigdbr(species = "zebrafish")
msigdbr_list <- split(x = gene_sets$gene_symbol, f = gene_sets$gs_name)
fgsea <- fgsea(pathways = msigdbr_list, ranked_vector, minSize = 15, maxSize = 500)
data.table::fwrite(fgsea, file="fgsea.tsv", sep="\t", sep2=c("", " ", ""))

fgsea <- read_tsv("fgsea.tsv") %>% as.data.frame()

