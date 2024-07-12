# packages
#install.packages('BiocManager')
library(BiocManager)

library(TCGAbiolinks)
install.packages('magrittr')
library(magrittr)
#install('org.Hs.eg.db')
library("org.Hs.eg.db", character.only = TRUE)
library(SummarizedExperiment)


# UTF-8

# downloading data
luad_query <- GDCquery(project = "TCGA-LUAD",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification"
)
GDCdownload(luad_query)

# check count of sample types
as.vector(as.data.frame(luad_query[[1]][[1]])$sample_type) %>%
  table()

# data preparing
luad_data <- GDCprepare(luad_query)
luad_expMat2 <- assay(luad_data)

luad_expMat2 <- read.csv('C:/Users/user/Desktop/Учёба/2 курс/Практика/Датасет/luad_expMat.csv')

# save the matrix
write.csv(luad_expMat2, 'luad_expMat.csv')

# Attantion!!!
# luad_expMat has not normalized data!

# converting gene names
genes_ensg <- gsub("\\..*","", rownames(luad_expMat2))
genes_symb <- mapIds(org.Hs.eg.db, keys = genes_ensg,
                     keytype ="ENSEMBL", column ="SYMBOL") %>%
  as.vector()

# our set of the genes
set_gen <- c("STAT2", "STAT1", "JAK2", "STAT5A", "JAK1", "SRC", "E2F8", "E2F7", "E2F6", "E2F5", "E2F4", "E2F2", "E2F1", "CDKN2A",
            "CDKN2B", "CDK4", "MYC", "CCND1", "CDKN1B", "CCND2", "CDK6", "CDKN1A", "CDC25A", "CDK2", "CCNE1", "STAT3", "STAT5B",
            "E2F3", "CDK1", "CCNB1", "CCNA1", "RBL2", "RBL1", "RB1")


# to take information only about genes from our set
luad_expMat2 %>%
  as.data.frame() %>%
  .[genes_symb %in% genes_ensg, ] %>%
  na.omit()

# luad_expMat2 %>%
#   as.data.frame() %>%
#   .[genes_symb %in% gene_names, ] %>%
#   na.omit()
# Let's start!

