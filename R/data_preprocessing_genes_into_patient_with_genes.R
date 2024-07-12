# packages
#install.packages('BiocManager')
library(BiocManager)
#install('TCGAbiolinks')
library(TCGAbiolinks)
#install.packages('magrittr')
library(magrittr)
#install('org.Hs.eg.db')
library("org.Hs.eg.db", character.only = TRUE)


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
luad_expMat <- assay(luad_data)

# save the matrix
write.csv(luad_expMat, 'luad_expMat.csv')

# Attantion!!!
# luad_expMat has not normalized data!

# converting gene names
genes_ensg <- gsub("\\..*","", rownames(luad_expMat))
genes_symb <- mapIds(org.Hs.eg.db, keys = genes_ensg,
                     keytype ="ENSEMBL", column ="SYMBOL") %>%
  as.vector()

# our set of the genes
set_gen <- c("STAT2", "STAT1", "JAK2", "STAT5A", "JAK1", "SRC", "E2F8", "E2F7", "E2F6", "E2F5", "E2F4", "E2F2", "E2F1", "CDKN2A",
            "CDKN2B", "CDK4", "MYC", "CCND1", "CDKN1B", "CCND2", "CDK6", "CDKN1A", "CDC25A", "CDK2", "CCNE1", "STAT3", "STAT5B",
            "E2F3", "CDK1", "CCNB1", "CCNA1", "RBL2", "RBL1", "RB1")


#My code

# to take information only about genes from our set
result <- luad_expMat %>%
  as.data.frame() %>%
  .[genes_symb %in% set_gen, ] %>%
  na.omit() %>% t()

# write.csv(result, 'genes_without_T.csv')

extract_feautures <- function(sample_id, column) {
  parts <- strsplit(sample_id, "-")[[1]]
  if (column == "id"){
    res <- paste(parts[1], parts[2], parts[3], sep = "-")
  }else if (column == "type_sample"){
    res <- substr(parts[4],1,2)
  }else if (column == "vial"){
    res <- substr(parts[4],3,3)
  }else if (column == "portion_mg"){
    res <- substr(parts[5],1,2)
  }else if (column == "type_analyte"){
    res <- substr(parts[5],3,3)
  }else if (column == "plate"){
    res <- parts[6]
  }else if (column == "centre"){
    res <- parts[7]
  }
  return (res)
}

result <- transform(result,id=mapply(extract_feautures,row.names(result),"id"))
result <- transform(result,type_sample=mapply(extract_feautures,row.names(result),"type_sample"))
result <- transform(result,vial=mapply(extract_feautures,row.names(result),"vial"))
result <- transform(result,portion_mg=mapply(extract_feautures,row.names(result),"portion_mg"))
result <- transform(result,type_analyte=mapply(extract_feautures,row.names(result),"type_analyte"))
result <- transform(result,plate=mapply(extract_feautures,row.names(result),"plate"))
result <- transform(result,centre=mapply(extract_feautures,row.names(result),"centre"))

#information about patients
clinical_data <- GDCquery_clinic(project = "TCGA-LUAD", type = "clinical")

result <- merge(clinical_data, result, by.x = "submitter_id", by.y = "id")

print(colnames(result))

write.csv(result, 'patients_with_genes.csv')
