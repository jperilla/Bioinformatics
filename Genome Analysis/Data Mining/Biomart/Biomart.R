source("https://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")
y
biocLite("biomaRt")
a
library(biomaRt)
listEnsembl()
listDatasets(mart=ensembl)
listMarts()
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
head(listDatasets(ensembl))


myentrez <- c("3043", "3045", "3046", "3047", "4151")
listAttributes(ensembl)
data <- getBM(attributes=c("entrezgene", "hgnc_symbol", "percentage_gene_gc_content"), 
                filters="entrezgene", values=myentrez, mart=ensembl)
 
data

affyids <- c("202763_at", "209310_s_at", "207500_at")
data2<- getBM(attributes=c("affy_hg_u133_plus_2", "hgnc_symbol", "chromosome_name",
                   "start_position", "end_position", 'band'), 
      filters="affy_hg_u133_plus_2", values=affyids, mart=ensembl)
data2

entrez <- c("673", "837")
data3 <- getBM(attributes=c("entrezgene", "go_id"), filters='entrezgene', values=entrez, mart=ensembl)
data3

go <- c("GO:0000114", "GO:0000082")
chrom <- c(17,20,"Y")
data4 <- getBM(attributes=c("hgnc_symbol"), filters=c('go', 'chromosome_name'), 
               values=list(go,chrom), mart=ensembl)
data4

refseqids <- c("NP_005350", "NP_000537")
data5 <- getBM(attributes=c("refseq_peptide", "interpro", "interpro_description"),
                filters="refseq_peptide", values=refseqids, mart=ensembl)
data5

data6 <- getBM(attributes=c("entrezgene", "hgnc_symbol"), filters="go", values="GO:0004707", mart=ensembl)
data6

data7 <- getBM(attributes="affy_hg_u133_plus_2", filters=c("chromosome_name", "start", "end"), 
               values=list(16,1100000,1250000), mart=ensembl)
data7

diabetes <- c("125853","222100")
ungradedhw <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), 
                    filters=c("mim_morbid_accession"), 
                    values=diabetes, mart=ensembl)
ungradedhw

lymphoma <- c("605027")
lymphoma_genes <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), 
                    filters=c("mim_morbid_accession"), 
                    values=lymphoma, mart=ensembl)
lymphoma_genes

casp10 <- c("ENSG00000003400")
listFilters(ensembl)
listAttributes(ensembl)
casp10_transcripts <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol","transcript_count"),
                                filters="ensembl_gene_id", values=casp10, mart=ensembl)
casp10_transcripts

huntingtons <- c("143100", "613004", "604802",  "606483", "603218")
huntingtons_table <- getBM(attributes=c("entrezgene", "hgnc_symbol", "ensembl_gene_id"),
                         filters="mim_morbid_accession", values=huntingtons, mart=ensembl)
huntingtons_table

huntingtons <- c("143100", "613004", "604802",  "606483", "603218")
huntingtons_transcripts <- getBM(attributes=c("hgnc_symbol", "ensembl_gene_id", "ensembl_transcript_id"),
                           filters="mim_morbid_accession", values=huntingtons, mart=ensembl)
huntingtons_transcripts


huntingtons_transcripts <- getBM(attributes=c("hgnc_symbol", "ensembl_gene_id", 
                                              "ensembl_transcript_id"),
                                              filters="mim_morbid_accession", 
                                              values=huntingtons, mart=ensembl)



# FUNGI
ensemblFungi <- useMart("fungi_mart", host="fungi.ensembl.org")

# Save the list of datasets available 
mart.datasets <- listDatasets(mart=ensemblFungi)
# Extract the "dataset" column which is the value to access the mart
mart.dataset <- as.character(mart.datasets$dataset)

genes_0006623_set <- data.frame(ensembl_gene_id=as.character()) 
for(dataset in mart.dataset) {
  ensemblFungi <- useMart("fungi_mart", host="fungi.ensembl.org", dataset=dataset)
  genes_0006623  <- getBM(attributes=c("ensembl_gene_id"),
                                   filters="go", 
                                   values=c("GO:0006623"), mart=ensemblFungi)
  genes_0006623_set <- rbind(genes_0006623_set, genes_0006623)
}

genes_0006623_set
