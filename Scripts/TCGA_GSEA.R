TCGA.TPM.file = "E:/Users/Claris_Gu/ComputationalBiologyProjects/GeneExpSig/Data/TCGA_SKCM/KALLISTO_V42_table_SKCM.tsv"
TPM = data.table::fread(TCGA.TPM.file)
selcols = colnames(TPM)[-1]
X = as.matrix(TPM.genes[,..selcols])
rownames(X) = TPM.genes$gene.name

# get samples for groups: BRAF, NRAS, NF1
BRAF.samples = CK$clinical_sample_id[CK$BRAF]
BRAF.samples = paste0(BRAF.samples, "A")
BRAF_cols = colnames(X)[colnames(X) %in% BRAF.samples]
BRAF.samples = X[,BRAF_cols]

NRAS.samples = CK$clinical_sample_id[CK$NRAS]
NRAS.samples = paste0(NRAS.samples, "A")
NRAS.cols = colnames(X)[colnames(X) %in% NRAS.samples]
NRAS.samples = X[,NRAS.cols]

NF1.Lof.samples = CK$clinical_sample_id[CK$NF1.Lof]
NF1.Lof.samples = paste0(NF1.Lof.samples, "A")
NF1.Lof.cols = colnames(X)[colnames(X) %in% NF1.Lof.samples]
NF1.Lof.samples = X[,NF1.Lof.cols]

# Run GSEA
library(xCell)
BRAF.ex = xCell::xCellAnalysis(BRAF.samples, rnaseq=TRUE)
NRAS.ex = xCell::xCellAnalysis(NRAS.samples, rnaseq=TRUE)
NF1.Lof.ex = xCell::xCellAnalysis(NF1.Lof.samples, rnaseq=TRUE)