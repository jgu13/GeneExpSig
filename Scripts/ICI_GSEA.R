ICI.TPM.file = "E:/Users/Claris_Gu/ComputationalBiologyProjects/GeneExpSig/Data/ICI/kallisto_ICI_V42.tsv.gz"
ICI.TPM = data.table::fread(ICI.TPM.file)
selcols = colnames(ICI.TPM)[-2:-1]
ICI.X = as.matrix(ICI.TPM[,..selcols])
rownames(ICI.X) = ICI.TPM$gene.name

# get samples for groups: BRAF, NRAS, NF1
ICI.BRAF.samples = ICI.CK$SRR.RNA[ICI.CK$BRAF]
ICI.BRAF.cols = colnames(ICI.X)[colnames(ICI.X) %in% ICI.BRAF.samples]
ICI.BRAF.samples = ICI.X[,ICI.BRAF.cols]

ICI.NRAS.samples = ICI.CK$SRR.RNA[ICI.CK$NRAS]
ICI.NRAS.cols = colnames(ICI.X)[colnames(ICI.X) %in% ICI.NRAS.samples]
ICI.NRAS.samples = ICI.X[,ICI.NRAS.cols]

ICI.NF1.Lof.samples = ICI.CK$SRR.RNA[ICI.CK$NF1.Lof]
ICI.NF1.Lof.cols = colnames(ICI.X)[colnames(ICI.X) %in% ICI.NF1.Lof.samples]
ICI.NF1.Lof.samples = ICI.X[,ICI.NF1.Lof.cols]

# Run GSEA
library(xCell)
ICI.BRAF.ex = xCell::xCellAnalysis(ICI.BRAF.samples, rnaseq=TRUE)
ICI.NRAS.ex = xCell::xCellAnalysis(ICI.NRAS.samples, rnaseq=TRUE)
ICI.NF1.Lof.ex = xCell::xCellAnalysis(ICI.NF1.Lof.samples, rnaseq=TRUE)
