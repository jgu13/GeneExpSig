# load mutations
ICI_mutations_MNV_snpEff_SNV3_HG19_V6_file = "E:/Users/Claris_Gu/ComputationalBiologyProjects/GeneExpSig/Data/ICI/ICI_mutations_MNV_snpEff_SNV3_HG19_V6.rds"
ICI_mutations_MNV_snpEff_SNV3_HG19_V6 <- readRDS(ICI_mutations_MNV_snpEff_SNV3_HG19_V6_file)
ICI.MK = data.frame(ICI_mutations_MNV_snpEff_SNV3_HG19_V6)

# grouping by mutation
ICI.MK = ICI.MK[ICI.MK$study == "Liu2019",]

ann.LoF = c("stop_gained","start_lost","splice_acceptor_variant","splice_donor_variant","frameshift_variant")
patt.BRAF = "BRAF.p.V600[A-Z]"
patt.NRAS = "NRAS.p.(Q61|G12|G13)[A-Z]"
ICI.MK.NF1.sample = ICI.MK$sample.id[ICI.MK$Gene_Name == "NF1" & ICI.MK$Annotation1 %in% ann.LoF]
ICI.MK.NRAS.sample = ICI.MK$sample.id[grepl(patt.NRAS, ICI.MK$mutAA)]
ICI.MK.BRAF.sample = ICI.MK$sample.id[grepl(patt.BRAF, ICI.MK$mutAA)]

# Load clinical data ####
clin.file = "E:/Users/Claris_Gu/ComputationalBiologyProjects/GeneExpSig/Data/ICI/clinical_ICI_RNA_dec10_2022.tsv.gz"
ICI.CK = data.table::fread(clin.file)

# Subset clinical
ICI.CK = ICI.CK[ICI.CK$study == "Liu2019" & ICI.CK$sample.id %in% unique(ICI.MK$sample.id),]

# group clinical samples
ICI.CK$NF1.Lof = ICI.CK$sample.id %in% ICI.MK.NF1.sample
ICI.CK$BRAF = ICI.CK$sample.id %in% ICI.MK.BRAF.sample
ICI.CK$NRAS = ICI.CK$sample.id %in% ICI.MK.NRAS.sample

# Add group ####
ICI.CK$Group = "3ICI.treated"
ICI.CK$Group[ICI.CK$BRAF] = "BRAF"
ICI.CK$Group[ICI.CK$NRAS] = "NRAS"
ICI.CK$Group[ICI.CK$NF1.Lof] = "NF1.Lof"
ICI.CK$subtype = ICI.CK$Group

