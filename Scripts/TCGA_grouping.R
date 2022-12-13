# LOAD TCGA ####
# Load M1K (TCGA) mutations ####
MK = data.table::fread("E:/Users/Claris_Gu/ComputationalBiologyProjects/GeneExpSig/Data/TCGA_SKCM/melanoma_combined_variants_public_only_V1.tsv.gz")
MK = MK[MK$study == "TCGA-SKCM",]
MK$mutAA = paste0(MK$Gene_Name,".",MK$HGVS.p)
MK$Annotation1 = gsub("&.*$","", MK$Annotation)
MK$VAF = (MK$alt_count/(MK$ref_count+MK$alt_count))

# Load clinical data ####
clin.file.skcm = "E:/Users/Claris_Gu/ComputationalBiologyProjects/GeneExpSig/Data/TCGA_SKCM/clinical_new_supp_tables_232_2_supp_20255_q8v87y.txt"
CK = data.table::fread(clin.file.skcm, skip = 2)

# Subset clinical
CK = CK[CK$study == "TCGA" & CK$MAF_sample_id %in% unique(MK$sample),]

# START ###
ann.LoF = c("stop_gained","start_lost","splice_acceptor_variant","splice_donor_variant","frameshift_variant")

patt.BRAF = "BRAF.p.V600[A-Z]"
patt.NRAS = "NRAS.p.(Q61|G12|G13)[A-Z]"
CK$NF1.Lof = CK$MAF_sample_id %in% MK$sample[MK$Gene_Name == "NF1" & MK$Annotation1 %in% ann.LoF]
CK$BRAF = CK$MAF_sample_id %in% MK$sample[grepl(patt.BRAF, MK$mutAA)]
CK$NRAS = CK$MAF_sample_id %in% MK$sample[grepl(patt.NRAS, MK$mutAA)]

# Add group ####
CK$Group = "3WT"
CK$Group[CK$BRAF] = "BRAF"
CK$Group[CK$NRAS] = "NRAS"
CK$Group[CK$NF1.Lof] = "NF1.Lof"
CK$subtype = CK$Group

