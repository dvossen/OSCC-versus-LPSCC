rm( list = ls( all = TRUE ) )

### Somatic mutations
# NKI OSCC
load("~/Documents/Projects/Analyses_SNV/VarScan/Somatic_prediction/Mondholte_Somatic_mutations.Rdata")
idx.2e.primaire = which(rownames(binary.matrix.mondholte.variants.somatic) == "D032-KOb")
idx.adenocarci = which(rownames(binary.matrix.mondholte.variants.somatic) == "D016-EE1")
nki.oscc.mut = binary.matrix.mondholte.variants.somatic[-c(idx.2e.primaire, idx.adenocarci),]

# NKI L/P-SCC
source("~/Documents/Projects/Analyses_SNV/MSigDB/RadPlat_clinical_data.R")
load("~/Documents/Projects/Analyses_SNV/VarScan/Somatic_prediction/RadPlat_Somatic_mutations_minfreq_10.Rdata")
radplat.hpv.neg.samples = clinical.data.radplat.excl.minor.sites$DNA_lab[which(!clinical.data.radplat.excl.minor.sites$HPV.status.DNAseq)]
radplat.hpv.neg.somatic.mutations = binary.matrix.radplat.variants.somatic[which(rownames(binary.matrix.radplat.variants.somatic) %in% radplat.hpv.neg.samples),]
nki.lpscc.mut = radplat.hpv.neg.somatic.mutations


### Copy number alterations
# NKI OSCC
load("~/Documents/Projects/Analyses_CNA/PureCN/Mondholte/CN_Matrix_PureCN_Mondholte.Rdata")
nki.oscc.cn = oscc.cn.matrix
idx.2e.primaire = which(rownames(nki.oscc.cn) == "D032-KOb1")
idx.adenocarci = which(rownames(nki.oscc.cn) == "D016-EE1")
nki.oscc.cn = nki.oscc.cn[-c(idx.2e.primaire, idx.adenocarci),]
# Change one rowname
rownames(nki.oscc.cn)[which(rownames(nki.oscc.cn) == "D032-KOa1")] = "D032-KOa"
# Check identical dimnames
identical(dimnames(nki.oscc.mut), dimnames(nki.oscc.cn))
# Make gain and loss matrices
nki.oscc.gain = ifelse(nki.oscc.cn >= 7, 1, 0)
nki.oscc.loss = ifelse(nki.oscc.cn <= 0.5, 1, 0)

# NKI L/P-SCC
source("~/Documents/Projects/Analyses_SNV/MSigDB/RadPlat_clinical_data.R")
load("~/Documents/Projects/Analyses_CNA/PureCN/RadPlat/CN_Matrix_PureCN_RadPlat.Rdata")
radplat.hpv.neg.samples = clinical.data.radplat.excl.minor.sites$DNA_lab[which(!clinical.data.radplat.excl.minor.sites$HPV.status.DNAseq)]
radplat.hpv.neg.cn = radplat.cn.matrix[which(rownames(radplat.cn.matrix) %in% radplat.hpv.neg.samples),]
# Match rownames to mutations
nki.lpscc.cn = radplat.hpv.neg.cn[pmatch(rownames(nki.lpscc.mut), rownames(radplat.hpv.neg.cn)),]
# Check identical dimnames
identical(dimnames(nki.lpscc.mut), dimnames(nki.lpscc.cn))

nki.lpscc.gain = ifelse(nki.lpscc.cn >= 7, 1, 0)
nki.lpscc.loss = ifelse(nki.lpscc.cn <= 0.5, 1, 0)

save(nki.oscc.mut, nki.oscc.cn, nki.oscc.gain, nki.oscc.loss,
     nki.lpscc.mut, nki.lpscc.cn, nki.lpscc.gain, nki.lpscc.loss,
     file = "~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Data/NKI_Mut_and_CNA.Rdata")