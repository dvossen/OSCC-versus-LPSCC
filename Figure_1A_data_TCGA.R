rm( list = ls( all = TRUE ) )

### Somatic mutations
# TCGA OSCC
load("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Data/TCGA_Mut_and_CNA_small.Rdata")
tcga.oscc.mut = tcga.mut.matrix[which(rownames(tcga.mut.matrix) %in% tcga.clinical.data.hpv.neg$bcr_patient_barcode[which(tcga.clinical.data.hpv.neg$Location == "Oral Cavity")]),]
# TCGA L/P-SCC
tcga.lpscc.mut = tcga.mut.matrix[which(rownames(tcga.mut.matrix) %in% tcga.clinical.data.hpv.neg$bcr_patient_barcode[which(tcga.clinical.data.hpv.neg$Location != "Oral Cavity")]),]
identical(colnames(tcga.oscc.mut), colnames(tcga.lpscc.mut))


### Copy number alterations
# TCGA OSCC
load("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Data/TCGA_Mut_and_CNA_small.Rdata")
tcga.oscc.cn = tcga.cna[which(rownames(tcga.cna) %in% tcga.clinical.data.hpv.neg$bcr_patient_barcode[which(tcga.clinical.data.hpv.neg$Location == "Oral Cavity")]),]
tcga.oscc.gain = ifelse(tcga.oscc.cn == 2, 1, 0)
tcga.oscc.loss = ifelse(tcga.oscc.cn == -2, 1, 0)
# TCGA L/P-SCC
tcga.lpscc.cn = tcga.cna[which(rownames(tcga.cna) %in% tcga.clinical.data.hpv.neg$bcr_patient_barcode[which(tcga.clinical.data.hpv.neg$Location != "Oral Cavity")]),]
tcga.lpscc.gain = ifelse(tcga.lpscc.cn == 2, 1, 0)
tcga.lpscc.loss = ifelse(tcga.lpscc.cn == -2, 1, 0)
identical(colnames(tcga.oscc.cn), colnames(tcga.lpscc.cn))

identical(dimnames(tcga.oscc.mut), dimnames(tcga.oscc.cn))
identical(dimnames(tcga.lpscc.mut), dimnames(tcga.lpscc.cn))

save(tcga.oscc.mut, tcga.oscc.gain, tcga.oscc.loss,
     tcga.lpscc.mut, tcga.lpscc.gain, tcga.lpscc.loss,
     file = "~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Data/Figure_1A_data_TCGA.Rdata")