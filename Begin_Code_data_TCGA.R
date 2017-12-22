rm( list = ls( all = TRUE ) )

library(Hmisc)
source("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Code/Begin_Code.R")

### Somatic mutations
# Load mutation file
tcga.mut = read.table("~/Documents/Projects/TCGA_Data/Mutation data/gdac.broadinstitute.org_HNSC-TP.MutSigNozzleReportCV.Level_4.2016012800.0.0/HNSC-TP.final_analysis_set.maf", header = T, sep = '\t', stringsAsFactors = F, quote = "", fill = T)
# Select samples from primary tumors
tcga.mut = tcga.mut[which(substring(tcga.mut$Tumor_Sample_Barcode, 14, 15) == "01"),]
# Rename samples
tcga.mut$Tumor_Sample_Barcode = substring(tcga.mut$Tumor_Sample_Barcode, 0, 12)
# Remove Silent and RNA mutations
tcga.mut.non.silent = tcga.mut[which(tcga.mut$is_coding | tcga.mut$is_del | tcga.mut$is_indel | tcga.mut$is_ins | tcga.mut$is_missense | tcga.mut$is_nonsense | tcga.mut$is_splice),]

tcga.mut.non.silent.diminished = cbind(tcga.mut.non.silent$Hugo_Symbol, tcga.mut.non.silent$Tumor_Sample_Barcode)
colnames(tcga.mut.non.silent.diminished) = c("annovar.gene", "sample")
tcga.mut.non.silent.diminished = as.data.frame(tcga.mut.non.silent.diminished)
tcga.mut.matrix = return.binary.varscan.variant.matrix(input.df.variants = tcga.mut.non.silent.diminished,
                                                       input.vector.target.genes = sort(unique(tcga.mut.non.silent.diminished$annovar.gene)),
                                                       input.vector.samples = sort(unique(tcga.mut.non.silent.diminished$sample)))



### CNAs
tcga.cna = read.table("~/Downloads/gdac.broadinstitute.org_HNSC-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/all_thresholded.by_genes.txt", header = T, sep = "\t", stringsAsFactors = F, quote = "")
# Replace . by - in sample names, to facilitate matching with mutation data
rownames(tcga.cna) = tcga.cna$Gene.Symbol
colnames(tcga.cna) = substring(gsub(pattern = ".", replacement = "-", x = colnames(tcga.cna), fixed = T), 0, 12)
tcga.cna = t(tcga.cna[,grep("TCGA", colnames(tcga.cna))])

### Integration
# Sort the matrices by samples/rownames
samples.with.mut.and.cna = sort(intersect(rownames(tcga.cna), rownames(tcga.mut.matrix)))
tcga.mut.matrix = tcga.mut.matrix[pmatch(samples.with.mut.and.cna, rownames(tcga.mut.matrix)),]
tcga.cna = tcga.cna[pmatch(samples.with.mut.and.cna, rownames(tcga.cna)),]

# Sort the matrices by genes/colnames
genes.with.mut.and.cna = sort(intersect(colnames(tcga.cna), colnames(tcga.mut.matrix)))
tcga.mut.matrix = tcga.mut.matrix[,pmatch(genes.with.mut.and.cna, colnames(tcga.mut.matrix))]
tcga.cna = tcga.cna[,pmatch(genes.with.mut.and.cna, colnames(tcga.cna))]

identical(dimnames(tcga.mut.matrix), dimnames(tcga.cna))



### Clinical data
# Load clinical data and remove additional header rows
tcga.clinical.data = read.table('~/Documents/Projects/TCGA_Data/stddata__2016_01_28/HNSC/20160128/gdac.broadinstitute.org_HNSC.Merge_Clinical.Level_1.2016012800.0.0/HNSC.clin.merged.txt', header = F, sep = "\t", stringsAsFactors = F, quote = "", fill = T)
tcga.clinical.data = t(tcga.clinical.data)
colnames(tcga.clinical.data) = tcga.clinical.data[1,]
tcga.clinical.data = tcga.clinical.data[-1,]
tcga.clinical.data = data.frame(tcga.clinical.data)

# Classify subsite/anatomic location
oral.cavity.locations = tolower(c("Alveolar Ridge", "Buccal Mucosa", "Floor of mouth", "Hard Palate", "Lip", "Oral Cavity", "Oral Tongue"))
oropharynx.locations = tolower(c("Base of tongue", "Oropharynx", "Tonsil"))
tcga.clinical.data$Location = as.character(tcga.clinical.data$patient.anatomic_neoplasm_subdivision)
tcga.clinical.data$Location[which(tcga.clinical.data$Location %in% oral.cavity.locations)] = "Oral Cavity"
tcga.clinical.data$Location[which(tcga.clinical.data$Location %in% oropharynx.locations)] = "Oropharynx"
tcga.clinical.data$Location[which(tcga.clinical.data$Location == "larynx")] = "Larynx"
tcga.clinical.data$Location[which(tcga.clinical.data$Location == "hypopharynx")] = "Hypopharynx"

# Select HPV neg samples with mut and CNA data
tcga.clinical.data$bcr_patient_barcode = toupper(as.character(tcga.clinical.data$patient.bcr_patient_barcode))
tcga.clinical.data$hpv_status = capitalize(as.character(tcga.clinical.data$patient.hpv_test_results.hpv_test_result.hpv_status))
tcga.clinical.data.hpv.neg = tcga.clinical.data[which(tcga.clinical.data$hpv_status == "Negative" & tcga.clinical.data$bcr_patient_barcode %in% rownames(tcga.mut.matrix)),]



save(tcga.mut.matrix, tcga.cna, tcga.clinical.data, tcga.clinical.data.hpv.neg,
     file = "~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Data/TCGA_Mut_and_CNA_small.Rdata")