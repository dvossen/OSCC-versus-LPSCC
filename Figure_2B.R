rm( list = ls( all = TRUE ) )

library(ggplot2)
library(GenVisR)
library(plyr)
library(reshape2)

source("~/Documents/Projects/External_R_Code/GenVisR-master/R/TvTi_alignPlot.R")
source("~/Documents/Projects/External_R_Code/GenVisR-master/R/TvTi_annoTransTranv.R")
source("~/Documents/Projects/External_R_Code/GenVisR-master/R/TvTi_buildMain.R")
source("~/Documents/Projects/External_R_Code/GenVisR-master/R/TvTi_calcTransTranvFreq.R")
source("~/Documents/Projects/External_R_Code/GenVisR-master/R/TvTi_convMaf.R")
source("~/Documents/Projects/External_R_Code/GenVisR-master/R/TvTi_qual.R")
source("~/Documents/Projects/External_R_Code/GenVisR-master/R/TvTi_rmIndel.R")
source("~/Documents/Projects/External_R_Code/GenVisR-master/R/TvTi_rmMnuc.R")
source("~/Documents/Projects/External_R_Code/GenVisR-master/R/TvTi.R")
source("~/Documents/Projects/External_R_Code/GenVisR-master/R/multi_selectOut.R")
source("~/Documents/Projects/External_R_Code/GenVisR-master/R/multi_buildClin.R")
load("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Data/TCGA_Mut_and_CNA_small.Rdata")

### Mutations
# Load mutation file
tcga.mut = read.table("~/Documents/Projects/TCGA_GDC_legacy_archive/HNSC_pairs.aggregated.capture.tcga.uuid.automated.somatic.maf", header = T, sep = '\t', stringsAsFactors = F, quote = "", skip = 4)
# Select samples from primary tumors
tcga.mut = tcga.mut[which(substring(tcga.mut$Tumor_Sample_Barcode, 14, 15) == "01"),]
# Rename samples
tcga.mut$Tumor_Sample_Barcode = substring(tcga.mut$Tumor_Sample_Barcode, 0, 12)
# Select only HPV-negative samples
tcga.mut.hpv.neg = tcga.mut[which(tcga.mut$Tumor_Sample_Barcode %in% tcga.clinical.data.hpv.neg$bcr_patient_barcode),]
# Remove RNA mutations
tcga.mut.hpv.neg.non.RNA = tcga.mut.hpv.neg[-which(tcga.mut.hpv.neg$Variant_Classification %in% c("RNA")),]

df.mut.sum = data.frame(Tumor_Sample_Barcode = names(table(tcga.mut.hpv.neg.non.RNA$Tumor_Sample_Barcode)), 
                 Sum = as.numeric(unname(table(tcga.mut.hpv.neg.non.RNA$Tumor_Sample_Barcode))))
df.mut.sum$Subsite = tcga.clinical.data.hpv.neg$Location[pmatch(df.mut.sum$Tumor_Sample_Barcode, tcga.clinical.data.hpv.neg$bcr_patient_barcode)]
df.mut.sum$Subsite = ifelse(df.mut.sum$Subsite == "Oral Cavity", "OSCC", "L/P-SCC")
df.mut.sum = arrange(df = df.mut.sum, Sum)
df.mut.sum.oscc = subset(df.mut.sum, Subsite == "OSCC")
df.mut.sum.lpscc = subset(df.mut.sum, Subsite == "L/P-SCC")
df.mut.sum2 = rbind(df.mut.sum.oscc, df.mut.sum.lpscc)



tcga.mut.hpv.neg.non.RNA$Order1 = pmatch(tcga.mut.hpv.neg.non.RNA$Tumor_Sample_Barcode, df.mut.sum$Tumor_Sample_Barcode, duplicates.ok = T)
tcga.mut.hpv.neg.non.RNA$Order2 = pmatch(tcga.mut.hpv.neg.non.RNA$Tumor_Sample_Barcode, df.mut.sum2$Tumor_Sample_Barcode, duplicates.ok = T)
tcga.mut.hpv.neg.non.RNA1 = arrange(df = tcga.mut.hpv.neg.non.RNA, Order1)
tcga.mut.hpv.neg.non.RNA2 = arrange(df = tcga.mut.hpv.neg.non.RNA, Order2)



input.clin.data = data.frame(sample = df.mut.sum$Tumor_Sample_Barcode,
                             variable = rep("Subsite", nrow(df.mut.sum)),
                             value = df.mut.sum$Subsite)
idx.more.than.1000.mut = which(df.mut.sum$Sum >= 1e3)
samples.more.than.1000.mut = df.mut.sum$Tumor_Sample_Barcode[idx.more.than.1000.mut]
color.vector = grey.colors(2)
names(color.vector) =  c("OSCC", "L/P-SCC")
pdf("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Figures/Figure_2B.pdf")
TvTi(tcga.mut.hpv.neg.non.RNA1[-which(tcga.mut.hpv.neg.non.RNA1$Tumor_Sample_Barcode %in% samples.more.than.1000.mut),], 
     lab_txtAngle=75, fileType="MAF", type = "Frequency", 
     clinData = input.clin.data[-idx.more.than.1000.mut,], lab_Xaxis = F, clinVarCol = color.vector)
dev.off()



median(df.mut.sum.oscc$Sum)
median(df.mut.sum.lpscc$Sum)
df.mut.sum$Subsite = relevel(as.factor(df.mut.sum$Subsite), "OSCC")
pdf("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Figures/Figure_2B-inlet.pdf")
par(mgp = c(0, 1, 0))
par(mar = c(6, 8, 3, 2))
boxplot(Sum ~ Subsite, data = df.mut.sum, outline = F, col = grey.colors(2),
        axes=FALSE, cex.lab = 1, cex.axis = 1, lwd = 1, las = 1, frame = F, ylim = c(0, 800))
axis(2, at = c(0, 2e2, 4e2, 6e2, 8e2), labels = c(0, 2e2, 4e2, 6e2, 8e2),
     cex.lab = 1, cex.axis = 1, lwd = 1, las = 1)
dev.off()
wilcox.test(Sum ~ Subsite, data = df.mut.sum)