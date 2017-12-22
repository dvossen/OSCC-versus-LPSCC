rm( list = ls( all = TRUE ) )

### Mutations
# Load mutation file
tcga.mut = read.table("~/Documents/Projects/TCGA_GDC_legacy_archive/HNSC_pairs.aggregated.capture.tcga.uuid.automated.somatic.maf", header = T, sep = '\t', stringsAsFactors = F, quote = "", skip = 4)
# Attach "chr" to chromosome
tcga.mut$Chromosome_2 = paste("chr", tcga.mut$Chromosome, sep = "")
# Select samples from primary tumors
tcga.mut = tcga.mut[which(substring(tcga.mut$Tumor_Sample_Barcode, 14, 15) == "01"),]
# Rename samples
tcga.mut$Tumor_Sample_Barcode = substring(tcga.mut$Tumor_Sample_Barcode, 0, 12)
# Remove RNA mutations
tcga.mut.non.RNA = tcga.mut[-which(tcga.mut$Variant_Classification %in% c("RNA")),]
tcga.mut.non.RNA.non.chrMT = tcga.mut.non.RNA[-which(tcga.mut.non.RNA$Chromosome == "MT"),]

library(deconstructSigs)
sigs.input <- mut.to.sigs.input(mut.ref = tcga.mut.non.RNA.non.chrMT,
                                sample.id = "Tumor_Sample_Barcode",
                                chr = "Chromosome_2",
                                pos = "Start_position",
                                ref = "Reference_Allele",
                                alt = "Tumor_Seq_Allele2")
result.cosmic = NULL
for(i in 1:nrow(sigs.input)){
  temp.result.cosmic = whichSignatures(tumor.ref = sigs.input,
                                           signatures.ref = signatures.cosmic,
                                           sample.id = rownames(sigs.input)[i],
                                           contexts.needed = TRUE,
                                           tri.counts.method = 'default')
  result.cosmic = rbind(result.cosmic, temp.result.cosmic$weights)
}
samples.fewer.than.50.mut = names(which(rowSums(sigs.input) < 50))

### cosmic
cosmic.binary = result.cosmic
cosmic.binary = ifelse(cosmic.binary == 0, 0, 1)
cosmic.binary = cosmic.binary[-which(rownames(cosmic.binary) %in% samples.fewer.than.50.mut),]
cosmic.binary = as.data.frame(cosmic.binary)
cosmic.binary$SampleID = rownames(cosmic.binary)
load("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Data/TCGA_Mut_and_CNA_small.Rdata")
cosmic.binary$Subsite = tcga.clinical.data.hpv.neg$Location[pmatch(cosmic.binary$SampleID, tcga.clinical.data.hpv.neg$bcr_patient_barcode)]
cosmic.binary$Subsite = as.factor(ifelse(cosmic.binary$Subsite == "Oral Cavity", "OSCC", "L/P-SCC"))
cosmic.binary$Subsite = with(cosmic.binary, relevel(Subsite, "OSCC"))
age.diag = tcga.clinical.data.hpv.neg$patient.age_at_initial_pathologic_diagnosis
tcga.clinical.data.hpv.neg$Age = as.numeric(levels(age.diag))[age.diag]
cosmic.binary$Age = tcga.clinical.data.hpv.neg$Age[pmatch(cosmic.binary$SampleID, tcga.clinical.data.hpv.neg$bcr_patient_barcode)]
cosmic.binary.hpvneg = cosmic.binary[-which(is.na(cosmic.binary$Subsite)),]

# Fisher's exact test
cosmic.signatures.pval = rep(NA, 30)
names(cosmic.signatures.pval) = colnames(cosmic.binary.hpvneg)[1:30]
for(i in 1:30){
  temp.table = table(cosmic.binary.hpvneg$Subsite, cosmic.binary.hpvneg[,i])
  if(ncol(temp.table) != 1){
    cosmic.signatures.pval[i] = fisher.test(temp.table)$p 
  } else{
    next()
  }
}



# "We then focused on signatures that made a contribution to at least 25% of OSCC or L/P-SCC samples"
cosmic.binary.hpvneg.oscc = subset(cosmic.binary.hpvneg, Subsite == "OSCC")
cosmic.binary.hpvneg.lpscc = subset(cosmic.binary.hpvneg, Subsite == "L/P-SCC")
signatures.25pct = names(which(round(colSums(cosmic.binary.hpvneg.oscc[,1:30]) / nrow(cosmic.binary.hpvneg.oscc) * 100, 0) >= 25 | 
                                                      round(colSums(cosmic.binary.hpvneg.lpscc[,1:30]) / nrow(cosmic.binary.hpvneg.lpscc) * 100, 0) >= 25))
names.significant.sign.25pct = names(which(p.adjust(cosmic.signatures.pval[signatures.25pct]) < .1))


matrix.fraction.signature = rbind(colSums(cosmic.binary.hpvneg.oscc[,1:30]) / nrow(cosmic.binary.hpvneg.oscc),
                                  colSums(cosmic.binary.hpvneg.lpscc[,1:30]) / nrow(cosmic.binary.hpvneg.lpscc)) * 100
rownames(matrix.fraction.signature) = c("OSCC", "L/P-SCC")
### Select significant signatures
matrix.fraction.signature = matrix.fraction.signature[,names.significant.sign.25pct]

pdf("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Figures/Figure_2C.pdf")
par(mar = c(5, 5.5, 0.5, 0.5))
barplot.tcga = barplot(matrix.fraction.signature,
                       ylim = c(0, 120), names = sapply(colnames(matrix.fraction.signature), function(x){unlist(strsplit(x, "Signature.", fixed = T))[2]}),
                       xlab = "COSMIC Signature", ylab = "Present in TCGA (%)", width = 1, cex.main = 3,
                       horiz = F, beside = T, cex.names = 2, cex.axis = 1.5, cex.lab = 2, las = 1, yaxt = "n")
text(barplot.tcga, c(matrix.fraction.signature) + 3, labels = paste(round(c(matrix.fraction.signature), 0), "%", sep = ""), cex = 1)
legend("topright", legend = c("OSCC", "L/P-SCC"), fill = gray.colors(2), bty = "n", cex = 1.3)
axis(side = 2, at = seq(from = 0, to = 100, by = 20), 
     labels = seq(from = 0, to = 100, by = 20), cex.axis = 1.5, cex.lab = 2, las = 1)
dev.off()



### Association mutational and genomic scar signatures
load("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Data/TCGA_Mut_and_CNA_small.Rdata")
tcga.oscc.samples = tcga.clinical.data.hpv.neg$bcr_patient_barcode[which(tcga.clinical.data.hpv.neg$Location == "Oral Cavity")]
tcga.lpscc.samples = tcga.clinical.data.hpv.neg$bcr_patient_barcode[which(tcga.clinical.data.hpv.neg$Location != "Oral Cavity")]
# https://static-content.springer.com/esm/art%3A10.1186%2Fs40364-015-0033-4/MediaObjects/40364_2015_33_MOESM8_ESM.txt
table.scars = read.table(file = "~/Documents/Projects/Article_Comparative_genomics_OSCC/Figure_RACSs/Genomic_scar/40364_2015_33_MOESM8_ESM.txt", header = T, sep = "\t", stringsAsFactors = F)
table.scars.hnscc = table.scars[which(table.scars$Tumor %in% tcga.clinical.data.hpv.neg$bcr_patient_barcode),]
table.scars.hnscc$Subsite = "L/P-SCC"
table.scars.hnscc$Subsite[which(table.scars.hnscc$Tumor %in% tcga.oscc.samples)] = "OSCC"
table.scars.hnscc$Subsite <- with(table.scars.hnscc, relevel(Subsite, "OSCC"))
table.scars.hnscc$COSMIC.Signature.3 = cosmic.binary.hpvneg$Signature.3[pmatch(table.scars.hnscc$Tumor, rownames(cosmic.binary.hpvneg))]
table.scars.hnscc$COSMIC.Signature.3.char = ifelse(table.scars.hnscc$COSMIC.Signature.3 == 1, "Present", "Absent")
table.scars.hnscc.with.sign = table.scars.hnscc[-which(is.na(table.scars.hnscc$COSMIC.Signature.3)),]



data.x = melt(table.scars.hnscc.with.sign, id.vars = c(1:2, 5:11), measure.vars = 2:4)
pdf("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Figures/Sup_Figure_4.pdf", family = "ArialMT", useDingbats = FALSE)
gg <- ggplot(data.x, aes(x=COSMIC.Signature.3.char, y=value))
#gg <- gg + geom_boxplot(aes(fill=Subsite))
gg <- gg + geom_boxplot(outlier.shape=NA, aes(fill=COSMIC.Signature.3.char))
gg <- gg + geom_dotplot(binaxis='y', stackdir='center',  binwidth = 1/2, dotsize = .8)
#gg <- gg + geom_jitter(position=position_jitter(width=.2, height=0), size = .5)
gg <- gg + facet_wrap(~variable)
gg <- gg + labs(x="COSMIC Signature 3")
gg <- gg + labs(y="Scar signature score")
gg <- gg + theme_bw()
gg <- gg + theme(strip.background=element_rect(fill="white"))
gg <- gg + theme(strip.text=element_text(color="black", face="bold"))
gg <- gg + theme(strip.text.x = element_text(size = 12))
#gg <- gg + scale_fill_manual(values = cols) 
gg <- gg + theme(legend.position="none")
gg <- gg + theme(axis.text=element_text(size=12, face="bold"))
gg <- gg + theme(axis.title=element_text(size=24))
gg <- gg + scale_y_continuous(limits = c(0, 50))
plot(gg)
dev.off()

wilcox.test(table.scars.hnscc$NtAI ~ table.scars.hnscc$COSMIC.Signature.3)
wilcox.test(table.scars.hnscc$LST ~ table.scars.hnscc$COSMIC.Signature.3)
wilcox.test(table.scars.hnscc$HRD.LOH ~ table.scars.hnscc$COSMIC.Signature.3)