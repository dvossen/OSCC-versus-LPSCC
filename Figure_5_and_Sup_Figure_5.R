rm( list = ls( all = TRUE ) )

load("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Data/Figure_1A_data_TCGA.Rdata")

### Our canonical FA/HR genes
candidate.genes = read.csv("~/Documents/Projects/Old/Caroline/Documentation/List candidate genes.csv", stringsAsFactors = F)
candidate.genes.vector = unique(c(candidate.genes$geneNameMart, candidate.genes$Alternative.name))
our.canonical.fa.hr = candidate.genes.vector[-which(candidate.genes.vector == "")]
rm(candidate.genes)
rm(candidate.genes.vector)
# TCGA
tcga.oscc.canonical.fa.hr = ifelse(apply(tcga.oscc.mut[,which(colnames(tcga.oscc.mut) %in% our.canonical.fa.hr)], 1, sum) != 0, "Mut", "Wt")
tcga.lpscc.canonical.fa.hr = ifelse(apply(tcga.lpscc.mut[,which(colnames(tcga.lpscc.mut) %in% our.canonical.fa.hr)], 1, sum) != 0, "Mut", "Wt")
tcga.canonical.fa.hr = c(tcga.oscc.canonical.fa.hr, tcga.lpscc.canonical.fa.hr)
#rm(tcga.oscc.canonical.fa.hr)
#rm(tcga.lpscc.canonical.fa.hr)

load("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Data/TCGA_Mut_and_CNA_small.Rdata")
tcga.oscc.samples = tcga.clinical.data.hpv.neg$bcr_patient_barcode[which(tcga.clinical.data.hpv.neg$Location == "Oral Cavity")]
tcga.lpscc.samples = tcga.clinical.data.hpv.neg$bcr_patient_barcode[which(tcga.clinical.data.hpv.neg$Location != "Oral Cavity")]

# https://static-content.springer.com/esm/art%3A10.1186%2Fs40364-015-0033-4/MediaObjects/40364_2015_33_MOESM8_ESM.txt
table.scars = read.table(file = "~/Documents/Projects/Article_Comparative_genomics_OSCC/Figure_RACSs/Genomic_scar/40364_2015_33_MOESM8_ESM.txt", header = T, sep = "\t", stringsAsFactors = F)
table.scars.hnscc = table.scars[which(table.scars$Tumor %in% tcga.clinical.data.hpv.neg$bcr_patient_barcode),]
table.scars.hnscc$Subsite = "L/P-SCC"
table.scars.hnscc$Subsite[which(table.scars.hnscc$Tumor %in% tcga.oscc.samples)] = "OSCC"

table.scars.hnscc$FA_HR = tcga.canonical.fa.hr[pmatch(table.scars.hnscc$Tumor, names(tcga.canonical.fa.hr))]
table.scars.hnscc$FA_HR = relevel(as.factor(table.scars.hnscc$FA_HR), "Wt")
table(table.scars.hnscc$FA_HR)



library(ggplot2)
library(reshape2)
data.all = cbind(table.scars.hnscc, "Type" = rep("All samples", nrow(table.scars.hnscc)))
data.oscc = cbind(subset(table.scars.hnscc, Subsite == "OSCC"))
data.oscc$Type = data.oscc$Subsite
data.lpscc = cbind(subset(table.scars.hnscc, Subsite == "L/P-SCC"))
data.lpscc$Type = data.lpscc$Subsite
data = rbind(data.all, data.oscc, data.lpscc)
data$Subsite = as.factor(data$Subsite)
data$Subsite <- with(data, relevel(Subsite, "OSCC"))

wilcox.test(table.scars.hnscc$NtAI ~ table.scars.hnscc$FA_HR)
wilcox.test(table.scars.hnscc$LST ~ table.scars.hnscc$FA_HR)
wilcox.test(table.scars.hnscc$HRD.LOH ~ table.scars.hnscc$FA_HR)



pdf("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Figures/Figure_5.pdf", family = "ArialMT", useDingbats = FALSE)
data.x = melt(data, id.vars = c(1:2, 5:11), measure.vars = 2:4)
gg <- ggplot(data.x, aes(x=FA_HR, y=value, fill = "1"))
#gg <- gg + geom_boxplot(aes(fill=FA_HR))
gg <- gg + geom_boxplot(outlier.shape=NA)
gg <- gg + scale_fill_manual(values = '#A4A4A4')
gg <- gg + geom_dotplot(binaxis='y', stackdir='center',  binwidth = 1/2, dotsize = 0.8)
#gg <- gg + geom_jitter(position=position_jitter(width=.2, height=0), size = .5)
gg <- gg + facet_wrap(~variable*Type)
gg <- gg + labs(x="Canonical FA/HR gene set status")
gg <- gg + labs(y="Genomic scar signature score")
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

wilcox.test(data.oscc$NtAI ~ data.oscc$FA_HR)
wilcox.test(data.oscc$LST ~ data.oscc$FA_HR)
wilcox.test(data.oscc$HRD.LOH ~ data.oscc$FA_HR)

wilcox.test(data.lpscc$NtAI ~ data.lpscc$FA_HR)
wilcox.test(data.lpscc$LST ~ data.lpscc$FA_HR)
wilcox.test(data.lpscc$HRD.LOH ~ data.lpscc$FA_HR)



table(tcga.oscc.canonical.fa.hr)
table(tcga.lpscc.canonical.fa.hr)
tcga.fa.hr.table = matrix(c(69, 207, 47, 87), 2, 2, dimnames = list(c("Mut", "Wt"), c("OSCC", "L/P-SCC")))
tcga.fa.hr.prop.table = prop.table(tcga.fa.hr.table, 2)
pdf("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Figures/Sup_Figure_5.pdf", family = "ArialMT", useDingbats = FALSE)
par(mar = c(3, 6, 3, 1), mgp = c(3.5, 1, 0))
barplot(tcga.fa.hr.prop.table, beside = T, ylim = c(0, 1), ylab = "Proportion", main = "HR/FA Pathway",
        cex.lab = 3, cex.axis = 1.5, lwd = 3, las = 1, cex.main = 3, cex.names = 1.5)
legend("topright", c("Somatic Point Mutation", "Wildtype"), fill = grey.colors(2), cex = 1.5, bty = "n")
legend("topleft", paste("P-value =", round(fisher.test(tcga.fa.hr.table)$p, 3)), cex = 1.5, bty = "n")
dev.off()