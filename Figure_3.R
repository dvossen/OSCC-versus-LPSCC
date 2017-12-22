rm( list = ls( all = TRUE ) )

library(ggplot2)
library(reshape2)

load("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Data/TCGA_Mut_and_CNA_small.Rdata")
tcga.oscc.barcodes = tcga.clinical.data.hpv.neg$bcr_patient_barcode[which(tcga.clinical.data.hpv.neg$Location == "Oral Cavity")]
tcga.lpscc.barcodes = tcga.clinical.data.hpv.neg$bcr_patient_barcode[which(tcga.clinical.data.hpv.neg$Location != "Oral Cavity")]

tcga.gistic.lesions = read.table("~/Downloads/gdac.broadinstitute.org_HNSC-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/all_lesions.conf_99.txt", header = T, sep = "\t", stringsAsFactors = F, quote = "")
tcga.gistic.lesions.2 = tcga.gistic.lesions[-grep("Actual", tcga.gistic.lesions$Amplitude.Threshold), -c(1:9)]
colnames(tcga.gistic.lesions.2) = substring(gsub(pattern = ".", replacement = "-", x = colnames(tcga.gistic.lesions.2), fixed = T), 0, 12)
tcga.gistic.lesion.sum = apply(tcga.gistic.lesions.2[, which(colnames(tcga.gistic.lesions.2) %in% c(tcga.oscc.barcodes, tcga.lpscc.barcodes))], 2, function(x){length(which(x != 0))})

tcga.tissue = tcga.clinical.data.hpv.neg$Location[pmatch(names(tcga.gistic.lesion.sum), tcga.clinical.data.hpv.neg$bcr_patient_barcode)]
tcga.tissue = ifelse(tcga.tissue == "Oral Cavity", "OSCC", "L/P-SCC")

tcga.casp8 = tcga.mut.matrix[pmatch(names(tcga.gistic.lesion.sum), rownames(tcga.mut.matrix)), "CASP8"]
tcga.casp8 = ifelse(tcga.casp8 == 1, "mut", "wt")

tcga.tissue.casp8 = paste(tcga.tissue, tcga.casp8)

idx.casp8.mut = which(tcga.casp8 == "mut")
idx.casp8.wt = which(tcga.casp8 == "wt")

wilcox.test(tcga.gistic.lesion.sum ~ tcga.tissue)$p.val
wilcox.test(tcga.gistic.lesion.sum[-idx.casp8.mut] ~ tcga.tissue[-idx.casp8.mut])$p.val
wilcox.test(tcga.gistic.lesion.sum[idx.casp8.mut] ~ tcga.tissue[idx.casp8.mut])$p.val
wilcox.test(tcga.gistic.lesion.sum[-idx.casp8.mut] ~ tcga.tissue[-idx.casp8.mut])$p.val
wilcox.test(tcga.gistic.lesion.sum[which(tcga.tissue.casp8 == "L/P-SCC mut")],
            tcga.gistic.lesion.sum[which(tcga.tissue.casp8 == "OSCC wt")])$p.val
wilcox.test(tcga.gistic.lesion.sum[which(tcga.tissue.casp8 == "L/P-SCC mut")],
            tcga.gistic.lesion.sum[which(tcga.tissue.casp8 == "L/P-SCC wt")])$p.val



data.all = data.frame(count = tcga.gistic.lesion.sum, subsite = tcga.tissue, type = rep("CASP8 mut and wt", length(tcga.tissue)))
data.all = rbind(data.all, data.all)
data.all$Extra = rep(letters[1:2], each = nrow(data.all) / 2)
data.all = melt(data.all, id.vars = 2:4, measure.vars = 1)
data.all$subsite <- with(data.all, relevel(subsite, "OSCC"))

cols = gray.colors(2)
names(cols) = c("OSCC", "L/P-SCC")

gg <- ggplot(data.all, aes(x=subsite, y=value))
gg <- gg + geom_boxplot(outlier.shape=NA, aes(fill=subsite))
gg <- gg + geom_dotplot(binaxis='y', stackdir='center',  binwidth = 1/2, dotsize = 0.8)
gg <- gg + facet_wrap(~Extra)
gg <- gg + labs(x="")
gg <- gg + theme_bw()
gg <- gg + theme(strip.background=element_rect(fill="white"))
gg <- gg + theme(strip.text=element_text(color="black", face="bold"))
gg <- gg + theme(strip.text.x = element_text(size = 12))
gg <- gg + scale_fill_manual(values = cols) 
gg <- gg + theme(legend.position="none")
gg <- gg + theme(axis.text=element_text(size=12, face="bold"))
gg <- gg + theme(axis.title=element_text(size=24))
gg <- gg + ylab("Regions with CNA")
gg <- gg + scale_y_continuous(limits = c(0, 60))

pdf("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Figures/Figure_3A.pdf", family = "ArialMT", useDingbats=FALSE)
gg
dev.off()



### Figure 3B
data.stratified = data.frame(count = tcga.gistic.lesion.sum, subsite = paste(tcga.tissue, tcga.casp8, sep = "-"), type = rep("CASP8 stratified", length(tcga.tissue)))
data.stratified = rbind(data.stratified, data.stratified)
data.stratified$Extra = rep(letters[1:2], each = nrow(data.stratified) / 2)
data.stratified = melt(data.stratified, id.vars = 2:4, measure.vars = 1)
#data.stratified$subsite = factor(data.stratified$subsite, levels = c("L/P-SCC-mut", "OSCC-mut", "OSCC-wt", "L/P-SCC-wt"))
data.stratified$subsite = factor(data.stratified$subsite, levels = c("OSCC-mut", "L/P-SCC-mut", "OSCC-wt", "L/P-SCC-wt"))

cols = rep(gray.colors(2), each = 2)
names(cols) = c("OSCC-mut", "OSCC-wt", "L/P-SCC-mut", "L/P-SCC-wt")

gg <- ggplot(data.stratified, aes(x=subsite, y=value))
gg <- gg + geom_boxplot(outlier.shape=NA, aes(fill=subsite))
gg <- gg + geom_dotplot(binaxis='y', stackdir='center',  binwidth = 1/2, dotsize = 0.8)
gg <- gg + facet_wrap(~Extra)
gg <- gg + labs(x="")
gg <- gg + theme_bw()
gg <- gg + theme(strip.background=element_rect(fill="white"))
gg <- gg + theme(strip.text=element_text(color="black", face="bold"))
gg <- gg + theme(strip.text.x = element_text(size = 12))
gg <- gg + scale_fill_manual(values = cols) 
gg <- gg + theme(legend.position="none")
gg <- gg + theme(axis.text=element_text(size=12, face="bold"))
gg <- gg + theme(axis.title=element_text(size=24))
gg <- gg + ylab("Regions with CNA")
gg <- gg + scale_y_continuous(limits = c(0, 60))

pdf("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Figures/Figure_3B.pdf", family = "ArialMT", useDingbats=FALSE)
gg
dev.off()


### Figure 3C
tcga.HRAS = tcga.mut.matrix[pmatch(names(tcga.gistic.lesion.sum), rownames(tcga.mut.matrix)), "HRAS"]
tcga.HRAS = ifelse(tcga.HRAS == 1, "mut", "wt")

tcga.tissue.HRAS = paste(tcga.tissue, tcga.HRAS)

idx.HRAS.mut = which(tcga.HRAS == "mut")
idx.HRAS.wt = which(tcga.HRAS == "wt")

wilcox.test(tcga.gistic.lesion.sum ~ tcga.tissue)$p.val
wilcox.test(tcga.gistic.lesion.sum[-idx.HRAS.mut] ~ tcga.tissue[-idx.HRAS.mut])$p.val
wilcox.test(tcga.gistic.lesion.sum[idx.HRAS.mut] ~ tcga.tissue[idx.HRAS.mut])$p.val
wilcox.test(tcga.gistic.lesion.sum[which(tcga.tissue.HRAS == "L/P-SCC mut")],
            tcga.gistic.lesion.sum[which(tcga.tissue.HRAS == "OSCC wt")])$p.val
wilcox.test(tcga.gistic.lesion.sum[which(tcga.tissue.HRAS == "L/P-SCC mut")],
            tcga.gistic.lesion.sum[which(tcga.tissue.HRAS == "L/P-SCC wt")])$p.val
wilcox.test(tcga.gistic.lesion.sum[which(tcga.tissue.HRAS == "OSCC mut")],
            tcga.gistic.lesion.sum[which(tcga.tissue.HRAS == "OSCC wt")])$p.val


data.all = data.frame(count = tcga.gistic.lesion.sum, subsite = tcga.tissue, type = rep("HRAS mut and wt", length(tcga.tissue)))
data.all = rbind(data.all, data.all)
data.all$Extra = rep(letters[1:2], each = nrow(data.all) / 2)
data.all = melt(data.all, id.vars = 2:4, measure.vars = 1)
data.all$subsite <- with(data.all, relevel(subsite, "OSCC"))

data.stratified = data.frame(count = tcga.gistic.lesion.sum, subsite = paste(tcga.tissue, tcga.HRAS, sep = "-"), type = rep("HRAS stratified", length(tcga.tissue)))
data.stratified = rbind(data.stratified, data.stratified)
data.stratified$Extra = rep(letters[1:2], each = nrow(data.stratified) / 2)
data.stratified = melt(data.stratified, id.vars = 2:4, measure.vars = 1)
#data.stratified$subsite = factor(data.stratified$subsite, levels = c("L/P-SCC-mut", "OSCC-mut", "OSCC-wt", "L/P-SCC-wt"))
data.stratified$subsite = factor(data.stratified$subsite, levels = c("OSCC-mut", "L/P-SCC-mut", "OSCC-wt", "L/P-SCC-wt"))

cols = rep(gray.colors(2), each = 2)
names(cols) = c("OSCC-mut", "OSCC-wt", "L/P-SCC-mut", "L/P-SCC-wt")

gg <- ggplot(data.stratified, aes(x=subsite, y=value))
#gg <- gg + geom_boxplot(aes(fill=subsite))
gg <- gg + geom_boxplot(outlier.shape=NA, aes(fill=subsite))
gg <- gg + geom_dotplot(binaxis='y', stackdir='center',  binwidth = 1/2, dotsize = 0.8)
#gg <- gg + geom_jitter(position=position_jitter(width=.2, height=0), size = .5)
gg <- gg + facet_wrap(~Extra)
gg <- gg + labs(x="")
gg <- gg + theme_bw()
gg <- gg + theme(strip.background=element_rect(fill="white"))
gg <- gg + theme(strip.text=element_text(color="black", face="bold"))
gg <- gg + theme(strip.text.x = element_text(size = 12))
gg <- gg + scale_fill_manual(values = cols) 
gg <- gg + theme(legend.position="none")
gg <- gg + theme(axis.text=element_text(size=12, face="bold"))
gg <- gg + theme(axis.title=element_text(size=24))
gg <- gg + ylab("Regions with CNA")
gg <- gg + scale_y_continuous(limits = c(0, 60))

pdf("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Figures/Figure_3C.pdf", family = "ArialMT", useDingbats=FALSE)
gg
dev.off()