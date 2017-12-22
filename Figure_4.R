rm( list = ls( all = TRUE ) )

library(ggplot2)

load("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Data/TCGA_Mut_and_CNA_small.Rdata")
tcga.oscc.samples = tcga.clinical.data.hpv.neg$bcr_patient_barcode[which(tcga.clinical.data.hpv.neg$Location == "Oral Cavity")]
tcga.lpscc.samples = tcga.clinical.data.hpv.neg$bcr_patient_barcode[which(tcga.clinical.data.hpv.neg$Location != "Oral Cavity")]

# https://static-content.springer.com/esm/art%3A10.1186%2Fs40364-015-0033-4/MediaObjects/40364_2015_33_MOESM8_ESM.txt
table.scars = read.table(file = "~/Documents/Projects/Article_Comparative_genomics_OSCC/Figure_RACSs/Genomic_scar/40364_2015_33_MOESM8_ESM.txt", header = T, sep = "\t", stringsAsFactors = F)
table.scars.hnscc = table.scars[which(table.scars$Tumor %in% tcga.clinical.data.hpv.neg$bcr_patient_barcode),]
table.scars.hnscc$Subsite = "L/P-SCC"
table.scars.hnscc$Subsite[which(table.scars.hnscc$Tumor %in% tcga.oscc.samples)] = "OSCC"

table.scars.hnscc$CASP8 = tcga.mut.matrix[pmatch(table.scars.hnscc$Tumor, rownames(tcga.mut.matrix)), "CASP8"]
table.scars.hnscc$CASP8 = ifelse(table.scars.hnscc$CASP8 == 1, "mut", "wt")
idx.casp8.mut = which(table.scars.hnscc$CASP8 == "mut")

table.scars.hnscc$HRD.classification = ifelse(table.scars.hnscc$HRD.LOH>= median(table.scars.hnscc$HRD.LOH), "High HRD", "Low HRD")



table(table.scars.hnscc$Subsite)

wilcox.test(table.scars.hnscc$NtAI ~ table.scars.hnscc$Subsite)
wilcox.test(table.scars.hnscc$LST ~ table.scars.hnscc$Subsite)
wilcox.test(table.scars.hnscc$HRD.LOH ~ table.scars.hnscc$Subsite)

wilcox.test(table.scars.hnscc$NtAI[-idx.casp8.mut] ~ table.scars.hnscc$Subsite[-idx.casp8.mut])
wilcox.test(table.scars.hnscc$LST[-idx.casp8.mut] ~ table.scars.hnscc$Subsite[-idx.casp8.mut])
wilcox.test(table.scars.hnscc$HRD.LOH[-idx.casp8.mut] ~ table.scars.hnscc$Subsite[-idx.casp8.mut])

idx.oscc = which(table.scars.hnscc$Subsite == "OSCC")
wilcox.test(table.scars.hnscc$NtAI[idx.oscc] ~ table.scars.hnscc$CASP8[idx.oscc])
wilcox.test(table.scars.hnscc$LST[idx.oscc] ~ table.scars.hnscc$CASP8[idx.oscc])
wilcox.test(table.scars.hnscc$NtAI[idx.oscc] ~ table.scars.hnscc$CASP8[idx.oscc])


### Compare proportion of samples with a high HRD score
# All samples
oscc.scars = subset(table.scars.hnscc, Subsite == 'OSCC')
lpscc.scars = subset(table.scars.hnscc, Subsite == 'L/P-SCC')
round(prop.table(table(oscc.scars$HRD.classification)) * 100, 0)
round(prop.table(table(lpscc.scars$HRD.classification)) * 100, 0)
fisher.test(table(table.scars.hnscc$Subsite, table.scars.hnscc$HRD.classification))

# CASP8 wildtype samples
hnscc.scars.casp8.wt = subset(table.scars.hnscc, CASP8 == 'wt')
oscc.scars.casp8.wt = subset(oscc.scars, CASP8 == 'wt')
lpscc.scars.casp8.wt = subset(lpscc.scars, CASP8 == 'wt')
round(prop.table(table(oscc.scars.casp8.wt$HRD.classification)) * 100, 0)
round(prop.table(table(lpscc.scars.casp8.wt$HRD.classification)) * 100, 0)
fisher.test(table(hnscc.scars.casp8.wt$Subsite, hnscc.scars.casp8.wt$HRD.classification))



### Figures
data.all = cbind(table.scars.hnscc, "Type" = rep("CASP8 mut and wt", nrow(table.scars.hnscc)))
data.wt = cbind(table.scars.hnscc[-idx.casp8.mut,], "Type" = rep("CASP8 wt", nrow(table.scars.hnscc[-idx.casp8.mut,])))
data.mut = cbind(table.scars.hnscc[idx.casp8.mut,], "Type" = rep("CASP8 mut", nrow(table.scars.hnscc[idx.casp8.mut,])))
data = rbind(data.all, data.wt, data.mut)
data$Subsite = as.factor(data$Subsite)
data$Subsite <- with(data, relevel(Subsite, "OSCC"))

cols = gray.colors(2)
names(cols) = c("OSCC", "L/P-SCC")

pdf("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Figures/Figure_4A.pdf", family = "ArialMT", useDingbats = FALSE)
library(reshape2)

data.x = melt(data, id.vars = c(1:2, 5:12), measure.vars = 2:4)

gg <- ggplot(data.x, aes(x=Subsite, y=value))
#gg <- gg + geom_boxplot(aes(fill=Subsite))
gg <- gg + geom_boxplot(outlier.shape=NA, aes(fill=Subsite))
gg <- gg + geom_dotplot(binaxis='y', stackdir='center',  binwidth = 1/2, dotsize = 2)
#gg <- gg + geom_jitter(position=position_jitter(width=.2, height=0), size = .5)
gg <- gg + facet_wrap(~variable*Type)
gg <- gg + labs(x="")
gg <- gg + theme_bw()
gg <- gg + theme(strip.background=element_rect(fill="white"))
gg <- gg + theme(strip.text=element_text(color="black", face="bold"))
gg <- gg + theme(strip.text.x = element_text(size = 12))
gg <- gg + scale_fill_manual(values = cols) 
gg <- gg + theme(legend.position="none")
gg <- gg + theme(axis.text=element_text(size=12, face="bold"))
gg <- gg + theme(axis.title=element_text(size=24))
gg <- gg + scale_y_continuous(limits = c(0, 50))
plot(gg)

dev.off()



data.splitted = table.scars.hnscc
data.splitted$Subsite = paste(data.splitted$Subsite, data.splitted$CASP8, sep = "-")
data.splitted$Subsite = ifelse(data.splitted$Subsite == "OSCC-mut", "Om", ifelse(data.splitted$Subsite == "OSCC-wt", "Owt", "Lwt"))
data.splitted$Subsite = factor(data.splitted$Subsite, levels = c("Om", "Owt", "Lwt"))
data.splitted = rbind(data.splitted, data.splitted, data.splitted)
data.splitted$Extra = rep(letters[1:3], each = nrow(data.splitted) / 3)
data.y = melt(data.splitted, id.vars = c(1, 5:12), measure.vars = 2:4)

cols = gray.colors(2)[c(2, 1, 1)]
names(cols) = c("Lwt", "Om", "Owt")

gg <- ggplot(data.y, aes(x=Subsite, y=value))
#gg <- gg + geom_boxplot(aes(fill=Subsite))
gg <- gg + geom_boxplot(outlier.shape=NA, aes(fill=Subsite))
gg <- gg + geom_dotplot(binaxis='y', stackdir='center',  binwidth = 1/2, dotsize = 2)
#gg <- gg + geom_jitter(position=position_jitter(width=.2, height=0), size = .5)
gg <- gg + facet_wrap(~variable*Extra)
gg <- gg + labs(x="")
gg <- gg + theme_bw()
gg <- gg + theme(strip.background=element_rect(fill="white"))
gg <- gg + theme(strip.text=element_text(color="black", face="bold"))
gg <- gg + theme(strip.text.x = element_text(size = 12))
gg <- gg + scale_fill_manual(values = cols) 
gg <- gg + theme(legend.position="none")
gg <- gg + theme(axis.text=element_text(size=12, face="bold"))
gg <- gg + theme(axis.title=element_text(size=24))
gg <- gg + scale_y_continuous(limits = c(0, 50))

pdf("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Figures/Figure_4B.pdf", family = "ArialMT", useDingbats = FALSE)
plot(gg)
dev.off()






table.scars.hnscc$HRAS = tcga.mut.matrix[pmatch(table.scars.hnscc$Tumor, rownames(tcga.mut.matrix)), "HRAS"]
table.scars.hnscc$HRAS = ifelse(table.scars.hnscc$HRAS == 1, "mut", "wt")
idx.HRAS.mut = which(table.scars.hnscc$HRAS == "mut")

data.splitted = table.scars.hnscc
data.splitted$Subsite = paste(data.splitted$Subsite, data.splitted$HRAS, sep = "-")
data.splitted$Subsite = ifelse(data.splitted$Subsite == "OSCC-mut", "Om", ifelse(data.splitted$Subsite == "OSCC-wt", "Owt", "Lwt"))
data.splitted$Subsite = factor(data.splitted$Subsite, levels = c("Om", "Owt", "Lwt"))
data.splitted = rbind(data.splitted, data.splitted, data.splitted)
data.splitted$Extra = rep(letters[1:3], each = nrow(data.splitted) / 3)
data.y = melt(data.splitted, id.vars = c(1, 5:13), measure.vars = 2:4)

cols = gray.colors(2)[c(2, 1, 1)]
names(cols) = c("Lwt", "Om", "Owt")

gg <- ggplot(data.y, aes(x=Subsite, y=value))
#gg <- gg + geom_boxplot(aes(fill=Subsite))
gg <- gg + geom_boxplot(outlier.shape=NA, aes(fill=Subsite))
gg <- gg + geom_dotplot(binaxis='y', stackdir='center',  binwidth = 1/2, dotsize = 2)
#gg <- gg + geom_jitter(position=position_jitter(width=.2, height=0), size = .5)
gg <- gg + facet_wrap(~variable*Extra)
gg <- gg + labs(x="")
gg <- gg + theme_bw()
gg <- gg + theme(strip.background=element_rect(fill="white"))
gg <- gg + theme(strip.text=element_text(color="black", face="bold"))
gg <- gg + theme(strip.text.x = element_text(size = 12))
gg <- gg + scale_fill_manual(values = cols) 
gg <- gg + theme(legend.position="none")
gg <- gg + theme(axis.text=element_text(size=12, face="bold"))
gg <- gg + theme(axis.title=element_text(size=24))
gg <- gg + scale_y_continuous(limits = c(0, 50))

pdf("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Figures/Figure_4C.pdf", family = "ArialMT", useDingbats = FALSE)
plot(gg)
dev.off()