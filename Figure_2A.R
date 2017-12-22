rm( list = ls( all = TRUE ) )

library(limma)
source("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Code/Begin_Code_functions.R")

# Load GenesSPM
load("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Data/Iorio_genes.Rdata")

### Somatic mutations
# NKI OSCC
load("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Data/NKI_Mut_and_CNA.Rdata")
nki.oscc.mut = nki.oscc.mut[,hnscc.cg]
nki.oscc.frac.mut = apply(nki.oscc.mut, 2, function(x){sum(x) / length(x) * 100})

# NKI L/P-SCC
nki.lpscc.mut = nki.lpscc.mut[,hnscc.cg]
nki.lpscc.frac.mut = apply(nki.lpscc.mut, 2, function(x){sum(x) / length(x) * 100})

# TCGA OSCC
load("~/Documents/Projects/Article_Comparative_genomics_OSCC/Figures_CGs/TCGA_Mut_and_CNA_small.Rdata")
tcga.oscc.mut = tcga.mut.matrix[which(rownames(tcga.mut.matrix) %in% tcga.clinical.data.hpv.neg$bcr_patient_barcode[which(tcga.clinical.data.hpv.neg$Location == "Oral Cavity")]),
                                which(colnames(tcga.mut.matrix) %in% hnscc.cg)]
tcga.oscc.frac.mut = apply(tcga.oscc.mut, 2, function(x){sum(x) / length(x) * 100})

tcga.oscc.mut.full = tcga.mut.matrix[which(rownames(tcga.mut.matrix) %in% tcga.clinical.data.hpv.neg$bcr_patient_barcode[which(tcga.clinical.data.hpv.neg$Location == "Oral Cavity")]),
                                which(colnames(tcga.mut.matrix) %in% iorio.hnscc.cancer.genes)]
tcga.oscc.frac.mut.full = apply(tcga.oscc.mut.full, 2, function(x){sum(x) / length(x) * 100})

# TCGA L/P-SCC
tcga.lpscc.mut = tcga.mut.matrix[which(rownames(tcga.mut.matrix) %in% tcga.clinical.data.hpv.neg$bcr_patient_barcode[which(tcga.clinical.data.hpv.neg$Location != "Oral Cavity")]),
                                which(colnames(tcga.mut.matrix) %in% hnscc.cg)]
tcga.lpscc.frac.mut = apply(tcga.lpscc.mut, 2, function(x){sum(x) / length(x) * 100})

tcga.lpscc.mut.full = tcga.mut.matrix[which(rownames(tcga.mut.matrix) %in% tcga.clinical.data.hpv.neg$bcr_patient_barcode[which(tcga.clinical.data.hpv.neg$Location != "Oral Cavity")]),
                                 which(colnames(tcga.mut.matrix) %in% iorio.hnscc.cancer.genes)]
tcga.lpscc.frac.mut.full = apply(tcga.lpscc.mut.full, 2, function(x){sum(x) / length(x) * 100})

identical(colnames(tcga.oscc.mut), colnames(nki.oscc.mut))
identical(colnames(tcga.lpscc.mut), colnames(nki.lpscc.mut))



### Plot
# Subsite comparison between datasets (NKI vs TCGA) - Log-scale
pdf("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Figures/Sup_Figure_2A.pdf", family = "ArialMT", useDingbats=FALSE)

nki.oscc.frac.mut.log = nki.oscc.frac.mut
nki.oscc.frac.mut.log[which(nki.oscc.frac.mut.log == 0)] = 0.1
nki.oscc.frac.mut.log = log10(nki.oscc.frac.mut.log)

nki.lpscc.frac.mut.log = nki.lpscc.frac.mut
nki.lpscc.frac.mut.log[which(nki.lpscc.frac.mut.log == 0)] = 0.1
nki.lpscc.frac.mut.log = log10(nki.lpscc.frac.mut.log)

tcga.oscc.frac.mut.log = log10(tcga.oscc.frac.mut)

tcga.lpscc.frac.mut.log = tcga.lpscc.frac.mut
tcga.lpscc.frac.mut.log[which(tcga.lpscc.frac.mut.log == 0)] = 0.1
tcga.lpscc.frac.mut.log = log10(tcga.lpscc.frac.mut.log)

par(mar = c(5, 6, 3, 2), mgp = c(3.5, 1, 0))
plot(x = tcga.oscc.frac.mut.log, y = nki.oscc.frac.mut.log,
     xlim = c(-1.1, 2), ylim = c(-1.1, 2),
     xlab = "% TCGA Mutated", ylab = "% NKI Mutated", main = expression("Genes"[SPM]), pch = 1,
     xaxs = "i", yaxs = "i", bty = "n",
     cex.main = 2, cex.lab = 3, cex.axis = 1.5, lwd = 1, las = 1, axes = F)
axis(side = 1, at = log10(c(0.1, 1, 10, 100)), 
     labels = c("0", "1", "10", "100"),
     cex.main = .6, cex.lab = 3, cex.axis = 1.5)
axis(side = 2, at = log10(c(0.1, 1, 10, 100)), 
     labels = c("0", "1", "10", "100"),
     cex.main = .6, cex.lab = 3, cex.axis = 1.5)
abline(a = 0, b = 1, lty = 2, lwd = 1)
points(x = tcga.lpscc.frac.mut.log, y = nki.lpscc.frac.mut.log, pch = 2, lwd = 1)

legend("topleft", 
       c(paste("OSCC (cor = ", round(cor(tcga.oscc.frac.mut, nki.oscc.frac.mut, method = "spearman"), 2), ")", sep = ""),
         paste("L/P-SCC (cor = ", round(cor(tcga.lpscc.frac.mut, nki.lpscc.frac.mut, method = "spearman"), 2), ")", sep = "")), 
       pch = c(1, 2), bty = "n", cex = 1, pt.lwd = 1)
dev.off()




# Barplot of significant genes
pdf("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Figures/Figure_2A.pdf", family = "ArialMT", useDingbats=FALSE)

p.nki = matrix.fisher(nki.oscc.mut, nki.lpscc.mut)
colnames(nki.oscc.mut)[which(p.adjust(p.nki, method = "fdr") < .1)]
idx.nki.significant = which(p.adjust(p.nki, method = "fdr") < .1)

p.tcga = matrix.fisher(tcga.oscc.mut.full, tcga.lpscc.mut.full)
colnames(tcga.oscc.mut.full)[which(p.adjust(p.tcga, method = "fdr") < .1)]
idx.tcga.significant = which(p.adjust(p.tcga, method = "fdr") < .1)
p.adjust(p.tcga, method = "fdr")[idx.tcga.significant]

names.significant = names(tcga.lpscc.frac.mut.full)[idx.tcga.significant]

m.tcga = matrix(c(tcga.oscc.frac.mut.full[names.significant], tcga.lpscc.frac.mut.full[names.significant]), byrow = T, nrow = 2)
m.nki = matrix(c(nki.oscc.frac.mut[names.significant], nki.lpscc.frac.mut[names.significant]), byrow = T, nrow = 2)

par(mar = c(2, 5.5, 3, 0.5), cex = 1, mfrow = c(2, 1))

barplot.tcga = barplot(m.tcga, ylim = c(0, 30), names = names.significant, 
                       ylab = "% TCGA Mutated", width = 1, main = expression("Genes"[SPM]), cex.main = 3,
                       horiz = F, beside = T, cex.names = 2, cex.axis = 1.5, cex.lab = 2, las = 1) #xaxt = "n",
text(barplot.tcga, c(m.tcga) + 1.8, labels = paste(round(c(m.tcga), 0), "%", sep = ""), cex = 1.3)



legend("topleft", legend = c("OSCC", "L/P-SCC"), fill = gray.colors(2), bty = "n", cex = 1.3)



barplot.nki = barplot(-m.nki, ylim = c(-30, 0), 
        ylab = "% NKI Mutated", width = 1, #names = names.significant, #main = "NKI",
        horiz = F, yaxt = "n", beside = T, cex.axis = 1.5, cex.lab = 2, las = 1)
text(barplot.nki, -c(m.nki) - 1.8, labels = paste(round(c(m.nki), 0), "%", sep = ""), cex = 1.3)

axis(side = 2, at = seq(from = -30, to = 0, by = 5), 
     labels = rev(seq(from = 0, to = 30, by = 5)), cex.axis = 1.5, cex.lab = 2, las = 1)

dev.off()