rm( list = ls( all = TRUE ) )

library(limma)
source("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Code/Begin_Code_functions.R")

# Load GenesCNA
load("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Data/Iorio_genes.Rdata")

### Copy number alterations
# NKI OSCC
load("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Data/NKI_Mut_and_CNA.Rdata")
nki.oscc.gain = ifelse(nki.oscc.cn[,hnscc.racg.gain.captured] >= 7, 1, 0)
nki.oscc.loss = ifelse(nki.oscc.cn[,hnscc.racg.loss.captured] <= 0.5, 1, 0)
nki.oscc.racg = cbind(nki.oscc.gain, nki.oscc.loss)

nki.oscc.frac = apply(nki.oscc.racg, 2, function(x){round(sum(x) / length(x) * 100, 2)})

# NKI L/P-SCC
nki.lpscc.gain = ifelse(nki.lpscc.cn[,hnscc.racg.gain.captured] >= 7, 1, 0)
nki.lpscc.loss = ifelse(nki.lpscc.cn[,hnscc.racg.loss.captured] <= 0.5, 1, 0)
nki.lpscc.racg = cbind(nki.lpscc.gain, nki.lpscc.loss)

nki.lpscc.frac = apply(nki.lpscc.racg, 2, function(x){round(sum(x) / length(x) * 100, 2)})

# TCGA OSCC
load("~/Documents/Projects/Article_Comparative_genomics_OSCC/Figures_CGs/TCGA_Mut_and_CNA_small.Rdata")
tcga.oscc.cn = tcga.cna[which(rownames(tcga.cna) %in% tcga.clinical.data.hpv.neg$bcr_patient_barcode[which(tcga.clinical.data.hpv.neg$Location == "Oral Cavity")]),]

tcga.oscc.gain = ifelse(tcga.oscc.cn[,hnscc.racg.gain.captured] == 2, 1, 0)
tcga.oscc.loss = ifelse(tcga.oscc.cn[,hnscc.racg.loss.captured] == -2, 1, 0)
tcga.oscc.racg = cbind(tcga.oscc.gain, tcga.oscc.loss)

tcga.oscc.frac = apply(tcga.oscc.racg, 2, function(x){round(sum(x) / length(x) * 100, 2)})


tcga.oscc.gain.full = ifelse(tcga.oscc.cn[,hnscc.racg.gain] == 2, 1, 0)
tcga.oscc.loss.full = ifelse(tcga.oscc.cn[,hnscc.racg.loss] == -2, 1, 0)
tcga.oscc.racg.full = cbind(tcga.oscc.gain.full, tcga.oscc.loss.full)

tcga.oscc.frac.full = apply(tcga.oscc.racg.full, 2, function(x){round(sum(x) / length(x) * 100, 2)})

# TCGA L/P-SCC
tcga.lpscc.cn = tcga.cna[which(rownames(tcga.cna) %in% tcga.clinical.data.hpv.neg$bcr_patient_barcode[which(tcga.clinical.data.hpv.neg$Location != "Oral Cavity")]),]

tcga.lpscc.frac.gain = ifelse(tcga.lpscc.cn[,hnscc.racg.gain.captured] == 2, 1, 0)
tcga.lpscc.frac.loss = ifelse(tcga.lpscc.cn[,hnscc.racg.loss.captured] == -2, 1, 0)
tcga.lpscc.racg = cbind(tcga.lpscc.frac.gain, tcga.lpscc.frac.loss)

tcga.lpscc.frac = apply(tcga.lpscc.racg, 2, function(x){round(sum(x) / length(x) * 100, 2)})


tcga.lpscc.gain.full = ifelse(tcga.lpscc.cn[,hnscc.racg.gain] == 2, 1, 0)
tcga.lpscc.loss.full = ifelse(tcga.lpscc.cn[,hnscc.racg.loss] == -2, 1, 0)
tcga.lpscc.racg.full = cbind(tcga.lpscc.gain.full, tcga.lpscc.loss.full)

tcga.lpscc.frac.full = apply(tcga.lpscc.racg.full, 2, function(x){round(sum(x) / length(x) * 100, 2)})

identical(names(tcga.oscc.frac), names(nki.oscc.frac))
identical(names(tcga.oscc.frac), names(nki.lpscc.frac))
identical(names(tcga.lpscc.frac), names(nki.lpscc.frac))
identical(names(tcga.lpscc.frac), names(nki.oscc.frac))



### Plot
pdf("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Figures/Sup_Figure_2B.pdf", family = "ArialMT", useDingbats=FALSE)

nki.oscc.frac.log = nki.oscc.frac
nki.oscc.frac.log[which(nki.oscc.frac.log == 0)] = 0.1
nki.oscc.frac.log = log10(nki.oscc.frac.log)

nki.lpscc.frac.log = nki.lpscc.frac
nki.lpscc.frac.log[which(nki.lpscc.frac == 0)] = 0.1
nki.lpscc.frac.log = log10(nki.lpscc.frac.log)

tcga.oscc.frac.log = log10(tcga.oscc.frac)

tcga.lpscc.frac.log = log10(tcga.lpscc.frac)

par(mar = c(5, 6, 3, 2), mgp = c(3.5, 1, 0))
plot(x = tcga.oscc.frac.log, 
     y = nki.oscc.frac.log,
     xlim = c(-1.1, 2), ylim = c(-1.1, 2),
     xlab = "% TCGA Aberrant", ylab = "% NKI Aberrant", main = expression("Genes"[CNA]), pch = 1,
     xaxs = "i", yaxs = "i", bty = "n",
     cex.main = 2, cex.lab = 3, cex.axis = 1.5, lwd = 1, las = 1, axes = F)
axis(side = 1, at = log10(c(0.1, 1, 10, 100)), 
     labels = c("0", "1", "10", "100"),
     cex.main = .6, cex.lab = 3, cex.axis = 1.5)
axis(side = 2, at = log10(c(0.1, 1, 10, 100)), 
     labels = c("0", "1", "10", "100"),
     cex.main = .6, cex.lab = 3, cex.axis = 1.5)
abline(a = 0, b = 1, lty = 2, lwd = 1)
points(x = tcga.lpscc.frac.log, 
       y = nki.lpscc.frac.log, pch = 2, lwd = 1)

legend("topleft", 
       c(paste("OSCC (cor = ", round(cor(tcga.oscc.frac, nki.oscc.frac, method = "spearman"), 2), ")", sep = ""),
         paste("L/P-SCC (cor = ", round(cor(tcga.lpscc.frac, nki.lpscc.frac, method = "spearman"), 2), ")", sep = "")), 
       pch = c(1, 2), bty = "n", cex = 1, pt.lwd = 1)

dev.off()



# Barplot of significant genes
pdf("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Figures/Figure_2D.pdf", family = "ArialMT", useDingbats=FALSE)

p.nki = matrix.fisher(nki.oscc.racg, nki.lpscc.racg)
idx.nki.significant = which(p.adjust(p.nki, method = "fdr") < .1)
colnames(nki.oscc.racg)[idx.nki.significant]
p.adjust(p.nki, method = "fdr")[idx.nki.significant]

p.tcga = matrix.fisher(tcga.oscc.racg.full, tcga.lpscc.racg.full)
idx.tcga.significant = which(p.adjust(p.tcga, method = "fdr") < .1)
colnames(tcga.oscc.racg.full)[idx.tcga.significant]
p.adjust(p.tcga, method = "fdr")[idx.tcga.significant]

names.significant = sort(unique(c(names(nki.oscc.frac)[idx.nki.significant], names(tcga.oscc.frac.full)[idx.tcga.significant])))

m.tcga = matrix(c(tcga.oscc.frac.full[names.significant], 0, 0, tcga.lpscc.frac.full[names.significant], 0, 0), byrow = T, nrow = 2)
m.nki = matrix(c(nki.oscc.frac[names.significant], 0, 0, nki.lpscc.frac[names.significant], 0, 0), byrow = T, nrow = 2)

par(mar = c(2, 5.5, 3, 0.5), cex = 1, mfrow = c(2, 1))

barplot.tcga = barplot(m.tcga, ylim = c(0, 80), names = c(names.significant, "", ""),
                       ylab = "% TCGA Aberrant", main = expression("Genes"[CNA]), cex.main = 3,
                       horiz = F, beside = T, cex.names = 2, cex.axis = 1.5, cex.lab = 2, yaxt = "n", las = 1) #xaxt = "n",
axis(side = 2, at = c(0, 20, 40, 60, 80), 
     labels = c(0, 20, 40, 60, 80), cex.axis = 1.5, cex.lab = 2, las = 1)
text(barplot.tcga, c(m.tcga) + 5, labels = paste(round(c(m.tcga), 0), "%", sep = ""), cex = 1.3)

legend("topleft", legend = c("OSCC", "L/P-SCC"), fill = gray.colors(2), bty = "n", cex = 1.3)

barplot.nki = barplot(-m.nki, ylim = c(-80, 0), 
                      ylab = "% NKI Aberrant", #names = names.significant, #main = "NKI",
                      horiz = F, yaxt = "n", beside = T, cex.axis = 1.5, cex.lab = 2, las = 1)
text(barplot.nki, -c(m.nki) - 5, labels = paste(round(c(m.nki), 0), "%", sep = ""), cex = 1.3)

axis(side = 2, at = -c(0, 20, 40, 60, 80), 
     labels = c(0, 20, 40, 60, 80), cex.axis = 1.5, cex.lab = 2, las = 1)
dev.off()


### Sup Table 5
# NKI total RACG/GenesCNA
x = c(rowSums(nki.oscc.racg), rowSums(nki.lpscc.racg))
y = c(rep("OSCC", nrow(nki.oscc.racg)), rep("L/P-SCC", nrow(nki.lpscc.racg)))
# Statistal test
wilcox.test(x ~ y)$p.val
# T-test only to get group means
t.test(x ~ y)

# NKI total RACG - excluding CASP8 tumors
oscc.casp8 = nki.oscc.mut[pmatch(rownames(nki.oscc.racg), rownames(nki.oscc.mut)), "CASP8"]
lpscc.casp8 = nki.lpscc.mut[pmatch(rownames(nki.lpscc.racg), rownames(nki.lpscc.mut)), "CASP8"]
nki.casp8 = c(oscc.casp8, lpscc.casp8)
nki.casp8.idx = which(nki.casp8 == 1)

wilcox.test(x[-nki.casp8.idx] ~ y[-nki.casp8.idx])
t.test(x[-nki.casp8.idx] ~ y[-nki.casp8.idx])

# TCGA total RACG
x = c(rowSums(tcga.oscc.racg.full), rowSums(tcga.lpscc.racg.full))
y = c(rep("OSCC", nrow(tcga.oscc.racg)), rep("L/P-SCC", nrow(tcga.lpscc.racg)))
wilcox.test(x ~ y)
t.test(x ~ y)

# TCGA total RACG - excluding CASP8 tumors
tcga.oscc.casp8 = tcga.mut.matrix[pmatch(rownames(tcga.oscc.racg.full), rownames(tcga.mut.matrix)), "CASP8"]
tcga.lpscc.casp8 = tcga.mut.matrix[pmatch(rownames(tcga.lpscc.racg.full), rownames(tcga.mut.matrix)), "CASP8"]
tcga.casp8 = c(tcga.oscc.casp8, tcga.lpscc.casp8)
tcga.casp8.idx = which(tcga.casp8 == 1)

wilcox.test(x[-tcga.casp8.idx] ~ y[-tcga.casp8.idx])
t.test(x[-tcga.casp8.idx] ~ y[-tcga.casp8.idx])