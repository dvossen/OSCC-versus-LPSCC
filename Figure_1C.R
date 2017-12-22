rm( list = ls( all = TRUE ) )

library(limma)
library(binom)

load("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Data/Iorio_genes.Rdata")

### Copy number alterations
# NKI OSCC
load("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Data/NKI_Mut_and_CNA.Rdata")
nki.oscc.gain = ifelse(nki.oscc.cn[,hnscc.racg.gain.captured] >= 7, 1, 0)
nki.oscc.loss = ifelse(nki.oscc.cn[,hnscc.racg.loss.captured] <= 0.5, 1, 0)
nki.oscc.racg = cbind(nki.oscc.gain, nki.oscc.loss)

nki.oscc.racg.binom = binom.confint(x = colSums(nki.oscc.racg), n = nrow(nki.oscc.racg), methods = "wilson")
rownames(nki.oscc.racg.binom) = colnames(nki.oscc.racg)

# NKI L/P-SCC
nki.lpscc.gain = ifelse(nki.lpscc.cn[,hnscc.racg.gain.captured] >= 7, 1, 0)
nki.lpscc.loss = ifelse(nki.lpscc.cn[,hnscc.racg.loss.captured] <= 0.5, 1, 0)
nki.lpscc.racg = cbind(nki.lpscc.gain, nki.lpscc.loss)

nki.lpscc.racg.binom = binom.confint(x = colSums(nki.lpscc.racg), n = nrow(nki.lpscc.racg), methods = "wilson")
rownames(nki.lpscc.racg.binom) = colnames(nki.lpscc.racg)

# TCGA OSCC
load("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Data/TCGA_Mut_and_CNA_small.Rdata")
tcga.oscc.cn = tcga.cna[which(rownames(tcga.cna) %in% tcga.clinical.data.hpv.neg$bcr_patient_barcode[which(tcga.clinical.data.hpv.neg$Location == "Oral Cavity")]),]

tcga.oscc.gain = ifelse(tcga.oscc.cn[,hnscc.racg.gain.captured] == 2, 1, 0)
tcga.oscc.loss = ifelse(tcga.oscc.cn[,hnscc.racg.loss.captured] == -2, 1, 0)
tcga.oscc.racg = cbind(tcga.oscc.gain, tcga.oscc.loss)

tcga.oscc.racg.binom = binom.confint(x = colSums(tcga.oscc.racg), n = nrow(tcga.oscc.racg), methods = "wilson")
rownames(tcga.oscc.racg.binom) = colnames(tcga.oscc.racg)

# TCGA L/P-SCC
tcga.lpscc.cn = tcga.cna[which(rownames(tcga.cna) %in% tcga.clinical.data.hpv.neg$bcr_patient_barcode[which(tcga.clinical.data.hpv.neg$Location != "Oral Cavity")]),]

tcga.lpscc.frac.gain = ifelse(tcga.lpscc.cn[,hnscc.racg.gain.captured] == 2, 1, 0)
tcga.lpscc.frac.loss = ifelse(tcga.lpscc.cn[,hnscc.racg.loss.captured] == -2, 1, 0)
tcga.lpscc.racg = cbind(tcga.lpscc.frac.gain, tcga.lpscc.frac.loss)

tcga.lpscc.racg.binom = binom.confint(x = colSums(tcga.lpscc.racg), n = nrow(tcga.lpscc.racg), methods = "wilson")
rownames(tcga.lpscc.racg.binom) = colnames(tcga.lpscc.racg)




matrix.racg.mean = rbind(nki.oscc.racg.binom$mean, tcga.oscc.racg.binom$mean, nki.lpscc.racg.binom$mean, tcga.lpscc.racg.binom$mean) * 100
matrix.racg.lower = rbind(nki.oscc.racg.binom$lower, tcga.oscc.racg.binom$lower, nki.lpscc.racg.binom$lower, tcga.lpscc.racg.binom$lower) * 100
matrix.racg.upper = rbind(nki.oscc.racg.binom$upper, tcga.oscc.racg.binom$upper, nki.lpscc.racg.binom$upper, tcga.lpscc.racg.binom$upper) * 100
dimnames(matrix.racg.mean) = dimnames(matrix.racg.lower) = dimnames(matrix.racg.upper) = list(c("NKI OSCC", "TCGA OSCC", "NKI L/P-SCC", "TCGA L/P-SCC"), colnames(nki.oscc.racg))

idx.ordered = order(apply(matrix.racg.mean, 2, mean))
matrix.racg.mean.ordered = matrix.racg.mean[,idx.ordered]
matrix.racg.lower.ordered = matrix.racg.lower[,idx.ordered]
matrix.racg.upper.ordered = matrix.racg.upper[,idx.ordered]

idx.mean.5.pct = which(apply(matrix.racg.mean.ordered, 2, function(x){max(x) >= 5}))
matrix.racg.mean.ordered.min5pct = matrix.racg.mean.ordered[,idx.mean.5.pct]
matrix.racg.lower.ordered.min5pct = matrix.racg.lower.ordered[,idx.mean.5.pct]
matrix.racg.upper.ordered.min5pct = matrix.racg.upper.ordered[,idx.mean.5.pct]



pdf("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Figures/Figure_1C.pdf", family = "ArialMT", useDingbats=FALSE)
par(mar = c(5, 8, 3, 2), mgp = c(3.5, 1, 0))
barCenters = barplot(matrix.racg.mean.ordered.min5pct, beside = T, horiz = T, las = 1, 
                     col = rep(grey.colors(2), each = 2),
                     xlab = "Samples aberrant (%)", main = expression("Genes"[CNA]),
                     cex.names = 2, cex.axis = 1.5, cex.lab = 3, xlim = c(0, 100), cex.main = 3)
barplot(matrix.racg.mean.ordered.min5pct, beside = T, horiz = T, las = 1,
        col = c("white", "white", 1, 1), density = c(0, 10, 0, 10), add = T,
        cex.names = 2, cex.axis = 1.5, cex.lab = 3, xlim = c(0, 100))
segments(matrix.racg.lower.ordered.min5pct, barCenters, 
         matrix.racg.upper.ordered.min5pct, barCenters,
         lwd = 1)
legend("bottomright", legend = rownames(matrix.racg.mean.ordered.min5pct), bty = "n",
       fill = rep(grey.colors(2), each = 2), cex = 1.5)
legend("bottomright", legend = rownames(matrix.racg.mean.ordered.min5pct), bty = "n",
       col = c("white", "white", 1, 1), density = c(0, 10, 0, 10), cex = 1.5)
dev.off()