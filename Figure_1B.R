rm( list = ls( all = TRUE ) )

library(binom)

load("~/Documents/Projects/Article_Comparative_genomics_OSCC/Figures_CGs/Questions_Conchita/Mutation_heatmap/Iorio_genes.Rdata")
load("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Data/Figure_1A_data_TCGA.Rdata")
load("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Data/NKI_Mut_and_CNA.Rdata")

# Select the 27 GenesSPM that were captured
nki.oscc.mut = nki.oscc.mut[,hnscc.cg]
nki.lpscc.mut = nki.lpscc.mut[,hnscc.cg]
tcga.oscc.mut = tcga.oscc.mut[,hnscc.cg]
tcga.lpscc.mut = tcga.lpscc.mut[,hnscc.cg]

nki.oscc.mut.binom = binom.confint(x = colSums(nki.oscc.mut), n = nrow(nki.oscc.mut), methods = "wilson")
rownames(nki.oscc.mut.binom) = colnames(nki.oscc.mut)
nki.lpscc.mut.binom = binom.confint(x = colSums(nki.lpscc.mut), n = nrow(nki.lpscc.mut), methods = "wilson")
rownames(nki.lpscc.mut.binom) = colnames(nki.lpscc.mut)
tcga.oscc.mut.binom = binom.confint(x = colSums(tcga.oscc.mut), n = nrow(tcga.oscc.mut), methods = "wilson")
rownames(tcga.oscc.mut.binom) = colnames(tcga.oscc.mut)
tcga.lpscc.mut.binom = binom.confint(x = colSums(tcga.lpscc.mut), n = nrow(tcga.lpscc.mut), methods = "wilson")
rownames(tcga.lpscc.mut.binom) = colnames(tcga.lpscc.mut)

matrix.mut.mean = rbind(nki.oscc.mut.binom$mean, tcga.oscc.mut.binom$mean, nki.lpscc.mut.binom$mean, tcga.lpscc.mut.binom$mean) * 100
matrix.mut.lower = rbind(nki.oscc.mut.binom$lower, tcga.oscc.mut.binom$lower, nki.lpscc.mut.binom$lower, tcga.lpscc.mut.binom$lower) * 100
matrix.mut.upper = rbind(nki.oscc.mut.binom$upper, tcga.oscc.mut.binom$upper, nki.lpscc.mut.binom$upper, tcga.lpscc.mut.binom$upper) * 100
dimnames(matrix.mut.mean) = dimnames(matrix.mut.lower) = dimnames(matrix.mut.upper) = list(c("NKI OSCC", "TCGA OSCC", "NKI L/P-SCC", "TCGA L/P-SCC"), colnames(nki.oscc.mut))

idx.ordered = order(apply(matrix.mut.mean, 2, mean))
matrix.mut.mean.ordered = matrix.mut.mean[,idx.ordered]
matrix.mut.lower.ordered = matrix.mut.lower[,idx.ordered]
matrix.mut.upper.ordered = matrix.mut.upper[,idx.ordered]

idx.mean.5.pct = which(apply(matrix.mut.mean.ordered, 2, function(x){max(x) >= 5}))
idx.brca1 = which(colnames(matrix.mut.mean.ordered) == "BRCA1")
idx.mean.5.pct = sort(c(idx.mean.5.pct, idx.brca1))
matrix.mut.mean.ordered.min5pct = matrix.mut.mean.ordered[,idx.mean.5.pct]
matrix.mut.lower.ordered.min5pct = matrix.mut.lower.ordered[,idx.mean.5.pct]
matrix.mut.upper.ordered.min5pct = matrix.mut.upper.ordered[,idx.mean.5.pct]



pdf("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Figures/Figure_1B.pdf", family = "ArialMT", useDingbats=FALSE)
layout(matrix(1))
par(mar = c(5, 8, 3, 2), mgp = c(3.5, 1, 0))
barCenters = barplot(matrix.mut.mean.ordered.min5pct, beside = T, horiz = T, las = 1, 
        col = rep(grey.colors(2), each = 2),
        xlab = "Samples mutated (%)", main = expression("Genes"[SPM]),
        cex.names = 2, cex.axis = 1.5, cex.lab = 3, xlim = c(0, 100), cex.main = 3)
barplot(matrix.mut.mean.ordered.min5pct, beside = T, horiz = T, las = 1,
        col = c("white", "white", 1, 1), density = c(0, 10, 0, 10), add = T,
        cex.names = 2, cex.axis = 1.5, cex.lab = 3, xlim = c(0, 100))
segments(matrix.mut.lower.ordered.min5pct, barCenters, 
         matrix.mut.upper.ordered.min5pct, barCenters,
         lwd = 1)
legend("bottomright", legend = rownames(matrix.mut.mean.ordered.min5pct), bty = "n",
       fill = rep(grey.colors(2), each = 2), cex = 1.5)
legend("bottomright", legend = rownames(matrix.mut.mean.ordered.min5pct), bty = "n",
       col = c("white", "white", 1, 1), density = c(0, 10, 0, 10), cex = 1.5)
dev.off()