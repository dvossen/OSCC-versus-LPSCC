rm( list = ls( all = TRUE ) )

library(ComplexHeatmap)

# Stuff for the plotting
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  },
  HOMDEL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "blue", col = NA))
  },
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "red", col = NA))
  },
  MUT = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "#008000", col = NA))
  }
)
col = c("MUT" = "#008000", "AMP" = "red", "HOMDEL" = "blue")

remove.zero.columns = function(x){
  x.col.sums = colSums(x)
  idx.remove = c(which(x.col.sums == 0), which(is.na(x.col.sums)))
  if(length(idx.remove) == 0){
    return(x)
  } else{
    return(x[,-idx.remove]) 
  }
}

# Actual data
load("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Data/Figure_1A_data_TCGA.Rdata")
load("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Data/NKI_Mut_and_CNA.Rdata")

all.selected.genes = unique(c("EGFR", "BRCA1", "HRAS", "ATRX", "ATM",  "ATR",  "NF1",  "FBXW7",  "APC",  
                       "CASP8",  "NOTCH1", "PIK3CA", "CDKN2A", "TP53",
                       "FANCG", "MDM2", "MYC", "EGFR", "TP63", "CDKN2A", "CCND1"))

# Plot combined data
combined.mut2 = remove.zero.columns(rbind(nki.oscc.mut[,all.selected.genes],
                                          nki.lpscc.mut[,all.selected.genes],
                                          tcga.oscc.mut[,all.selected.genes],
                                          tcga.lpscc.mut[,all.selected.genes]))
combined.loss2 = remove.zero.columns(rbind(nki.oscc.loss[,all.selected.genes],
                                           nki.lpscc.loss[,all.selected.genes],
                                           tcga.oscc.loss[,all.selected.genes],
                                           tcga.lpscc.loss[,all.selected.genes]))
combined.gain2 = remove.zero.columns(rbind(nki.oscc.gain[,all.selected.genes],
                                           nki.lpscc.gain[,all.selected.genes],
                                           tcga.oscc.gain[,all.selected.genes],
                                           tcga.lpscc.gain[,all.selected.genes]))

mat_list_combined = list(MUT = t(combined.mut2),
                         HOMDEL = t(combined.loss2),
                         AMP = t(combined.gain2))
mat_list_combined2 = unify_mat_list(mat_list_combined)

pdf("~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Figures/Figure_1A.pdf")
ha2 = HeatmapAnnotation(df = data.frame(Dataset = ifelse(grepl("TCGA", colnames(mat_list_combined2$MUT)), "TCGA", "NKI"),
                                        Subsite = c(rep("OSCC", nrow(nki.oscc.mut)), rep("L/P-SCC", nrow(nki.lpscc.mut)),
                                                    rep("OSCC", nrow(tcga.oscc.mut)), rep("L/P-SCC", nrow(tcga.lpscc.mut)))), 
                        col = list(Dataset = c("TCGA" =  "red", "NKI" = "blue"),
                                   Subsite = c("OSCC" =  gray.colors(2)[1], "L/P-SCC" = gray.colors(2)[2])))
oncoPrint(mat_list_combined2, alter_fun = alter_fun, col = col, bottom_annotation = ha2,
          show_row_barplot = F, top_annotation = NULL)
dev.off()