rm( list = ls( all = TRUE ) )

library(limma)

### Recurrently mutated and captured genes
# Iorio et al., 2016 cancer genes for which mutation patern in WES is consistent with positive selection in HNSCC
iorio.s2a = read.csv(file = "~/Documents/Projects/Capture_set/Pan-cancer/Iorio/S2A.csv", header = T, stringsAsFactors = F)

iorio.hnscc.cancer.genes = sort(unique(iorio.s2a$Gene[which(iorio.s2a$Cancer.Type == "HNSC")]))
iorio.hnscc.cancer.genes = gsub(pattern = "*", replacement = "", x = iorio.hnscc.cancer.genes, fixed = T)
iorio.hnscc.cancer.genes = unique(alias2Symbol(iorio.hnscc.cancer.genes))

### Cancer genes recurrently gained / lossed (overlap Iorio HNSCC RACS and Cancer Census genes) and captured genes
# Iorio et al., 2016 cancer genes for which mutation patern in WES is consistent with positive selection in HNSCC
load("~/Documents/Projects/Capture_set/Pan-cancer/Iorio/RACS/Iorio_HNSCC_RACSs_CensusGenes.Rdata")
hnscc.racg.loss = unique(alias2Symbol(as.character(iorio.hnscc.racg$Gene[which(iorio.hnscc.racg$Type == "Deletion")])))
hnscc.racg.gain = unique(alias2Symbol(as.character(iorio.hnscc.racg$Gene[which(iorio.hnscc.racg$Type == "Amplification")])))

# Captured genes
bed.file = read.table("~/Documents/Projects/General_Data/BED_file/fanc-baits-flat-annotated_11_04_2017.bed",
                      header = F, stringsAsFactors = F, sep = "\t")
targeted.genes = unique(alias2Symbol(sort(unique(bed.file$V5))))

hnscc.cg = targeted.genes[which(targeted.genes %in% iorio.hnscc.cancer.genes)]
hnscc.racg.loss.captured = targeted.genes[which(targeted.genes %in% hnscc.racg.loss)]
hnscc.racg.gain.captured = targeted.genes[which(targeted.genes %in% hnscc.racg.gain)]

all.iorio.genes = c(hnscc.cg, hnscc.racg.loss.captured, hnscc.racg.gain.captured)

save(all.iorio.genes, 
     iorio.hnscc.cancer.genes, hnscc.cg, 
     hnscc.racg.gain, hnscc.racg.loss,
     hnscc.racg.loss.captured, hnscc.racg.gain.captured,
     file = "~/Documents/Projects/Article_Comparative_genomics_OSCC/Manuscript_code/Data/Iorio_genes.Rdata")