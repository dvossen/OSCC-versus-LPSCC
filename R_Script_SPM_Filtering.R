rm( list = ls( all = TRUE ) )

# Define functions
filter.varscan.results <- function(input.varscan.results = NULL, input.sample.name = NULL,
                                   # Select only variants in particular genes
                                   input.vector.target.genes = NULL, remove.non.target.genes = FALSE,
                                   # Filters aimed at removing germline variants
                                   remove.in.dbsnp.not.in.cosm = FALSE,
                                   custom.snp.database = NULL, custom.snp.cutoff = NULL,
                                   # Filters aimed at removing variant unlikely to have a functional effect
                                   remove.non.exonic.splicing = FALSE,
                                   remove.synonymous.variants = FALSE, remove.splice.variants.not.in.cosm = FALSE,
                                   # Filters aimed at removing possible false positives
                                   remove.minimal.fraction = FALSE, minimal.fraction = NULL,
                                   # Somatic option
                                   somatic.annotation = FALSE){
  
  # Make variables that indicate if the variant is annotated in dbSNP, COSM, and ExAC_NFE
  input.varscan.results = input.varscan.results
  input.varscan.results$in.dbsnp142                = !is.na(input.varscan.results$snp142)
  input.varscan.results$in.cosmic70                = !is.na(input.varscan.results$cosmic70)
  input.varscan.results$n.in.cosmic70              = unname(sapply(input.varscan.results$cosmic70, function.return.n.cosmic.entries))
  input.varscan.results$n.in.cosmic70[which(is.na(input.varscan.results$n.in.cosmic70))] = 0
  
  input.varscan.results$sample = rep(input.sample.name, nrow(input.varscan.results))
  input.varscan.results$annovar.gene = unname(sapply(input.varscan.results$Gene.refGene, function(x){unlist(unname(strsplit(x, "(", fixed = T)))[1]}))
  input.varscan.results$sample.vaf = as.numeric(gsub("%", "", unname(sapply(input.varscan.results$Otherinfo, function(x){unlist(strsplit(unlist(strsplit(x, "\t"))[10], ":"))[7]}))))
  
  if(remove.non.exonic.splicing)
  {
    positions.exonic.or.splicing = grep("^exonic$|^splicing$|^exonic;splicing$", input.varscan.results$Func.refGene)
    
    if(length(positions.exonic.or.splicing) != 0)
    {
      input.varscan.results = input.varscan.results[positions.exonic.or.splicing,] 
    } else{
      stop("Input does not contain exonic or splicing variants")
    }
  }
  
  if(remove.non.target.genes)
  {
    positions.on.target             = which(input.varscan.results$annovar.gene %in% input.vector.target.genes)
    positions.off.target            = which(!input.varscan.results$annovar.gene %in% input.vector.target.genes)
    
    if(length(positions.off.target) != 0)
    {
      input.varscan.results                  = input.varscan.results[-positions.off.target,]
    }
    
    if(length(positions.on.target) == 0)
    {
      return(input.varscan.results)
    }
  }
  
  # Remove variants in dbSNP 142, except those also in COSMIC 70
  if(remove.in.dbsnp.not.in.cosm)
  {
    positions.in.dbsnp.not.in.cosm      = which(input.varscan.results$in.dbsnp142 & !input.varscan.results$in.cosmic70)
    
    if(length(positions.in.dbsnp.not.in.cosm) != 0)
    {
      input.varscan.results                  = input.varscan.results[-positions.in.dbsnp.not.in.cosm,]
    }
  }
  
  # Remove variants whose MAF exceeds the maximum MAF
  if(!is.null(custom.snp.database) & !is.null(custom.snp.cutoff))
  {
    positions.too.high.maf            = which(as.numeric(input.varscan.results[,custom.snp.database]) > custom.snp.cutoff)
    
    if(length(positions.too.high.maf) != 0)
    {
      input.varscan.results                  = input.varscan.results[-positions.too.high.maf,]
    }
  }
  
  if(remove.synonymous.variants)
  {
    positions.synonymous                = which(input.varscan.results$ExonicFunc.refGene == "synonymous SNV")
    
    if(length(positions.synonymous) != 0)
    {
      input.varscan.results                  = input.varscan.results[-positions.synonymous,]
    }
  }
  
  if(remove.splice.variants.not.in.cosm)
  {
    position.splicing.not.in.cosm      = which(input.varscan.results$Func.refGene == "splicing" & !input.varscan.results$in.cosmic70)
    
    if(length(position.splicing.not.in.cosm) != 0)
    {
      input.varscan.results                  = input.varscan.results[-position.splicing.not.in.cosm,]
    }
  }
  
  if(remove.minimal.fraction)
  {
    # Remove variants with a fraction < X (e.g., 20%)
    position.too.low.fraction         = which(input.varscan.results$sample.vaf < minimal.fraction)
    
    if(length(position.too.low.fraction) != 0)
    {
      input.varscan.results                  = input.varscan.results[-position.too.low.fraction,]
    }
  }
  
  if(somatic.annotation)
  {
    input.varscan.results$prediction = rep("germline", nrow(input.varscan.results))
    input.varscan.results$prediction[which(is.na(input.varscan.results$X1000g2015aug_eur) & is.na(input.varscan.results$esp6500siv2_ea) & is.na(input.varscan.results$ExAC_NFE))] = "somatic"
    input.varscan.results$prediction[which(input.varscan.results$ExonicFunc.refGene == "unknown")] = "germline"
    input.varscan.results$prediction[which(input.varscan.results$n.in.cosmic70 > 20)] = "somatic"
  }
  
  return(input.varscan.results)
}

function.return.n.cosmic.entries = function(x){
  if(!is.na(x)){
    x = unlist(strsplit(x, "OCCURENCE="))[2]
    x = unlist(strsplit(x, ",", fixed = T))
    x = sum(as.numeric(unname(sapply(x, function(x){unlist(strsplit(x, "(", fixed = T))[1]}))))
    return(x) 
  } else{
    return(NA)
  }
}

return.binary.varscan.variant.matrix = function(input.df.variants = NULL, input.vector.samples = NULL, input.vector.target.genes = NULL){
  n.samples = length(input.vector.samples)
  n.genes = length(input.vector.target.genes)
  
  result.matrix = matrix(0, nrow = n.samples, ncol = n.genes, 
                         dimnames = list(input.vector.samples, input.vector.target.genes))
  
  for(i in 1:n.samples){
    temp.sample = input.vector.samples[i]
    temp.sample.idx = grep(temp.sample, input.df.variants$sample)
    
    if(length(temp.sample.idx) == 0){
      next
    } else{
      temp.sample.variants = input.df.variants[temp.sample.idx,]
      
      result.matrix[i,] = ifelse(input.vector.target.genes %in% temp.sample.variants$annovar.gene, 1, 0)
    }
  }
  
  return(result.matrix)
}

# Load data
bed.file = read.table("~/...bed",
                      header = F, stringsAsFactors = F, sep = "\t")
targeted.genes = sort(unique(bed.file$V5))

varscan.files = dir("~/...", pattern = ".hg19_multianno.csv", full.names = TRUE, ignore.case = TRUE)

# Start
mondholte.file.names = NULL
for(i in 1:length(varscan.files)){
  temp.file.name = unlist(strsplit(basename(varscan.files[i]), "-uniq", fixed = T))[1]
  mondholte.file.names = c(mondholte.file.names, temp.file.name)
}
mondholte.file.names = unique(mondholte.file.names)

mondholte.variants = NULL
for(i in 1:length(varscan.files)){
  temp.file = read.csv(file = varscan.files[i], header = T, stringsAsFactors = F)
  temp.file.name = unlist(strsplit(basename(varscan.files[i]), "-uniq", fixed = T))[1]
  
  temp.file = filter.varscan.results(input.varscan.results = temp.file, input.sample.name = temp.file.name,
                                     input.vector.target.genes = targeted.genes, remove.non.target.genes = T,
                                     # Already remove SNPs to save space
                                     custom.snp.database = 'X1000g2015aug_eur', custom.snp.cutoff = 0.01,
                                     remove.synonymous.variants = T, remove.non.exonic.splicing = T, 
                                     # Annotate with somatic prediction
                                     somatic.annotation = T)
  
  mondholte.variants = rbind(mondholte.variants, temp.file)
}
mondholte.variants.somatic = mondholte.variants[which(mondholte.variants$prediction == "somatic"),]

# Binary variant matrix
binary.matrix.mondholte.variants.somatic = return.binary.varscan.variant.matrix(input.df.variants = mondholte.variants.somatic, 
                                                                                input.vector.target.genes = targeted.genes, 
                                                                                input.vector.samples = sort(unique(mondholte.variants.somatic$sample)))

missing.samples = mondholte.file.names[which(!mondholte.file.names %in% sort(unique(mondholte.variants.somatic$sample)))]
n.missing.samples = length(missing.samples)
if(n.missing.samples != 0){
  n.missing.samples.matrix = matrix(0, nrow = n.missing.samples, ncol = length(targeted.genes), dimnames = list(missing.samples, targeted.genes))
  binary.matrix.mondholte.variants.somatic = rbind(binary.matrix.mondholte.variants.somatic, n.missing.samples.matrix)
}
rownames(binary.matrix.mondholte.variants.somatic) = substring(rownames(binary.matrix.mondholte.variants.somatic), 6, 13)

save(mondholte.variants, mondholte.variants.somatic, binary.matrix.mondholte.variants.somatic,
     file = "~/...Rdata")