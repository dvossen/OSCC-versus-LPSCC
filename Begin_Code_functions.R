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

matrix.fisher = function(x, y){
  n = ncol(x)
  result = rep(NA, n)
  for(j in 1:n){
    m = matrix(c(length(which(x[,j] == 1)),
                 length(which(x[,j] == 0)),
                 length(which(y[,j] == 1)),
                 length(which(y[,j] == 0))), 2, 2)
    result[j] = fisher.test(m)$p.val
  }
  return(result)
}