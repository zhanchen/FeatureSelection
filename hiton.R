new_hiton = function(data, target, max_k = 1 , alpha = 0.05){
  cutoff <- qnorm(1 - alpha/2)
  n <- nrow(data)
  
  genes = setdiff(colnames(data), target)
  y = data[,target]
  dm = data[,genes]
  data = cbind(y,dm)
  colnames(data)[1] = target
  
  C = cor(data)
  
  order_genes = sort(abs(sapply(genes,ind_test,target,NULL,C,n)), decreasing = TRUE)
  open = names(order_genes)[order_genes > cutoff]
  pc = c()
  
  while (length(open) > 0){
    pc = c(pc, open[1])
    l = 1
    con_set = NULL
    while (l <= max_k){
      con_set = nextSet(setdiff(pc, open[1]),con_set,l)
      if ( is.null(con_set) ){
        l = l+1
      }
      else{
        if ( abs(ind_test(open[1], target, con_set, C, n)) <= cutoff ){
          pc = setdiff(pc, open[1])
          break
        }
      }
    }
    open  = setdiff(open, open[1])
    print(c(length(open),length(pc)))
  }
  
  for ( x  in pc ){
    l = 1
    con_set = NULL
    print(c(x,length(pc)))
    while (l <= max_k){
      if ((which(pc == x))<length(pc)){
        con_set = nextSet(pc[(which(pc == x)+1):length(pc)],con_set,l)
        }
      if ( is.null(con_set) ){
        l = l+1
      }
      else{
        if ( abs(ind_test(x, target, con_set, C, n)) <= cutoff){
          pc = setdiff(pc, x)
          print(c('remove',x))
          break
        }
      }
    }
  }
  pc
}