fs = function(data, target, alpha = 0.05, max_l = 1, threshold = 0.01){
  
  cutoff <- qnorm(1 - alpha/2)
  n <- nrow(data)
  
  genes = setdiff(colnames(data), target)
  m = length(genes)*threshold
  y = data[,target]
  dm = data[,genes]
  data = cbind(y,dm)
  colnames(data)[1] = target
  
  C = cor(data)
  
  direct_cor = genes[abs(sapply(genes,ind_test,target,NULL,C,n)) > cutoff]
  print(length(direct_cor)) #####
  res = c()
  for (x in direct_cor){
    print(c(x,length(direct_cor) - which(direct_cor == x))) #####
    all_con_set = setdiff(direct_cor, x)
    count = 0
    l = 1
    con_set = NULL
    time = 0
    
    while (l <= max_l){
      con_set = nextSet(all_con_set, con_set, l)
      time = time + 1
      print(c('test:',time)) #####
      if ( is.null(con_set) ){
        l = l+1
      }
      else{
        if ( abs(ind_test(x, target, con_set, C, n)) <= cutoff ){
          count = count+1
          print(c('count =', count)) #####
          if (count >= m) break
        }
      }
    }
    
    if ( count < m ) res = c(res, x)
    
  }
  res
}