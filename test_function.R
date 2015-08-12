getr = function (i, j, k, C, cut.at = 0.9999999) 
{
  if (length(k) == 0) {
    r <- C[i, j]
  }
  else if (length(k) == 1) {
    r <- (C[i, j] - C[i, k] * C[j, k])/sqrt((1 - C[j, k]^2) * 
                                              (1 - C[i, k]^2))
  }
  else {
    PM <- MASS::ginv(C[c(i, j, k), c(i, j, k)])
    r <- -PM[1, 2]/sqrt(PM[1, 1] * PM[2, 2])
  }
  if (is.na(r)) 
    0
  else min(cut.at, max(-cut.at, r))
}
log.q1pm <- function(r) log1p(2*r/(1-r))
ind_test = function (x, y, S, C, n) 
{
  r <- getr(x, y, S, C)
  res <- sqrt(n - length(S) - 3) * 0.5 * log.q1pm(r)
  if (is.na(res)) 
    0
  else res
}
nextSet = function(all_var, cur_set = NULL, l = -1){
  next_set = NULL
  if ( length(all_var) > 0 ){
    if (is.null(cur_set)){
      if ( l > 0 & l <= length(all_var)){
        next_set = all_var[c(1:l)]
      }
    }
    else{
      l = length(cur_set)
      index = c()
      for (e in cur_set) index = c(index, which(all_var == e))
      end = length(all_var)
      index[l] = index[l]+1
      for (i in (l+1 - 1:l)){
        if (index[i]>end){
          if (i > 1){
            index[i] = index[i]%%end
            index[i-1] = index[i-1]+1
          }
          else{
            index = NULL
          }
        }
        else{
          break
        }
      }
      if (length(index)>0) next_set = all_var[index]
    }
  }
  next_set
}
