contamDE <- 
function(data,R,n=NA,match=TRUE){
  
  data = as.matrix(data);

  
  if(!is.numeric(data[1])){
    stop('The input should be a numeric matrix');
  }
  
  if(match==FALSE){
    if(!is.list(n) | length(n[[1]]) <= 1){
      stop('The number of normal samples should be at least 2');
    }
    nncol=length(n[[1]])
    
    for (j in 2:R){
      nncol=nncol+length(n[[j]])
      if(length(n[[j]]) <= 1){
        stop('The number of tumor samples should be at least 2');
      }
    }  
    if(length(n)!=1){
      if(nncol!=ncol(data)){
        stop("The sum of the input list 'n' should be the same as the colomn number of 'data'");
      }
      else if(length(n)!=R){
        stop("The input list 'n' should have the same length as the condition number 'R'");
      }
    }
  }
  
  ms = as.numeric(apply(as.matrix(data), 2, "median"));
  m0 = exp(mean(log(ms)));
  data = sweep(data, 2, ms, "/") * m0;
  data = as.matrix(round(data));
  
  
  if(R==2){
    if(match){
      if(ncol(data)%%2 != 0){
        stop('Tumor samples should be paired with normal samples');
      }else{
        return(matched(data));
      }
    }else{
      In=length(n[[1]])
      return(unmatched(data,In));
    }
  }
  else{
    if(match){
      if(ncol(data)%%R != 0){
        stop('Tumor samples should be matched with normal samples');
      }else{
        return(multi_matched(data,R))
      }
    }
    else{
      return(multi_unmatched(data,R,n))
    }
  }
}
