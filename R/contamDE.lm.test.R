contamDE.lm.test <-
  function(d,H,b=NULL,covariate=NULL) {
   G <- nrow(d$counts)
   K <- ncol(d$log2FC)
   
   if(is.matrix(H)){
        r = nrow(H)
   } else {
       stop("H must be a numeric matrix")
   }
   if(is.null(b)){
     b <- rep(0,r)
   }else{
     if(!is.numeric(b) | length(b)!=r){
       stop('b should be a numeric vector with length = the row number of H')
     }
   }
   F.stat <-rep(0, G)
  
   for (g in 1:G){
    thetaghat <- as.matrix(d$log2FC[g,]);
    v.covghat <- matrix(d$log2FC.cov[g,],as.numeric(K),as.numeric(K));
    v.covghat <- matrix(d$log2FC.cov[g,],3,3);
    F.stat[g] <- t(H %*% thetaghat-b)%*%solve(H %*% v.covghat %*% t(H))%*%(H %*% thetaghat-b)/r
   }
   p.ftest <- 1 - pf(F.stat, df1 = r,df2 = d$df )
   d$p.contamDE.lm.test <- p.ftest

  return(d)
  }


