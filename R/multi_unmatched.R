


multi_unmatched = function(data,condition,n){ #In:Number of normal sample
  
  f11=function(mu,R=R,n=n,phij=phij,W.est=W.est,dataj=dataj){  #mu = c(alpha,(1...,I-1),muj,deltaj)
    ###par:
    mj= mu[1];
    deltaj = mu[-1];
    ######
    L22=0;
    yj=dataj[1:length(n[[1]])]
    mu1.phi = 1/(1 + mj * phij); ##vector of length I.
    L1 = 1/phij * log(mu1.phi) + yj * log(1-mu1.phi);
    for(i in 1:R){
      
      zj=dataj[n[[i+1]]]
      
      mu2.phi = 1/(1 + (mj + W.est[[i]]*deltaj[i]) * phij);
      
      
      L2 = 1/phij * log(mu2.phi) + zj * log(1-mu2.phi);
      L22=L22+sum(L2)
    }
    return(-sum(L1) - L22);
  }
  
  g11=function(mu,R=R,n=n,phij=phij,W.est=W.est,dataj=dataj){
    
    ###par:
    mj= mu[1];
    deltaj = mu[-1];
    yj=dataj[1:length(n[[1]])]
    #######
    mu1.phi = 1 + mj * phij; ##the inverse of the mu1.phi above in the f1 function.
    L1 = -1/mu1.phi + phij*yj/((mu1.phi-1)*mu1.phi);
    partial_deltaj=NULL;
    L22=0;
    for (i in 1:R)
    {
      zj=dataj[n[[i+1]]]
      
      mu2.phi = 1 + (mj+W.est[[i]]* deltaj[i]) * phij;
      
      
      L2 = -1/mu2.phi + phij*zj/((mu2.phi-1)*mu2.phi);
      L22=L22+sum(L2)
      
      partial_deltaj1 =  sum(-W.est[[i]]/mu2.phi+phij*W.est[[i]]*zj/((mu2.phi-1)*mu2.phi));
      partial_deltaj=cbind(partial_deltaj,partial_deltaj1)
    }
    partial_mj = sum(L1)+ L22;
    return (-c(partial_mj,partial_deltaj));
  }
  
  # minus loglikelihood function under null hypothesis
  f01=function(mu,R=R,n=n, phij=phij,dataj=dataj){  ##mu = muj
    ###par:
    mj= mu;
    #######
    yj=dataj[1:length(n[[1]])]
    mu.phi = 1/(1 + mj * phij);
    
    L1 = 1/phij * log(mu.phi) + yj * log(1-mu.phi);
    L22=0;
    for(i in 1:R)
    {
      zj=dataj[n[[i+1]]]
      L2 = 1/phij * log(mu.phi) + zj * log(1-mu.phi);
      L22=L22+sum(L2)
    }
    return(-sum(L1) - L22);
  }
  
  g01 = function(mu,R=R,n=n, phij=phij,dataj=dataj){
    ###par:
    mj= mu;
    #######
    yj=dataj[1:length(n[[1]])]
    
    mu.phi = 1 + mj * phij;  #the inverse of the mu.phi function in the f0 func above
    
    partial_mj=0
    for (i in 1:R)
    {
      zj=dataj[n[[i+1]]]
      partial_mj1 = sum(-1/mu.phi+ yj/(mu.phi*mj)) + sum(-1/mu.phi+ zj/(mu.phi*mj));
      partial_mj= partial_mj+ partial_mj1
    }
    return(-sum(partial_mj));
  }
   
  
  R=condition-1
  

  md = 'nloptr';
  
  
  data = as.matrix(data);
 
  nncol=rep(0,condition)
  for (j in 1:condition){
    nncol[j]=length(n[[j]])
  }   
  
  #  call edgeR first to estimate phi and find subset genes to estimate W.
  res.edgeR = edgeR_multi(data,condition,nncol=nncol,is.matched=F);
  Phi.est = res.edgeR$dispersion;

  data[data==0]=1;

  J = nrow(data);

  Y = data[,1:nncol[1]];  
  mu0 = rowMeans(Y);


  W.multi=NULL;
  W.est=NULL;

  W.max=NULL;
  Z=NULL;
  delta0=NULL;
  for(i in 1:R){
    
    Z=data[,n[[i+1]]]
    res=edgeR(cbind(Y,Z),In=nncol[1],is.matched=F)
    p.edgeR1=res$p.value
    FC.edgeR1=res$FC
    id = which(FC.edgeR1 > 0 & p.edgeR1 <1e-3);
  
    delta0[[i]] = rowMeans(Z) - mu0;
    dif = (Z[id,]-Y[id,]);
    W.estt = apply(dif,2,mean,trim=0.05);
    # normalized W so that they are identifiable
    W.est[[i]] = as.numeric(W.estt/mean(W.estt));
  
    W.max[[i]] = max( W.est[[i]]);
    
  }
  W1=list()
  W1=W.est;
  

  FC =matrix(NA,nrow=J,ncol=R)

  lrt =  p.new = rep(NA,J);
  t=1;#iteration time;
  Max.T=3;	

  opts = list("algorithm" = "NLOPT_LD_LBFGS", "xtol_rel"=1.0e-8);


  
	while(t< Max.T){
		### estimation of mu, delta, alpha; likelihood ratio test
	  if(t==1){
	   
	    mu   = rep(NA,J);    
	    delta  = matrix(NA,nrow=J,ncol=R)
	  }else{
	    W.est = W.est1;
	  }
    
		for(m in 1:J){

			phij = Phi.est[m];

      dataj=data[m,]
			# minus loglikelihood function under alternative hypothesis
	

			##### parameter estimation under H1
			mm = mu0[m];
      
      
			dd=rep(0,R)
			for (i in 1:R){
			  dd[i] =delta0[[i]][m]
			}
      
      for (i in 1:R){
  			if(dd[i]< 0){
  				if(-dd[i]*W.max[i] > mm){
  					dd[i] = - mm/W.max[i]*0.9;
  				}
  			}
      }
			par01=c(mm,dd);
			par00=mm;
      par01 = impute.par(par01); 
      par00 = impute.par(par00);

			res0 = nloptr(x0=par00, eval_f=f01,eval_grad_f=g01, opts=opts,R=R,n=n, phij=phij,dataj=dataj); # mle under null
			par0 = res0$solution;
			val0 = res0$objective;
			res11 = nloptr(x0=par01, eval_f=f11, eval_grad_f=g11, opts=opts,R=R,n=n,phij=phij,W.est=W.est,dataj=dataj); # mle under alternative
			val11 = res11$objective;
			res1 = res11;
		
			par1 = res1$solution;
			val1 = res1$objective;
		

			mu[m] = par1[1];
			delta[m,] = par1[-1];

			# likelihood ratio test statistic
			lrt[m] = 2*(val0-val1);    	
		}

		lrt[lrt < 0] = 1e-5;

		# p-value
		p.new = 1-pchisq(lrt,R);
		# fold change
		for (i in 1:R){
		  FC[,i]= log2((mu+delta[,i])/mu);
		  FC[is.na(FC)] = -100;
		}

    if(t==1){
    W.est1=list();
    for (i in 1:R){
      id = which(FC[,i] > 0 & p.new<1e-3);
      Z=data[,n[[i+1]]] 
      dif = Z[id,]-rowMeans(Y[id,]);	
      
      W.est1[[i]] = apply(dif,2,mean,trim=0.05);
      W.est1[[i]] = as.numeric(W.est1[[i]]/ mean( W.est1[[i]] ));
    
      W.max[[i]]= max(W.est1[[i]]);
        W.estfinal=list()
        W.estfinal=W.est1
      }
    }
		
		t = t+1;		
	}
  LR = cbind(lrt=lrt,p.value=p.new,logFC=FC);
  return (list(W=W.estfinal,LR=LR));
}
