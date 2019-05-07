


multi_matched = function(data,condition){
  
  f12=function(mu,R,I,phij,W.est,dataj){  #mu = c(alpha,(1...,I-1),muj,deltaj)
    ###par:
    alpha = c(mu[1:(I-1)],-sum(mu[1:(I-1)]));
    e.alpha = exp(alpha);
    mj= mu[I];
    deltaj = mu[-(1:I)];
    mu1.phi = 1/(1 + mj*e.alpha * phij); ##vector of length I.
    L1 = 1/phij * log(mu1.phi) +  dataj[1:I] * log(1-mu1.phi);
    
    L22=0;
    for (i in 1:R){
      mu2.phi = 1/(1 + (mj + W.est[[i]]*deltaj[i]) * e.alpha * phij);  
      L2 = 1/phij * log(mu2.phi) +  dataj[(i*I+1):((i+1)*I)] * log(1-mu2.phi);
      L22=L22+L2
    }
    return(-sum(L1+L22));
    
    
  }
  
  g12=function(mu,R,I,phij,W.est,dataj){
    
    yj = dataj[1:I];  
    yj1 = yj[1:(I-1)];
    yj2 = yj[I];
    #   partial_alpha=NULL;
    #   partial_mj=NULL;
    partial_deltaj=NULL;
    L2=0;
    L4=0;
    alpha = c(mu[1:(I-1)],-sum(mu[1:(I-1)]));
    e.alpha =  exp(c(mu[1:(I-1)],-sum(mu[1:(I-1)])));
    mj= mu[I];
    deltaj = mu[-(1:I)];
    mu1.phi = 1 + mj * phij*e.alpha; ##the inverse of the mu1.phi above in the f1 function.
    mu1.phi1 = mu1.phi[1:(I-1)]; mu1.phi2 = mu1.phi[I];
    mu1.alpha = mj*e.alpha;
    mu1.alpha1 = mu1.alpha[1:(I-1)];mu1.alpha2 = mu1.alpha[I];
    e.alpha1 = e.alpha[1:(I-1)]; e.alpha2 = e.alpha[I];
    
    L1 = -mu1.alpha1/mu1.phi1 + mu1.alpha2/mu1.phi2 + yj1/mu1.phi1 - yj2/mu1.phi2;
    L3 = -e.alpha/mu1.phi + phij*yj*e.alpha/((mu1.phi-1)*mu1.phi);
    
    for(i in 1:R){  
      zj=dataj[(i*I+1):((i+1)*I)] 
      zj1 = zj[1:(I-1)];
      zj2 = zj[I];
      ###par:
      
      #######
      
      mu2.phi = 1 + (mj+W.est[[i]]* deltaj[i]) * phij*e.alpha;
      
      mu2.phi1 = mu2.phi[1:(I-1)]; mu2.phi2 = mu2.phi[I];
      
      
      mu2.alpha = (mj+W.est[[i]]*deltaj[i])*e.alpha;
      
      mu2.alpha1 = mu2.alpha[1:(I-1)];mu2.alpha2 = mu2.alpha[I];
            
      L22 = -mu2.alpha1/mu2.phi1 + mu2.alpha2/mu2.phi2 + zj1/mu2.phi1 - zj2/mu2.phi2;
      
      L2 =L2+L22;  #vector of length I-1
      
      
      L44 = -e.alpha/mu2.phi + phij*zj*e.alpha/((mu2.phi-1)*mu2.phi);
      L4=L4+L44
      
      partial_deltaj1 =  sum(-W.est[[i]]*e.alpha/mu2.phi+phij*W.est[[i]]*zj*e.alpha/((mu2.phi-1)*mu2.phi));
      partial_deltaj=c(partial_deltaj,partial_deltaj1)
      
      
    }
    
    partial_alpha=L1+L2;
    partial_mj = sum(L3+L4);
    
    return (-c(partial_alpha,partial_mj,partial_deltaj));  
    
    
  }
  
  # minus loglikelihood function under null hypothesis
  
  f02=function(mu,R,I,phij,dataj){  ##mu = c(alpha,(1,...,I-1),muj)
    ###par:
    alpha = c(mu[1:(I-1)],-sum(mu[1:(I-1)]));
    e.alpha = exp(alpha);
    mj= mu[I];
    mu.phi = 1/(1 + e.alpha * mj * phij);
    L1 = 1/phij * log(mu.phi) + dataj[1:I] * log(1-mu.phi);
    L22=0;
    
    for (i in 1:R){      
      L2 = 1/phij * log(mu.phi) + dataj[(i*I+1):((i+1)*I)] * log(1-mu.phi);
      L22=L22+L2
    }
    return(-sum(L1+L22));
  }
  
  g02 = function(mu,R,I,phij,dataj){
    yj = dataj[1:I];  
    yj1 = yj[1:(I-1)];
    yj2 = yj[I];
    partial_alpha=0;
    partial_mj=0;
    alpha = c(mu[1:(I-1)],-sum(mu[1:(I-1)]));
    e.alpha = exp(alpha);
    mj= mu[I];
    mu.phi = 1 + mj * phij * e.alpha;  #the inverse of the mu.phi function in the f0 func above
    mu.phi1 = mu.phi[1:(I-1)];mu.phi2 = mu.phi[I];
    
    mu.alpha = mj*e.alpha;
    mu.alpha1 = mu.alpha[1:(I-1)];mu.alpha2 = mu.alpha[I];
    ####
    sumj=sumj1=sumj2=0
    for(i in 1:R){  
      zj=dataj[(i*I+1):((i+1)*I)] 
      zj1 = zj[1:(I-1)];
      zj2 = zj[I];
      
      sumj=sumj+zj
      sumj1=sumj1+zj1
      sumj2=sumj2+zj2
    }
    partial_alpha = -(R+1)*mu.alpha1/mu.phi1+(R+1) * mu.alpha2 / mu.phi2 + (sumj1+yj1)/mu.phi1 - (yj2+sumj2)/mu.phi2;  #length of I-1
    
    partial_mj = sum(-(R+1)*e.alpha/mu.phi+(yj+sumj)/(mu.phi*mj));
    
    return(-c(partial_alpha,partial_mj));
  }
  R=condition-1
  data = as.matrix(data);
	
	res.edgeR = edgeR_multi(data,(R+1));
	Phi.est = res.edgeR$dispersion;

	data[data==0]=1;
  
  I = ncol(data)/(R+1);   
  J = nrow(data);

	Y = data[,1:I];  
	mu0 = rowMeans(Y);


	W.est=NULL;
  W.max=NULL;
  Z=NULL;
  delta0=NULL;

  for(i in 1:R){
    Z=data[,(i*I+1):((i+1)*I)]
    res=edgeR(cbind(Y,Z),I,is.matched=T)
    p.edgeR1=res$p.value
    FC.edgeR1=res$FC
    id= which(FC.edgeR1 > 0 & p.edgeR1<1e-3);
    
    delta0[[i]] = rowMeans(Z) - mu0;
    dif = (Z[id,]-Y[id,]);
    W.estt = apply(dif,2,mean,trim=0.05);
    # normalized W so that they are identifiable
    W.est[[i]] = as.numeric(W.estt/mean(W.estt));
    W.max[[i]] = max(W.est[[i]])

  }
  
  
  W1=list()
  W1=W.est
  FC =matrix(NA,nrow=J,ncol=R)
	
  lrt =  p.new = rep(NA,J);
	t=1;#iteration time;
	Max.T=3;	
	while( t< Max.T ){
    opts = list("algorithm" = "NLOPT_LD_LBFGS", "xtol_rel"=1.0e-8);
    if(t==1){
      alpha.m.new1 =  matrix(NA,nrow = J,ncol = I); #alpha.matrix
      mu1  = rep(NA,J);    
      delta1= matrix(NA,nrow=J,ncol=R)
    }else{
    		W.est = W.est1;
    }

    for(m in 1:J){
      phij = Phi.est[m];
      dataj = data[m,];
      mm = mu0[m];
      ##### parameter estimation under H1
   
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

      par01=c(rep(1e-3,I-1),mm,dd);
      par00=c(rep(1e-3,I-1),mm);
			par01 = impute.par(par01); 
      par00 = impute.par(par00);
		
      res0 = nloptr(x0=par00, eval_f=f02,eval_grad_f=g02, opts=opts,R=R,I=I, phij=phij,dataj=dataj); # mle under null
      par0 = res0$solution;
      val0 = res0$objective;
      res11 = nloptr(x0=par01, eval_f=f12, eval_grad_f=g12, opts=opts,R=R,I=I,phij=phij,W.est=W.est,dataj=dataj); # mle under alternative
      res1 = res11;
      par1 = res1$solution;
      val1 = res1$objective;

    	alpha.m.new1[m,] = c(par1[1:I-1],-sum(par1[1:I-1]));
      mu1[m] = par1[I];
      delta1[m,] = par1[-(1:I)];

      lrt[m] = 2*(val0-val1);    	
    }
    	

    ###############	
    lrt[lrt < 0] = 1e-5;
    # p-value
    p.new = 1-pchisq(lrt,R);
    if(t==1){
      p.new0=p.new
    }
    # fold change
    for (i in 1:R){
      FC[,i]= log2((mu1+delta1[,i])/mu1);
      FC[is.na(FC)] = -100;
    }

    if (t==1){
      W.est1=list();
      for (i in 1:R){
        id = which(FC[,i] > 0 & p.new<1e-3);
        Z=data[,(i*I+1):((i+1)*I)] 
    		dif = Z[id,]-Y[id,];	
    		#update the estimated value of W
    		denomitor = delta1[id,R]*exp(alpha.m.new1[id,]);
    		dif = apply(dif,2,mean,trim = 0.05);
    		denomitor = apply(denomitor,2,mean,trim = 0.05);
    		W.est1[[i]] = dif/denomitor;		
    		W.est1[[i]]= as.numeric(W.est1[[i]]/ mean( W.est1[[i]] ));
        W.max[[i]] = max(W.est1[[i]])
        W.estfinal=list()
        W.estfinal=W.est1
      }
    }	
    t = t+1;
	}	
  LR = cbind(lrt=lrt,p.value=p.new,logFC=FC);
  return (list (W=W.estfinal,LR=LR));
}

