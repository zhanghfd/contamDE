
matched<-
function(data){
   
	W.all=NULL;

  data = as.matrix(data);
  
  ##################### 
  # begin calling edgeR

  count = round(data);

  dm = dim(count);

  R = 2;
  I = dm[2]/2;

  condition = as.factor(rep(1:R,each=I));
  sample = as.factor(rep(1:I,R));

  d = DGEList(counts=count, group=condition);
  d = calcNormFactors(d);
  design = model.matrix(~condition+sample);
  d = estimateGLMCommonDisp(d, design);
  d = estimateGLMTagwiseDisp(d,design);

  fit = glmFit(d, design);

  result = glmLRT(fit, coef=2:R);

  p.edgeR = result$table$PValue;
  FC.edgeR = result$table$logFC;
  Phi.est = as.numeric(result$dispersion);
  
  # end calling edgeR
  #########################
  
	data[data==0]=1;
  I = ncol(data)/2;
  J = nrow(data);
  Y = data[,1:I];
  Z = data[,-(1:I)];

  mu0 = rowMeans(Y);
  delta0 = rowMeans(Z) - mu0;

  # pick out some most significantly down-regulated genes
  id=which(p.edgeR<1e-3&FC.edgeR < 0);

  # estimate W
  dif = (Z[id,]-Y[id,]);
  W.est = apply(dif,2,mean,trim=0.05);

  # normalized W so that they are identifiable
  W.est = as.numeric(W.est/mean(W.est));
  W.all = rbind(W.all, W.est);
  inv.W.min = 1/min(W.est);
  W.max = max(W.est);
  ###
	
	FC =  lrt =  p.new = rep(NA,J);
	t=1;#iteration time;
	Max.T=3;	

  ##begin iteration
	while( t< Max.T ){
    opts = list("algorithm" = "NLOPT_LD_LBFGS", "xtol_rel"=1.0e-8);
     
    if(t==1){
      alpha.m.new = alpha.m.new1 =  matrix(NA,nrow = J,ncol = I); #alpha.matrix
      mu = mu1 = delta = delta1 = rep(NA,J);        	
    }else{
    	alpha.m.new = alpha.m.new1;
    	mu = mu1;
      delta = delta1;
    	W.est = W.est1;
    }
    ### estimation of mu, delta, alpha; and do likelihood ratio test

    for(m in 1:J){
      phij = Phi.est[m];
      yj = Y[m,];
      zj = Z[m,];
    		
      ##### parameter estimation under H1
      mm = mu0[m];
      dd = delta0[m];
          
      if(dd < 0){
        if(-dd*W.max > mm){
          dd = - mm/W.max*0.9;
        }
      }
          
      par01=c(rep(1e-3,I-1),mm,dd);
      par00=c(rep(1e-3,I-1),mm);
		  par01 = impute.par(par01); 
      par00 = impute.par(par00);
		
      res0 = nloptr(x0=par00, eval_f=f0,eval_grad_f=g0, opts=opts, I=I, phij=phij,yj=yj,zj=zj); # mle under null
      par0 = res0$solution;
      val0 = res0$objective;
      par02 = c(par0,0);
			par02 = impute.par(par02);
			res11 = nloptr(x0=par01, eval_f=f1, eval_grad_f=g1,opts=opts,W.est=W.est,I=I,phij=phij,yj=yj,zj=zj); # mle under alternative 
			res12 = nloptr(x0=par02, eval_f=f1, eval_grad_f=g1,opts=opts,W.est=W.est,I=I,phij=phij,yj=yj,zj=zj); # mle under alternative
			    	
      val11 = res11$objective;
      val12 = res12$objective;

      if(val11 < val12){
        res1 = res11;
      }else{
        res1 = res12;
      }

      par1 = res1$solution;
      val1 = res1$objective;

      alpha.m.new1[m,] = c(par1[1:I-1],-sum(par1[1:I-1]));
      mu1[m] = par1[I];
      delta1[m] = par1[I+1];
      lrt[m] = 2*(val0-val1);    	
    }
    	
    lrt[lrt < 0] = 1e-5;
    p.new = 1-pchisq(lrt,1);

    # fold change
    FC = log2((mu1+delta1)/mu1);
  	FC[is.na(FC)] = -100;
    id=which(p.new<1e-3&FC<0);
	  dif = Z[id,]-Y[id,];	

    #update the estimated value of W
		denomitor = delta1[id]*exp(alpha.m.new1[id,]);
  	dif = apply(dif,2,mean,trim = 0.05);
	  denomitor = apply(denomitor,2,mean,trim = 0.05);
	  W.est1 = dif/denomitor;		
	  W.est1 = as.numeric(W.est1/ mean( W.est1 ));
		W.est = W.est1;
	  inv.W.min = 1/min(W.est1);
	  W.max = max(W.est1);
    W.all = rbind(W.all, W.est1);

		t = t+1;
  }
  ##end iteration
  LR = cbind(lrt=lrt,p.value=p.new,logFC=FC);
	return(list(W=W.all[2,],LR=LR));
}
