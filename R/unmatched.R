unmatched <-
function(data,In=In){ #In:Number of normal sample

	W.all = NULL;

  data = as.matrix(data);

  In = In;
	Ic = ncol(data)-In;
  J = nrow(data);


  Y = data[,1:In];
  Z = data[,-(1:In)];

  mu0 = rowMeans(Y);
  delta0 = rowMeans(Z) - mu0;

  ###
  ## estimate phi using edgeR
  dge = DGEList(counts=data);
  dge = calcNormFactors(dge);

  Tissue = c(rep(1,In),rep(2,Ic));
  design = model.matrix(~Tissue); #unmatched design

  dge = estimateCommonDisp(dge,design); #estimate dispersion
  dge = estimateGLMTagwiseDisp(dge,design);

  fit = glmFit(dge,design);
  fit = glmLRT(fit,2);
  Phi.est = fit$dispersion;
  ###
  p.edgeR = fit$table$PValue;  ##report edgeR's p-value.

  # pick out some most significantly down-regulated genes
	id = which(fit$table$logFC < 0 & p.edgeR<1e-3);

  # estimate W
  dif = (Z[id,]-rowMeans(Y[id,]));
  W.est = apply(dif,2,mean,trim=0.05);
  # normalized W so that they are identifiable
  W.est = as.numeric(W.est/mean(W.est));

  W.all = rbind(W.all,W.est);
  inv.W.min = 1/min(W.est);
  W.max = max(W.est);
  ###
  data[data==0] = 1;
  opts = list("algorithm" = "NLOPT_LD_LBFGS", "xtol_rel"=1.0e-8);

  lrt = mu = delta = rep(NA,J);
	t=1;Max.T = 3;
	while(t< Max.T){
		###
		### estimation of mu, delta, alpha; likelihood ratio test
		for(m in 1:J){
			phij = Phi.est[m];
			yj = data[m,1:In];
			zj = data[m,-(1:In)];

			# minus loglikelihood function under alternative hypothesis
			f1=function(mu){  #mu = c(alpha,(1...,I-1),muj,deltaj)
				###par:
				mj= mu[1];
				deltaj = mu[2];
				######

				mu1.phi = 1/(1 + mj * phij); ##vector of length I.
				mu2.phi = 1/(1 + (mj + W.est*deltaj) * phij);

				L1 = 1/phij * log(mu1.phi) + yj * log(1-mu1.phi);
				L2 = 1/phij * log(mu2.phi) + zj * log(1-mu2.phi);
				return(-sum(L1) - sum(L2));
			}

			g1=function(mu){

				###par:
				mj= mu[1];
				deltaj = mu[2];
				#######

				mu1.phi = 1 + mj * phij; ##the inverse of the mu1.phi above in the f1 function.
				mu2.phi = 1 + (mj+W.est* deltaj) * phij;

				L1 = -1/mu1.phi + phij*yj/((mu1.phi-1)*mu1.phi);
				L2 = -1/mu2.phi + phij*zj/((mu2.phi-1)*mu2.phi);
				partial_mj = sum(L1)+ sum(L2);

				partial_deltaj =  sum(-W.est/mu2.phi+phij*W.est*zj/((mu2.phi-1)*mu2.phi));
				return (-c(partial_mj,partial_deltaj));
			}

			# minus loglikelihood function under null hypothesis
			f0=function(mu){  ##mu = muj
				###par:
				mj= mu;
				#######

				mu.phi = 1/(1 + mj * phij);

				L1 = 1/phij * log(mu.phi) + yj * log(1-mu.phi);
				L2 = 1/phij * log(mu.phi) + zj * log(1-mu.phi);

				return(-sum(L1) - sum(L2));
			}

			g0 = function(mu){
				###par:
				mj= mu;
				#######

				mu.phi = 1 + mj * phij;  #the inverse of the mu.phi function in the f0 func above

				partial_mj = sum(-1/mu.phi+ yj/(mu.phi*mj)) + sum(-1/mu.phi+ zj/(mu.phi*mj));

				return(-partial_mj);
			}

			##### parameter estimation under H1
			mm = mean(data[m,1:In]);
			dd = delta0[m];
			if(dd < 0){
				if(-dd*W.max > mm){
					dd = - mm/W.max*0.9;
				}
			}

			par01=c(mm,dd);
			par00=mm;

      res0 = nloptr(x0=par00, eval_f=f0,eval_grad_f=g0,opts=opts); # mle under null
			par0 = res0$solution;
			val0 = res0$objective;
			par02 = c(par0,0);
			res11 = nloptr(x0=par01, eval_f=f1, eval_grad_f=g1,opts=opts); # mle under alternative
			res12 = nloptr(x0=par02, eval_f=f1, eval_grad_f=g1,opts=opts); # mle under alternative

      val11 = res11$objective;
			val12 = res12$objective;

			if(val11 < val12){
				res1 = res11;
			}else{
				res1 = res12;
			}

			par1 = res1$solution;
			val1 = res1$objective;

			mu[m] = par1[1];
			delta[m] = par1[2];

			# likelihood ratio test statistic
			lrt[m] = 2*(val0-val1);
			if(lrt[m] < 0) cat(m,lrt[m],'\n');
		}

		lrt[lrt < 0] = 1e-5;

		# p-value
		p.new = 1-pchisq(lrt,1);
		# fold change
		FC = log2((mu+delta)/mu);
		FC[is.na(FC)] = -100;

		###update W

    id = which(FC < 0 & p.new<1e-3);

		dif = Z[id,]-rowMeans(Y[id,]);
		W.est = apply(dif,2,mean,trim=0.05);
		W.est = as.numeric(W.est/ mean( W.est ));
		inv.W.min = 1/min(W.est);
		W.max = max(W.est);
		W.all = rbind(W.all,W.est);

		t = t+1;
	}
	LR = cbind(lrt=lrt,p.value=p.new,logFC=FC);
	return(list(W=W.all[2,],LR=LR));
}
