f1 <-
function(mu,I,phij,W.est,yj,zj){  #mu = c(alpha,(1...,I-1),muj,deltaj)
	###par:
	alpha = c(mu[1:(I-1)],-sum(mu[1:(I-1)]));
	e.alpha = exp(alpha);

	mj= mu[I];
	deltaj = mu[I+1];
	######

	mu1.phi = 1/(1 + mj*e.alpha * phij); ##vector of length I.
	mu2.phi = 1/(1 + (mj + W.est*deltaj) * e.alpha * phij);
		
	L1 = 1/phij * log(mu1.phi) + yj * log(1-mu1.phi);
	L2 = 1/phij * log(mu2.phi) + zj * log(1-mu2.phi);
	return(-sum(L1+L2));
}
