f0 <-
function(mu,I,phij,yj,zj){  ##mu = c(alpha,(1,...,I-1),muj)
	###par:
	alpha = c(mu[1:(I-1)],-sum(mu[1:(I-1)]));
	e.alpha = exp(alpha);
	mj= mu[I];
	#######

	mu.phi = 1/(1 + e.alpha * mj * phij);

	L1 = 1/phij * log(mu.phi) + yj * log(1-mu.phi);
	L2 = 1/phij * log(mu.phi) + zj * log(1-mu.phi);

	return(-sum(L1+L2));
}
