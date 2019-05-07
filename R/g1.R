g1 <-
function(mu,I,phij,W.est,yj,zj){
	yj1 = yj[1:(I-1)];
	yj2 = yj[I];
	zj1 = zj[1:(I-1)];
	zj2 = zj[I];

	###par:
	e.alpha = exp(c(mu[1:(I-1)],-sum(mu[1:(I-1)])));
	mj= mu[I];
	deltaj = mu[I+1];
	#######

	mu1.phi = 1 + mj * phij*e.alpha; ##the inverse of the mu1.phi above in the f1 function.
	mu2.phi = 1 + (mj+W.est* deltaj) * phij*e.alpha;
	mu1.phi1 = mu1.phi[1:(I-1)]; mu1.phi2 = mu1.phi[I];
	mu2.phi1 = mu2.phi[1:(I-1)]; mu2.phi2 = mu2.phi[I];
	
	mu1.alpha = mj*e.alpha;
	mu2.alpha = (mj+W.est*deltaj)*e.alpha;
	mu1.alpha1 = mu1.alpha[1:(I-1)];mu1.alpha2 = mu1.alpha[I];
	mu2.alpha1 = mu2.alpha[1:(I-1)];mu2.alpha2 = mu2.alpha[I];
	
	e.alpha1 = e.alpha[1:(I-1)]; e.alpha2 = e.alpha[I];
	
	L1 = -mu1.alpha1/mu1.phi1 + mu1.alpha2/mu1.phi2 + yj1/mu1.phi1 - yj2/mu1.phi2;
	L2 = -mu2.alpha1/mu2.phi1 + mu2.alpha2/mu2.phi2 + zj1/mu2.phi1 - zj2/mu2.phi2;
	
	partial_alpha = L1+L2;  #vector of length I-1
	
	L3 = -e.alpha/mu1.phi + phij*yj*e.alpha/((mu1.phi-1)*mu1.phi);
	L4 = -e.alpha/mu2.phi + phij*zj*e.alpha/((mu2.phi-1)*mu2.phi);
	partial_mj = sum(L3+L4);
	
	partial_deltaj =  sum(-W.est*e.alpha/mu2.phi+phij*W.est*zj*e.alpha/((mu2.phi-1)*mu2.phi));
	return (-c(partial_alpha,partial_mj,partial_deltaj));
}
