g0 <-
function(mu,I,phij,yj,zj){
	yj1 = yj[1:(I-1)];
	yj2 = yj[I];
	zj1 = zj[1:(I-1)];
	zj2 = zj[I];
	###par:
	alpha = c(mu[1:(I-1)],-sum(mu[1:(I-1)]));
	e.alpha = exp(alpha);
	mj= mu[I];
	#######

	mu.phi = 1 + mj * phij * e.alpha;  #the inverse of the mu.phi function in the f0 func above
	mu.phi1 = mu.phi[1:(I-1)];mu.phi2 = mu.phi[I];
	
	mu.alpha = mj*e.alpha;
	mu.alpha1 = mu.alpha[1:(I-1)];mu.alpha2 = mu.alpha[I];
	####
	partial_alpha = -2*mu.alpha1/mu.phi1+2 * mu.alpha2 / mu.phi2 + (yj1+zj1)/mu.phi1 - (yj2+zj2)/mu.phi2;  #length of I-1
	
	partial_mj = sum(-2*e.alpha/mu.phi+(yj+zj)/(mu.phi*mj));
	
	return(-c(partial_alpha,partial_mj));
}
