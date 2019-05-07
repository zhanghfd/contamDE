impute.par <-
function(x){
	if(sum(is.na(x))){
	ind = which(is.na(x));
	x[ind] = 1e-3;
	}
	return(x);
}
