edgeR = function(count,In,is.matched=FALSE){
  
  count = round(count);
  if(is.matched){
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
    p.value = result$table$PValue;
		FC = result$table$logFC;
    dis = as.numeric(result$dispersion);
  }else{
    dm = dim(count);
    R = 2;
    I = dm[2];
    Ic=I-In;
    dge = DGEList(counts=count);
    dge = calcNormFactors(dge);
  
    Tissue = c(rep(1,In),rep(2,Ic));
    design = model.matrix(~Tissue); #unmatched design
  
    dge = estimateCommonDisp(dge,design); #estimate dispersion
    dge = estimateGLMTagwiseDisp(dge,design);
  
    fit = glmFit(dge,design);
    fit = glmLRT(fit,2);
    dis = fit$dispersion;
    ###
    p.value = fit$table$PValue;  ##report edgeR's p-value.
    FC = fit$table$logFC;  #report edgeR's Fold Change.
  }

  return(list(p.value=p.value,dispersion=dis,FC = FC));
}

