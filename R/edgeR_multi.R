

edgeR_multi = function(count,R,nncol=nncol,is.matched=T){

  #count = round(count);

  dm = dim(count);
  J=dm[1]


  FC =matrix(NA,nrow=J,ncol=(R-1))


  if (!is.matched){
    condition1=NULL
    for(i in 1:R){
      condition2=rep(i,nncol[i]);
      condition1=c(condition1,condition2)
    }
    condition1=as.factor(condition1)
      d = DGEList(counts=count, group=condition1);
      d = calcNormFactors(d);

      design = model.matrix(~condition1);
      d = estimateGLMCommonDisp(d, design);
      d = estimateGLMTagwiseDisp(d, design);

      fit = glmFit(d, design);
      result = glmLRT(fit, coef=2:R);

      p.value = result$table$PValue;
      dis = as.numeric(result$dispersion);
      for( i in 1:(R-1))
      {

        FC[,i]=result$table[[i]]
      }
     
  }
  else {
    
    I = dm[2]/R;
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
      
      dis = as.numeric(result$dispersion);
      for( i in 1:(R-1))
{
        FC[,i]=result$table[[i]]
      }
  
  }
#  write.table(p.value,file="p.edgeR.txt",col.names=FALSE,row.names=FALSE,quote=FALSE);
return(list(p.edgeR=p.value,FC=FC,dispersion=dis))  
}
