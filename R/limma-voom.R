
limma_voom <- function(counts){

  d <- DGEList(counts)
  d <- calcNormFactors(d)
  size <- d$samples[,2]*d$sample[,3]
  size <- size/mean(size)

  N <- ncol(counts)/2; # number of patients
  pair <- as.factor(rep(1:N,2))

  condition <- as.factor(c(rep(0,N),rep(1,N)))

  design <- model.matrix(~0+condition+pair)
  d <- voom(d, design, normalize.method="quantile")
  d <- lmFit(d, design)
  contr <- 'condition1-condition0'
  contrs <- makeContrasts(contrasts= contr,levels=colnames(design))
  d <- contrasts.fit(d,contrs)
  d <- eBayes(d,trend=TRUE)
  p.limma <- d$p.value
  log2FC.limma <- d$coefficient
  res <- list(counts=counts,size=size,p.limma=p.limma,log2FC.limma=log2FC.limma)

  return(res)

}
