\name{contamDE}
\alias{contamDE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   DE analysis using contaminated tumor samples
}
\description{
   Conduct differentially expression analysis between contaminated tumor samples and normal samples. This function provides estimates of relative proportions of pure tumor cells in the contaminated samples, and pesudo likelihood ratio statistics, p-values, and log2-fold change for detecting differentially expressed genes between tumor samples and normal samples. Here the reported proportions are scaled to have average value 1, so that some of them are greater than 1 and some of them are smaller than 1. As a disadvantage, the absulute proportions cannot be estimated without using extra information, refer to Shen et al. (2016) for details. A computationally very fast algorithm implemented in the R package `nloptr' is used to obtain the fold changes, which may report some warning messages that do not affect the final results so they can be safely ignored.
}
\usage{
contamDE(data,R,n=NA,match=TRUE)
}
\arguments{
  \item{data}{
     This is a G x N read count matrix, where G is the number of genes and N is the total sample size. The (g,i)th entry is the read count of the gth gene of the ith sample.}

  \item{R}{The number of conditions.}
   \item{n}{A list of sample sizes for all conditions(valid if match=FALSE).}
  \item{match}{
    TRUE if the tumor samples are matched with normal samples and FALSE otherwise.
    If matched=TRUE, then in `data', the ith, (i+N/R)th,...,(i+(R-1)N/R)th samples are matched with each other (i=1,...,N/R).
  }
}
\value{
    \item{W}{Scaled proportions (average value = 1) of pure tumor cells in  contaminated tumor samples.}
    \item{LR}{A J x (R+1) matrix. Column 1: pseudo likelihood ratio statistics for DE analysis; column 2: p-value for differential expression analysis;  column 3~(1+R): log2 fold changes (cancer vs. normal); }
}

\examples{
   data("prostate");

   ## Not run
   # d = contamDE(prostate[,-1],R=2,match=TRUE);

   data("drosophila");

   ## Not run
   # d = contamDE(drosophila[,-1],R=2,n=list(1:4,5:7),match=FALSE);
}

\references{
SHEN Qi, HU Jiyuan, JIANG Ning, HU Xiaohua, LUO Zewei, and ZHANG Hong (2016) contamDE: Differential expression analysis of RNA-seq data for contaminated tumor samples. Bioinformatics 32(5): 705-712.
}


\keyword{ differential expression }
\keyword{ contamination}
