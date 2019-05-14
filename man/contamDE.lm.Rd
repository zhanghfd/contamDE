\name{contamDE.lm}
\alias{contamDE.lm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   DE analysis using RNA-seq data from paired tumor samples
}
\description{
   This function implements the DE analysis method contamDE-lm developed in Ji et al. (2019). The method contamDE-lm uses RNA-seq data from paired tumor samples, which can account for multiple tumor sutypes and further covariates if available. contamDE-lm is computationally much faster than contamDE, especailly when the pair number is large.
}
\usage{
contamDE.lm(counts, subtype=NULL, covariate=NULL, is.contaminated=TRUE, robust = TRUE)
}
\arguments{
  \item{counts}{
    This is a G x (2N) read count matrix, where G is the number of genes and N is the number of patients, each patient having both normal and tumor samples. The first N samples are for normals and the last N samples are for tumors. The ith and the (N+i)th samples are from the same patient.
  }
  \item{subtype}{
    This argument specifies the tumor subtype information for each patient. The default value of this argument is NULL (i.e., the subtype information is ignored).
  }
  \item{covariate}{
    Further covariate information can be specified. This argument should be a vector of length N or a N by q matrix (q is the number of covariates) it is not NULL (default value).
  }
  \item{is.contaminated}{
    If is.contaminated = TRUE (default value), then scaled proportion of pure tumor cells of each pacient will be estimated. If is.contaminzed = FALSE, then contamination information will be ignored.
  }
  \item{robust}{
  If robust = TRUE (default value),then a robust strategy will be applied on residual variance estimation which may work on data with outliers. If robust = FALSE, then a common empirical bayes method will be used.
  }
}
\value{
  \item{p.contamDE.lm}{
    T-test p-values for DE analysis.
  }
  \item{design}{
    Design matrix.
  }
  \item{log2FC}{
    Estimated log2-fold changes.
  }
  \item{log2FC.cov}{
    Estimated ariance-covariance of the estimated log2-fold changes. Each row is a vectorized variance-covariance matrix for a gene.
  }
  \item{proportion}{
    Estimated pure proportions of contaminated tumor samples, which are scaled to have a maximal value 1.
  }
  \item{df}{
    Degree of freedom for t-tests.
  }
  \item{counts}{
    Original counts for analysis
  }
}


\examples{
\dontrun{
 library("contamDE");

 data("simcounts");

 subtype <- c(rep("Type1",4),rep("Type2",2),rep("Type3",2));

 d <- contamDE.lm(simcounts,subtype)
}
}

\references{
JI Yifan, YU Chang, and ZHANG Hong (2019) contamDE-lm: Linear model based differential gene expression analysis using next-generation RNA-seq data from contaminated tumor samples.
}


\keyword{ differential expression analysis }
\keyword{ t-test }

