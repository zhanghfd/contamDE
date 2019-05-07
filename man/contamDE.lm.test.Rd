\name{contamDE.lm.test}
\alias{contamDE.lm.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   F-test for testing linear hypothesis
}
\description{
This function conducted linear linear hypothesis for the method contamDE-lm developed in Ji et al. (2019). The null hypothesis can be generally specified in a linear form H beta = b, where H is a r x (p+q) matrix of rank r, b is a r-vector, and beta is the regression parameters in the linear model (Ji et al., 2019).
}
\usage{
contamDE.lm.test(d,H,b=NULL)
}
\arguments{
  \item{d}{
    A list returned by the function 'contamDE.lm'. When calling contamDE.lm, subtype information should be provided.
  }
   \item{H}{
    a r x (p+q) design matrix of rank r.}
  \item{b}{
    a p vector with the default value NULL (=0).
  }
}
\value{
  \item{p.contamDE.lm.test}{
    F-test p-values for DE analysis.
  }
}

\examples{
  \dontrun{
    library("contamDE");

    data("simcounts");

    subtype <- c(rep("Type1",4),rep("Type2",2),rep("Type3",2));

    d <- contamDE.lm(simcounts,subtype=subtype);

    # overall DE analysis between the tumors and the normal
    H <- diag(0,3);
    d1 <- contamDE.lm.test(d,H);

    # "subtype 2 versus normal" DE analysis
    H <- matrix(c(0,1,0),1,3);
    d1 <- contamDE.lm.test(d,H);
  }
}

\references{
JI Yifan, YU Chang, and ZHANG Hong (2019) contamDE-lm: A linear model for differential gene expression analysis using next generation RNA-seq data from contaminated tumor samples.
}

\keyword{ linear hypothesis }
\keyword{ F test}
