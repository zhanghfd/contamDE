\name{drosophila}
\alias{drosophila}
\docType{data}
\title{
Drosophila melanogaster dataset
}
\description{
The dataset consists of the RNA-seq read counts for 7 samples of Drosophila melanogaster S2 cells, of which 4 samples were untreated while 3 samples were treated with siRNA targeting the splicing factor pasilla (CG1844).
}
\usage{data("drosophila")}
\format{
  A data frame with 7196 observations (genes) on the following 8 variables.
  \describe{
    \item{\code{nameOfGene}}{Gene name}
    \item{\code{CT.PA.1}}{Read count from Untreated-3.}
    \item{\code{CT.PA.2}}{Read count from Untreated-4.}
    \item{\code{CT.SI.5}}{Read count from Untreated-1.}
    \item{\code{CT.SI.7}}{Read count from Untreated-6.}
    \item{\code{KD.PA.3}}{Read count from "CG8144_RNAi-3".}
    \item{\code{KD.SI.6}}{Read count from "CG8144_RNAi-1".}
    \item{\code{KD.PA.4}}{Read count from "CG8144_RNAi-4".}
  }
}
\source{Brooks, A.N., Yang, L., Duff, M.O., Hansen, K.D., Park, J.W., Dudoit, S., Brenner, S.E., and Graveley, B.R. (2011). Conservation of an RNA regulatory map between Drosophila and mammals. \emph{Genome Research}, 21(2), 193--202.}
\keyword{Drosophila melanogaster}
