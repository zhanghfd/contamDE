\name{lung}
\alias{lung}
\docType{data}
\title{Lung cancer data}
\description{
The data are from a study of the lung cancer. Six patients provided tissue samples and normal samples besides the lung tissues. The read counts were summarized by RefSeq transcript, and only those transcripts with at least 50 aligned reads for at
least one tissue in each condition were provided in the table. RefSeq identifiers were mapped to the latest official gene symbols by following the user guide of the Bioconductor package `edgeR' using the
Bioconductor annotation package `org.Hs.eg.db' (version 2.7.1). Those RefSeq identifiers not in the database were discarded, and each gene was represented by the RefSeq transcript with the
greatest number of exons and the other transcripts were removed. Altogether 11,597 transcripts (genes) were kept.}
\usage{data(lung)}
\format{
  A data frame with 11,597 observations on the following 13 variables.
  \describe{
    \item{\code{nameOfGene}}{Gene name}
    \item{\code{N4}}{Read count for normal sample of patient 4.}
    \item{\code{T4}}{Read count for normal sample of patient 4.}
    \item{\code{N12}}{Read count for normal sample of patient 12.}
    \item{\code{T12}}{Read count for tumor sample of patient 12.}
    \item{\code{N13}}{Read count for normal sample of patient 13.}
    \item{\code{T13}}{Read count for tumor sample of patient 13.}
    \item{\code{N14}}{Read count for normal sample of patient 14.}
    \item{\code{T14}}{Read count for tumor sample of patient 14.}
    \item{\code{N15}}{Read count for normal sample of patient 15.}
    \item{\code{T15}}{Read count for tumor sample of patient 15.}
    \item{\code{N16}}{Read count for normal sample of patient 16.}
    \item{\code{T16}}{Read count for tumor sample of patient 16.}
  }
}

\keyword{dataset}
