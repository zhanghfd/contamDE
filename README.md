# contamDE
contamDE: Differential expression analysis of RNA-seq data for contaminated tumor samples

The R package ‘contamDE’ conducts differential expression (DE) analysis using high throughput next-generation RNA-seq read count data generated from contaminated tumor samples that are either matched or unmatched with normal samples, which estimates the proportion of pure tumor cells in each contaminated tumor sample, and provides tumor vs. normal log2-fold change, likelihood ratio test statistic and p-value of DE analysis for each gene.

The manual file is "contamDE-manual.pdf". 

Installation of contamDE in R:

> library(‘devtools’);

> install_github(‘zhanghfd/contamDE’);
