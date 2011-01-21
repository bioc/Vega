\name{Vega-package}
\alias{Vega-package}
\alias{Vega}
\docType{package}
\title{An R package for copy number data segmentation}
\description{Vega (Variational Estimator for Genomic Aberrations) is an algorithm that adapts a very popular variational model (Mumford and Shah) used in image segmentation so that chromosomal aberrant regions can be efficiently detected.
}
\details{
\tabular{ll}{
Package: \tab Vega\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
biocViews: \tab aCGH, DNACopyNumber\cr
}
}

\references{ Morganella S. et al. (2010). VEGA: Variational segmentation for copy number detection, Bioinformatics.}
\author{Sandro Morganella, Luigi Cerulo, Giuseppe Viglietto, Michele Ceccarelli
Maintainer: Sandro Morganella <morganellaalx@gmail.com>}


\examples{

# Load tha data for Granta-519 Cell Line
data(G519);

# Start Vega segmentation on all chromosomes of G519 data and save the results in G519_segmentation.txt tab delimited file
seg <- vega(CNVdata=G519, chromosomes=c(1:22, "X", "Y"), out_file_name="G519_segmentation.txt", beta=0.5, min_region_size=2);
}
