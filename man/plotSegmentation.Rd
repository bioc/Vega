\name{plotSegmentation}
\alias{plotSegmentation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot observations and the respective segmentation.}
\description{
 This function allows to plot the observed data superimposing the respective segmentation. By the parameter 'opt' he user can plot the LRR mean values of each segment or the computed aberration kind. In plot window the gain and the loss are identified by a line having value of 1 and -1 respectively.
}
\usage{
plotSegmentation(CNVdata, segmentation, chromosomes, opt = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CNVdata}{The data matrix.}
  \item{segmentation}{The computed segmentation.}
  \item{chromosomes}{The chromosomes that have to be plotted.}
  \item{opt}{If opt=0 (default) then the LRR mean values are plotted, if opt=1 the computed labels are plotted.}
}
\references{ Morganella S. et al. (2010). VEGA: Variational segmentation for copy number detection, Bioinformatics.}
\author{Sandro Morganella, Luigi Cerulo, Giuseppe Viglietto, Michele Ceccarelli
Maintainer: Sandro Morganella <morganellaalx@gmail.com>}

\note{If the argument opt=1 then gains and losses are identified by 1 and -1 respectively.
}
\examples{
	
	# Import the data
	data(G519)
	
	# Compute the segmentation for all chromosomes
	seg <- vega(G519, c(1:22, "X", "Y"))
	
	# Plot the results for all chromosomes in terms of mean of LRRs 
	plotSegmentation(G519, seg, c(1:22, "X", "Y"), opt=0)

	# Plot the results for all chromosomes in terms of aberration kinds
	plotSegmentation(G519, seg, c(1:22, "X", "Y"), opt=1)
	
}
