\name{vega}
\alias{vega}

\title{This function computes, saves and returns the copy number segmentation on the aCGH data passed as argument.}

\description{
  This function computes the segmentation of the copy number data passed as argument. 'vega' function returns the computed segmentation and save it as a tab delimited file. Users need just to use this function.
}
\usage{
vega(CNVdata, chromosomes, out_file_name = "", beta = 0.5, min_region_size = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CNVdata}{This argument is matrix containing the data that have to be analyzed. This matrix must have 4 columns:\cr 
		- the first row indicates the chromosome;\cr
		- the second row indicates the start bp of the probe;\cr
		- the third row indicates the end bp of the probe;\cr
		- the fourth row reports the measured Log R Ratio;\cr

}

  \item{chromosomes}{This is a vector containing the chromosome sthat have to be analyzed. By using c(1:22, "X", "Y") the whole genome will be segmented.
  }
  \item{out_file_name}{This is the file name used to save the computed segmentation.}
  \item{beta}{This argument is used for the stop condition definition.}
  \item{min_region_size}{This argument specifies the minimum size for the segmented regions.}
}

\value{
  \item{segmentation}{This is a data frame containing the computed segmentation.}
}

\references{ Morganella S. et al. (2010). VEGA: Variational segmentation for copy number detection, Bioinformatics.}
\author{Sandro Morganella, Luigi Cerulo, Giuseppe Viglietto, Michele Ceccarelli
Maintainer: Sandro Morganella <morganellaalx@gmail.com>}

\examples{

# Load tha data for Granta-519 Cell Line
data(G519);

# Start Vega segmentation on all chromosomes of G519 data 
v <- vega(CNVdata=G519, chromosomes=c(1:22, "X", "Y"), beta=0.5, min_region_size=2);

# Start Vega segmentation on chromosomes 8 of G519 data and save the results in G519_segmentation_chr8.txt tab delimited file
v <- vega(CNVdata=G519, chromosomes=c(8), out_file_name="G519_segmentation_chr8.txt", beta=0.5, min_region_size=2);

}

