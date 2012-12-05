\name{seqtree}
\alias{seqtree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Tree structured analysis of a state sequence object.
}
\description{
	Facility for growing a regression tree for a state sequence object.
}
\usage{
seqtree(formula, data = NULL, weighted = TRUE, minSize = 0.05,
     maxdepth = 5, R = 1000, pval = 0.01,
	 weight.permutation = "replicate",
	 seqdist_arg = list(method = "LCS", norm = TRUE),
     diss = NULL, squared = FALSE, first = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a formula where the left hand side is a state sequence object (see \code{\link{seqdef}})  and the right hand specifies the candidate variables for partitioning the set of sequences.}
  \item{weighted}{Logical. If \code{TRUE}, use the weights of the state sequence object.}
  \item{data}{a data frame where variables in the formula will be searched}
  \item{minSize}{minimum number of cases in a node, in percentage if less than 1.}
  \item{maxdepth}{maximum depth of the tree.}
  \item{R}{Number of permutations used to assess the significance of the split.}
  \item{pval}{Maximum p-value, in percent.}
  \item{weight.permutation}{Weights permutation method: "diss" (attach weights to the dissimilarity matrix), "replicate" (replicate case according to the \code{weights} arguments), "rounded-replicate" (replicate case according to the rounded \code{weights} arguments), "random-sampling" (random assignment of covariate profiles to the objects using distributions defined by the weights.)}
  \item{seqdist_arg}{list of arguments directly passed to \code{\link{seqdist}}, only used if \code{diss=NULL}}
  \item{diss}{An optional dissimilarity matrix. If not provided, a dissimilarity matrix is computed using \code{\link{seqdist}} and \code{seqdist_arg}}
  \item{squared}{Logical. If \code{TRUE}, the dissimilarity matrix is squared}
  \item{first}{Character. An optional variable name to force the first split.}
}
\details{
	The function provides a simplified interface for applying \code{\link{disstree}} on state sequence objects.
	
	The \code{seqtree} objects can be "plotted" with \code{\link{seqtreedisplay}}. A print method is also available which prints the medoid sequence for each terminal node.
}
\value{
	a seqtree object.
}

\references{
Studer, M., G. Ritschard, A. Gabadinho and N. S. Müller (2011). Discrepancy analysis of state sequences, \emph{Sociological Methods and Research}, Vol. 40(3), 471-510.
}

\seealso{
	\code{\link{seqtreedisplay}}, \code{\link{disstree}}
}


\examples{
data(mvad)

## Defining a state sequence object
mvad.seq <- seqdef(mvad[, 17:86])

## Building a seqtree using hamming distance
## You should use much higher value for the R parameter,
## thus the results may be unstable.
## R should be at least 1000
## (R is keeped low to avoid too much computation in examples).
seqt <- seqtree(mvad.seq~ male + Grammar + funemp + gcse5eq + fmpr + livboth,
    data=mvad, R = 10, seqdist_arg=list(method="HAM", norm=TRUE))
print(seqt)

## Building a seqtree using a specified distance matrix
mvad.dhd <- seqdist(mvad.seq, method="DHD")
seqt <- seqtree(mvad.seq~ male + Grammar + funemp + gcse5eq + fmpr + livboth,
    data=mvad, R = 10, diss=mvad.dhd)
print(seqt)


### Will only work if GraphViz is properly installed
\dontrun{
seqtreedisplay(seqt, type="d", border=NA)
seqtreedisplay(seqt, type="I", sortv=cmdscale(mvad.dhd, k=1))
}

}

\author{Matthias Studer (with Gilbert Ritschard for the help page)}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Dissimilarity-based analysis}
