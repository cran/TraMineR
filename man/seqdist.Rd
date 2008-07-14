\name{seqdist}
\alias{seqdist}
\title{Compute distances between sequences}
\description{
Compute distances between sequences. Several metrics are available: optimal matching and other metrics proposed by Elzinga.
}
\usage{
seqdist(seqdata, method, refseq=NULL, norm=FALSE, indel=1, sm)
}
\arguments{
  \item{seqdata}{a sequence object (see \code{\link{seqdef}} function).}
  \item{method}{a character string indicating the metric to use for distance. One of "OM" (optimal matching),"LCP" (Longest Common Prefix), "LCS" (Longest Common Subsequence).}
  \item{refseq}{Optional reference sequence to compute the distances from. Can be the index of a sequence in the data set or 0 for the most frequent sequence in the data set. If refseq is specified, a vector with distances between the sequences in the data set and the reference sequence is returned. If refseq is not specified (default), the distance matrix containing the distance between all sequences in the data set is returned.}
  \item{norm}{if TRUE, OM, LCP and LCS distances are rescaled to be unit free, ie insensitive to sequences length. Default to FALSE.}
  \item{indel}{the insertion/delation cost if optimal matching ("OM") is choosed. Default to 1. Don't specify if other metric is used.}
  \item{sm}{substitution-cost matrix for the optimal matching method ("OM"). Default to NA. Don't specify if other method is used.}
}
\details{
	The seqdist function returns a matrix of distances between sequences or a vector of distances to a reference sequence. The available metrics (see 'method' option) are optimal matching ("OM"), longuest common prefix ("LCP") or longuest common subsequence ("LCS"). LCP and LCS distances can optionaly be normalized (see 'norm' option).
}
\value{a distance matrix or a vector containing distances to the specified reference sequence.}
\seealso{
 \code{\link{seqsubm}}.
}
\examples{
## optimal matching distances with substitution cost matrix 
## using transition rates
data(biofam)
biofam.seq <- seqdef(biofam,10:25,informat='STS')
costs <- seqsubm(biofam.seq,method="TRATE")
biofam.om <- seqdist(biofam.seq,method="OM",indel=3,sm=costs)

## normalized LCP distances
biofam.lcp <- seqdist(biofam.seq,method="LCP", norm=TRUE)

## normalized LCS distances to the most frequent sequence in the data set
biofam.lcs <- seqdist(biofam.seq,method="LCS", refseq=0, norm=TRUE)

## histogram of the normalized LCS distances
hist(biofam.lcs)
}
\keyword{misc}
