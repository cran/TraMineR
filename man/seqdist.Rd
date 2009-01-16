\encoding{latin1}
\name{seqdist}
\alias{seqdist}
\title{Distances between sequences}
\description{
Compute distances between sequences. Several metrics are available: optimal matching and other metrics proposed by \cite{Elzinga (2008)}.
}
\usage{
seqdist(seqdata, method, refseq=NULL, norm=FALSE, 
	indel=1, sm, with.miss = FALSE, full.matrix = TRUE)
}
\arguments{
\item{seqdata}{a sequence object as defined by the the \code{\link{seqdef}} function.}

  \item{method}{a character string indicating the metric to use for distance. One of "OM" (optimal matching),"LCP" (Longest Common Prefix), "LCS" (Longest Common Subsequence).}
  \item{refseq}{Optional reference sequence to compute the distances from. Can be the index of a sequence in the data set or 0 for the most frequent sequence in the data set. If refseq is specified, a vector with distances between the sequences in the data set and the reference sequence is returned. If refseq is not specified (default), the distance matrix containing the distance between all sequences in the data set is returned.}
  \item{norm}{if TRUE, OM, LCP and LCS distances are rescaled to be unit free, ie insensitive to sequences length. Default to FALSE.}
  \item{indel}{the insertion/delation cost if optimal matching ("OM") is choosed. Default to 1. Don't specify if other metric is used.}
  \item{sm}{substitution-cost matrix for the optimal matching method ("OM"). Default to NA. Don't specify if other method is used.}
  \item{with.miss}{If sequences contain gaps (missing values) (see \code{\link{seqdef}} for the available options for handling missing values), this option must be set to TRUE to compute distances, otherwise the function will stop. If optimal matching method is used, the substitution cost matrix must contain one entry for the missing state. See \cite{Gabadinho et al. (2008)} for more details on how to compute distances with sequences containing gaps.} 

  \item{full.matrix}{If TRUE (default), the full distance matrix is returned. This is for compatibility with the previous version of the seqdist function. If false, an object of class \code{\link{dist}} is returned, that is, a vector containing only the half distance matrix. Since the distance matrix is symetrical, no information is lost with this representation but the size is divided by 2. Objects of class dist can be passed directly as arguments to many clustering functions.}

}
\details{
The seqdist function returns a matrix of distances between sequences or a vector of distances to a reference sequence. The available metrics (see 'method' option) are optimal matching ("OM"), longuest common prefix ("LCP") or longuest common subsequence ("LCS"). Distances can optionaly be normalized (see 'norm' option). For more details, see \cite{Elzinga (2008)} and \cite{Gabadinho (2008)}.
}
\value{a distance matrix or a vector containing distances to the specified reference sequence.}
\seealso{
 \code{\link{seqsubm}}.
}

\references{
Elzinga, Cees H. (2008). Sequence analysis: Metric representations of categorical time
series. \emph{Sociological Methods and Research}, forthcoming.

Gabadinho, A., G. Ritschard, M. Studer and N. S. Müller (2008). Mining Sequence Data in \code{R} with \code{TraMineR}: A user's guide. \emph{Department of Econometrics and Laboratory of Demography, University of Geneva}.
}

\examples{
## optimal matching distances with substitution cost matrix 
## using transition rates
data(biofam)
biofam.seq <- seqdef(biofam, 10:25)
costs <- seqsubm(biofam.seq, method="TRATE")
biofam.om <- seqdist(biofam.seq, method="OM", indel=3, sm=costs)

## normalized LCP distances
biofam.lcp <- seqdist(biofam.seq, method="LCP", norm=TRUE)

## normalized LCS distances to the most frequent sequence in the data set
biofam.lcs <- seqdist(biofam.seq,method="LCS", refseq=0, norm=TRUE)

## histogram of the normalized LCS distances
hist(biofam.lcs)
}
\keyword{misc}
