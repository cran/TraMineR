\encoding{latin1}
\name{seqdist}
\alias{seqdist}
\title{Distances between sequences}
\description{
Compute pairwise distances between sequences or distances to a reference sequence. Several metrics are available: optimal matching (OM) and other metrics such as the longest common prefix (LCP), the longest common suffix (RLCP) and the longest common subsequence (LCS).
}
\usage{
seqdist(seqdata, method, refseq, norm=FALSE, 
	indel=1, sm, with.miss = FALSE, full.matrix = TRUE)
}
\arguments{
\item{seqdata}{a state sequence object defined with the \code{\link{seqdef}} function.}

  \item{method}{a character string indicating the metric to be used. One of "OM" (Optimal Matching), "LCP" (Longest Common Prefix), "RLCP" (reversed LCP, i.e. Longest Common Suffix), "LCS" (Longest Common Subsequence).}
  \item{refseq}{Optional reference sequence to compute the distances from. Can be the index of a sequence in the state sequence object or 0 for the most frequent sequence, or an external sequence passed as a sequence object with 1 row. If refseq is specified, a vector with distances between the sequences in sequence object and the reference sequence is returned. If refseq is not specified (default), the whole matrix of pairwise distances sequences in the data sequence object is returned.}
  \item{norm}{if TRUE, the computed OM, LCP, RLCP or LCS distances are normalized to account for differences in sequence lengths. Default is FALSE. See details}
  \item{indel}{the insertion/deletion cost ("OM" method). Default is 1. Ignored with non OM metrics.}
  \item{sm}{substitution-cost matrix ("OM" method). Default is NA. Ignored with non OM metrics.}
  \item{with.miss}{Must be set to TRUE when sequences contain non deleted gaps (missing values). See details.} 
  \item{full.matrix}{If TRUE (default), the full distance matrix is returned. This is for compatibility with the previous version of the \code{seqdist} function. If false, an object of class \code{\link{dist}} is returned, that is, a vector containing only the half distance matrix. Since the distance matrix is symmetrical, no information is lost with this representation but the size is divided by 2. Objects of class dist can be passed directly as arguments to many clustering functions.}

}
\details{
The seqdist function returns a matrix of distances between sequences or a vector of distances to a reference sequence. The available metrics (see 'method' option) are optimal matching ("OM"), longest common prefix ("LCP"), longest common suffix ("RLCP") and longest common subsequence ("LCS"). 

Distances can optionally be normalized by means of the \code{norm} option). If set to TRUE, Elzinga's normalization (similarity divided by geometrical mean of the two sequence lengths) is applied to LCP, RLCP and LCS distances, while Abbott's normalization (distance divided by length of the longer sequence) is used for OM. For more details, see \cite{Elzinga (2008)} and \cite{Gabadinho et al. (2009)}.


When sequences contain gaps and the \code{gaps=NA} option was passed to \code{\link{seqdef}}, i.e. when there are non deleted missing values, the \code{with.miss} argument should be set to TRUE. If left to FALSE the function stops when it encounters a gap. This is to make the user aware that there are gaps in his sequences. If "OM" method is selected, \code{seqdist} expects a substitution cost matrix with a row and a column entry for the missing state (symbol defined with the \code{nr} option of \code{\link{seqdef}}). This will be the case for substitution cost matrices returned by \code{\link{seqsubm}}. More details on how to compute distances with sequences containing gaps are given in \cite{Gabadinho et al. (2009)}.
}
\value{a distance matrix or a vector containing distances to the specified reference sequence.}
\seealso{
 \code{\link{seqsubm}}, \code{\link{seqdef}}.
}

\references{
Elzinga, Cees H. (2008). Sequence analysis: Metric representations of categorical time
series. \emph{Sociological Methods and Research}, forthcoming.

Gabadinho, A., G. Ritschard, M. Studer and N. S. Müller (2009). Mining Sequence Data in \code{R} with \code{TraMineR}: A user's guide for version 1.1. \emph{Department of Econometrics and Laboratory of Demography, University of Geneva}.
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
