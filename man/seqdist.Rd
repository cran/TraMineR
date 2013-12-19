\encoding{latin1}
\name{seqdist}
\alias{seqdist}
\title{Distances (dissimilarities) between sequences}

\description{
 Computes pairwise dissimilarities between sequences or
dissimilarities with a reference sequence. Several dissimilarities measures or metrics are available: optimal matching (OM), distance based on the longest common prefix (LCP), on the longest common suffix (RLCP), on the longest common subsequence (LCS), the Hamming distance (HAM) and the Dynamic Hamming Distance (DHD).
 }

\usage{
seqdist(seqdata, method, refseq=NULL, norm=FALSE,
     indel=1, sm=NA, with.missing=FALSE, full.matrix=TRUE)
}

\arguments{
  \item{seqdata}{a state sequence object defined with the \code{\link{seqdef}} function.}

  \item{method}{a character string indicating the metric to be used. One of \code{"OM"} (Optimal Matching), \code{"LCP"} (Longest Common Prefix), \code{"RLCP"} (reversed LCP, i.e. Longest Common Suffix), \code{"LCS"} (Longest Common Subsequence), \code{"HAM"} (Hamming distance), \code{"DHD"} (Dynamic Hamming distance).}

  \item{refseq}{Optional baseline sequence to compute the distances from. Can be the index of a sequence in the state sequence object, 0 for the most frequent sequence, or an external sequence passed as a sequence object with 1 row.}

  \item{norm}{if \code{TRUE}, the computed OM, LCP, RLCP or LCS distances are normalized to account for differences in sequence lengths, and the normalization method is automatically selected.  Default is \code{FALSE}. Can also be one of \code{"none"}, \code{"maxlength"}, \code{"gmean"}, \code{"maxdist"}, \code{"YujianBo"}. See details.}

  \item{indel}{the insertion/deletion cost (OM method). Default is 1. Ignored with non OM metrics.}

  \item{sm}{substitution-cost matrix (OM, HAM and DHD method). Can also be one of the \code{\link{seqsubm}} build methods \code{"TRATE"} or \code{"CONSTANT"}. Default is \code{NA}. Ignored with LCP, RLCP and LCS metrics.}

  \item{with.missing}{must be set to \code{TRUE} when sequences contain non deleted gaps (missing values). See details.}

  \item{full.matrix}{If \code{TRUE} (default), the full distance matrix is returned. This is for compatibility with  earlier versions of the \code{seqdist} function. If \code{FALSE}, an object of class \code{\link{dist}} is returned, that is, a vector containing only values from the upper triangle of the distance matrix. Since the distance matrix is symmetrical, no information is lost with this representation while size is divided by 2. Objects of class \code{dist} can be passed directly as arguments to most clustering functions. Ignored when \code{refseq} is set.}
}

\details{
The \code{seqdist} function returns a matrix of distances between sequences or a vector of distances to a reference sequence. The available metrics (see 'method' option) are optimal matching (\code{"OM"}), longest common prefix (\code{"LCP"}), longest common suffix (\code{"RLCP"}), longest common subsequence (\code{"LCS"}), Hamming distance (\code{"HAM"}) and Dynamic Hamming Distance (\code{"DHD"}). The Hamming distance is OM without indels and the
Dynamic Hamming Distance is HAM with specific substitution costs at each position as proposed by \cite{Lesnard (2006)}. Note that HAM and DHD apply only to sequences of equal length.

For OM, HAM and DHD, a user specified substitution cost matrix can be provided with the \code{sm} argument. For DHD, this should be a series of matrices grouped in a 3-dimensional matrix with the third index referring to the position in the sequence. When \code{sm} is not specified, a constant substitution cost of 1 is used with HAM, and \cite{Lesnard (2006)}'s proposal for DHD.

Distances can optionally be normalized by means of the \code{norm} argument. If set to \code{TRUE}, Elzinga's normalization (similarity divided by geometrical mean of the two sequence lengths) is applied to LCP, RLCP and LCS distances, while Abbott's normalization (distance divided by length of the longer sequence) is used for OM, HAM and DHD. Elzinga's method can be forced with \code{"gmean"} and Abbott's rule with \code{"maxlength"}. With \code{"maxdist"} the distance is normalized by its maximal possible value. For more details, see \cite{Elzinga (2008)} and \cite{Gabadinho et al. (2009)}.

When sequences contain gaps and the \code{gaps=NA} option was passed to \code{\link{seqdef}}, i.e. when there are non deleted missing values, the \code{with.missing} argument should be set to \code{TRUE}. If left to \code{FALSE} the function stops when it encounters a gap. This is to make the user aware that there are gaps in his sequences. If the OM method is selected, \code{seqdist} expects a substitution cost matrix with a row and a column entry for the missing state (symbol defined with the \code{nr} option of \code{\link{seqdef}}). This will be the case for substitution cost matrices returned by \code{\link{seqsubm}}. More details on how to compute distances with sequences containing gaps are given in \cite{Gabadinho et al. (2009)}.
}

\value{
 When \code{refseq} is specified, a vector with distances between the sequences in the data sequence object and the reference sequence is returned. When \code{refseq} is \code{NULL} (default), the whole matrix of pairwise distances between sequences is returned.
 }

\seealso{
 \code{\link{seqsubm}}, \code{\link{seqdef}}, and for multichannel distances \code{\link{seqdistmc}}.
}

\references{
Elzinga, Cees H. (2008). Sequence analysis: Metric representations of categorical time
series. \emph{Technical Report}, Department of Social Science Research Methods, Vrije Universiteit, Amsterdam.

Gabadinho, A., G. Ritschard, N. S. Müller and M. Studer (2011). Analyzing and Visualizing State Sequences in R with TraMineR. \emph{Journal of Statistical Software} \bold{40}(4), 1-37.

Gabadinho, A., G. Ritschard, M. Studer and N. S. Müller (2009). Mining Sequence Data in \code{R} with the \code{TraMineR} package: A user's guide. Department of Econometrics and Laboratory of Demography, University of Geneva

Lesnard, L. (2006) Optimal Matching and Social Sciences. \emph{Série des Documents de Travail du CREST},  Institut National de la Statistique et des Etudes Economiques, Paris.

}

\examples{
## optimal matching distances with substitution cost matrix
## derived from transition rates
data(biofam)
biofam.seq <- seqdef(biofam, 10:25)
costs <- seqsubm(biofam.seq, method="TRATE")
biofam.om <- seqdist(biofam.seq, method="OM", indel=3, sm=costs)

## normalized LCP distances
biofam.lcp <- seqdist(biofam.seq, method="LCP", norm=TRUE)

## normalized LCS distances to the most frequent sequence in the data set
biofam.lcs <- seqdist(biofam.seq, method="LCS", refseq=0, norm=TRUE)

## histogram of the normalized LCS distances
hist(biofam.lcs)

## =====================
## Example with missings
## =====================
data(ex1)
ex1.seq <- seqdef(ex1,1:13)

subm <- seqsubm(ex1.seq, method="TRATE", with.missing=TRUE)
ex1.om <- seqdist(ex1.seq, method="OM", sm=subm, with.missing=TRUE)
}

\author{Matthias Studer, Alexis Gabadinho and Nicolas S. Müller (first version) (with Gilbert Ritschard for the help page)}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Dissimilarity measures}
