\name{seqdistmc}
\alias{seqdistmc}
\title{Multichannel distances between sequences}
\description{Compute multichannel pairwise optimal matching (OM) distances between sequences by deriving the substitution costs from the costs of the single channels.
Works with OM and its following variants:
distance based on longest common subsequence (LCS),
Hamming distance (HAM), and Dynamic Hamming distance (DHD).}

\usage{
seqdistmc(channels, method, norm="none", indel="auto", sm=NULL,
     with.missing=FALSE, full.matrix=TRUE, link="sum", cval=2,
     miss.cost=2, cweight=NULL, what="diss", ch.sep="@@@@TraMineRSep@@@@") }

\arguments{
 \item{channels}{A list of state sequence objects defined
    with the \code{\link{seqdef}} function, each state sequence object corresponding
    to a "channel".}

 \item{method}{a character string indicating the metric to be used.
    One of \code{"OM"} (Optimal Matching), \code{"LCS"} (Longest Common Subsequence),
    \code{"HAM"} (Hamming distance), \code{"DHD"} (Dynamic Hamming distance).}

  \item{norm}{
    String.
    Default: \code{"none"}.
    The normalization method to use. See \code{\link{seqdist}}.
    }

 \item{indel}{Double, vector of doubles, or list with an insertion/deletion cost or a vector of state dependent indel costs for each channel. Can also be \code{"auto"} (default), in which case the indel cost of each channel is automatically set in accordance with the \code{sm} value of the channel. See \code{indel} argument of \code{\link{seqdist}}.
     }

 \item{sm}{A list with a substitution-cost matrix for each channel
   or a list of method names for generating the substitution-costs
   (see \code{\link{seqcost}}). Ignored when \code{method="LCS"}.
   }

 \item{with.missing}{Logical or vector of logical. Must be \code{TRUE} for channels with non deleted
   gaps (missing values). See details.
   }

 \item{full.matrix}{Logical. If \code{TRUE} (default), the full distance matrix is returned.
   If \code{FALSE}, an object of class \code{\link{dist}} is returned.
   }

 \item{link}{Character string. One of \code{"sum"} or \code{"mean"}. Method to compute the "link"
   between channels. Default is to sum the substitution costs.
   }

 \item{cval}{Double. Substitution cost for \code{"CONSTANT"} matrix, see \code{\link{seqcost}}.}

 \item{miss.cost}{Double. Missing values substitution cost, see \code{\link{seqcost}}.}

 \item{cweight}{A vector of channel weights. Default is 1 (same weight for each channel).}

 \item{what}{Character string. What output should be returned? One of \code{"diss", "sm", "seqmc"}.}

 \item{ch.sep}{Character string. Separator used for building state names of the expanded alphabet.}
}

\details{The \code{seqdistmc} function builds a state sequence by combining the channels, derives the multichannel indel and substitution costs from the indel and substitution costs of each channel (following the strategy proposed by \cite{Pollock, 2007}), and computes the multichannel distances using the multichannel costs. The available metrics (see \code{method} argument) are optimal matching (\code{"OM"}), longest common subsequence (\code{"LCS"}), Hamming distance (\code{"HAM"}), and Dynamic Hamming Distance (\code{"DHD"}). For other edit distances, extract the combined state sequence object (by setting \code{what="seqmc"}) and the multichannel substitution and indel costs (by setting \code{what="sm"}). Then use these outcomes as input in a call to \code{seqdist}.  See \code{\link{seqdist}} for more information about available distances measures.

Normalization may be useful when dealing with sequences that are not all of the same length. For details on the applied normalization, see \code{\link{seqdist}}. }

\value{When \code{what="diss"}, a matrix of pairwise distances between multichannel sequences.\cr
When \code{what="sm"}, the matrix of substitution costs with three attributes: \code{indel} the indel cost(s), \code{alphabet} the alphabet of the combined state sequences, and \code{cweight} the channel weights used.\cr
When \code{what="seqmc"}, the combined state sequence object.
}

\seealso{ \code{\link{seqcost}}, \code{\link{seqdef}}, \code{\link{seqdist}}. }

\references{ Pollock, Gary (2007) Holistic trajectories: a study of combined employment, housing and family careers by using multiple-sequence analysis. \emph{Journal of the Royal Statistical Society: Series A} \bold{170}, Part 1, 167--183. }

\examples{data(biofam)

## Building one channel per type of event left, children or married
bf <- as.matrix(biofam[, 10:25])
children <-  bf==4 | bf==5 | bf==6
married <- bf == 2 | bf== 3 | bf==6
left <- bf==1 | bf==3 | bf==5 | bf==6

## Building sequence objects
child.seq <- seqdef(children)
marr.seq <- seqdef(married)
left.seq <- seqdef(left)

## Using transition rates to compute substitution costs on each channel
mcdist <- seqdistmc(channels=list(child.seq, marr.seq, left.seq),
 	method="OM", sm =list("INDELSLOG", "INDELSLOG", "TRATE"))

## Using a weight of 2 for children channel and specifying substitution-cost
smatrix <- list()
smatrix[[1]] <- seqsubm(child.seq, method="CONSTANT")
smatrix[[2]] <- seqsubm(marr.seq, method="CONSTANT")
smatrix[[3]] <- seqsubm(left.seq, method="TRATE")
mcdist2 <- seqdistmc(channels=list(child.seq, marr.seq, left.seq),
	method="OM", sm =smatrix, cweight=c(2,1,1))

## Retrieving the multichannel sequences
mcseq <- seqdistmc(channels=list(child.seq, marr.seq, left.seq),
	method="OM", sm =smatrix, cweight=c(2,1,1), what="seqmc", ch.sep="+")
alphabet(mcseq)

## Retrieving the multichannel substitution costs
mcsm <- seqdistmc(channels=list(child.seq, marr.seq, left.seq),
	method="OM", sm =smatrix, cweight=c(2,1,1), what="sm", ch.sep="+")
}


\author{Matthias Studer and Gilbert Ritschard}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Dissimilarity measures}
