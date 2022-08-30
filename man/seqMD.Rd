\name{seqMD}
\alias{seqMD}
\alias{seqdistmc}
\title{Multidomain sequences}
\description{Build multidomain (MD) sequences of combined individual domain states (expanded alphabet), derive multidomain CAT indel and substitution costs from domain costs by means of an additive trick, and compute OM pairwise distances using CAT costs.
}

\usage{
seqMD(channels, method=NULL, norm="none", indel="auto", sm=NULL,
     with.missing=FALSE, full.matrix=TRUE, link="sum", cval=2,
     miss.cost=2, cweight=NULL, what="MDseq", ch.sep="+")

seqdistmc(channels, what="diss", ch.sep="@@@@TraMineRSep@@@@", ...)

}


\arguments{
 \item{channels}{A list of domain state sequence objects defined
    with the \code{\link{seqdef}} function, each state sequence object corresponding
    to a domain.}

 \item{method}{a character string indicating a dissimilarity measure between sequences.
    One of \code{"OM"} (Optimal Matching), \code{"LCS"} (Longest Common Subsequence),
    \code{"HAM"} (Hamming distance), \code{"DHD"} (Dynamic Hamming distance).}

  \item{norm}{
    String.
    Default: \code{"none"}.
    The normalization method to use. See \code{\link{seqdist}}.
    }

 \item{indel}{Double, vector of doubles, or list with an insertion/deletion cost or a vector of state dependent indel costs for each domain. Can also be \code{"auto"} (default), in which case the indel cost of each domain is automatically set in accordance with the \code{sm} value of the domain. See \code{indel} argument of \code{\link{seqdist}}.
     }

 \item{sm}{A list with a substitution-cost matrix for each domain
   or a list of method names for generating the domain substitution costs
   (see \code{\link{seqcost}}). Ignored when \code{method="LCS"}.
   }

 \item{with.missing}{Logical or vector of logical. Must be \code{TRUE} for channels with non deleted
   gaps (missing values). See details.
   }

 \item{full.matrix}{Logical. If \code{TRUE} (default), the full distance matrix between MD sequences is returned.
   If \code{FALSE}, an object of class \code{\link{dist}} is returned.
   }

 \item{link}{Character string. One of \code{"sum"} or \code{"mean"}. Method to compute the "link"
   between domains. Default is to sum substitution and indel costs.
   }

 \item{cval}{Double. Domain substitution cost for \code{"CONSTANT"} matrix, see \code{\link{seqcost}}.}

 \item{miss.cost}{Double. Cost to substitute missing values at domain level, see \code{\link{seqcost}}.}

 \item{cweight}{A vector of domain weights. Default is 1 (same weight for each domain).}

 \item{what}{Character string. What output should be returned? One of \code{"MDseq", "cost", "diss"}. The deprecated value \code{what="sm"} is treated as \code{what="cost"}.}

 \item{ch.sep}{Character string. Separator used for building state names of the expanded alphabet.}

 \item{...}{arguments passed to \code{seqMD}}
}

\details{The \code{seqMD} function builds MD sequences by combining the domain states. When \code{what="cost"}, it derives multidomain indel and substitution costs from the indel and substitution costs of each domain by means of the additive trick (AT) proposed by \cite{Pollock, 2007}). When \code{what="cost"}, it computes multidomain distances using the AT-multidomain costs. The available metrics (see \code{method} argument) are optimal matching (\code{"OM"}), longest common subsequence (\code{"LCS"}), Hamming distance (\code{"HAM"}), and Dynamic Hamming Distance (\code{"DHD"}). The selected metric is used to compute pairwise domain dissimilarities. It is also used to compute MD distances except when \code{"LCS"}, in which case MD distances are obtained with OM.   For other edit distances, extract the combined state sequence object (by setting \code{what="MDseq"}) and the AT-multidomain substitution and indel costs (by setting \code{what="cost"}). Then use these outcomes as input in a call to \code{seqdist}.  See \code{\link{seqdist}} for more information about available distance measures.

Normalization may be useful when dealing with sequences that are not all of the same length. For details on the applied normalization, see \code{\link{seqdist}}. }

\value{
When \code{what="MDseq"}, the MD sequences of combined states as a \code{stslist} sequence object.\cr
When \code{what="cost"}, the matrix of AT-substitution costs with three attributes: \code{indel} the AT-indel cost(s), \code{alphabet} the alphabet of the combined state sequences, and \code{cweight} the channel weights used.\cr
When \code{what="diss"}, a matrix of pairwise distances between MD sequences.
}

\seealso{ \code{\link{seqcost}}, \code{\link{seqdef}}, \code{\link{seqdist}}. }

\references{ Pollock, Gary (2007) Holistic trajectories: a study of combined employment, housing and family careers by using multiple-sequence analysis. \emph{Journal of the Royal Statistical Society: Series A} \bold{170}, Part 1, 167--183. }

\examples{
data(biofam)

## Building one channel per type of event left home, married, and child
cases <- 200
bf <- as.matrix(biofam[1:cases, 10:25])
left <- bf==1 | bf==3 | bf==5 | bf==6
married <- bf == 2 | bf== 3 | bf==6
children <-  bf==4 | bf==5 | bf==6

## Building sequence objects
left.seq <- seqdef(left)
marr.seq <- seqdef(married)
child.seq <- seqdef(children)
channels <- list(LeftHome=left.seq, Marr=marr.seq, Child=child.seq)

## AT-multidomain distances based on channel specific cost methods
MDdist <- seqMD(channels, method="OM",
    sm =list("INDELSLOG", "INDELSLOG", "TRATE"), what="diss")

## Providing channel specific substitution costs
smatrix <- list()
smatrix[[1]] <- seqsubm(left.seq, method="TRATE")
smatrix[[2]] <- seqsubm(marr.seq, method="CONSTANT")
smatrix[[3]] <- seqsubm(child.seq, method="CONSTANT")

## Retrieving the MD sequences
MDseq <- seqMD(channels)
alphabet(MDseq)

## Retrieving the AT-multidomain substitution costs
## Using a weight of 2 for domain "child"
MDcost <- seqMD(channels,
    sm=smatrix, cweight=c(1,1,2), what="cost")

## OMspell distances between MD sequences
MDdist2 <- seqdist(MDseq, method="OMspell",
    sm = MDcost, indel=attr(MDcost,"indel"))

}


\author{Gilbert Ritschard and Matthias Studer}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Dissimilarity measures}
