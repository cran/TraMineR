\encoding{latin1}

\name{seqsubm}

\alias{seqsubm}

\title{Create a substitution-cost matrix}

\description{
 The substitution-cost matrix is used when computing distances between
sequences by the method of optimal matching. The function creates the
substitution matrix using either a constant or the transition rates
computed from the sequence data or other methods to be implemented in
the future.
 }

\details{
 The substitution-cost matrix has dimension \eqn{ns*ns}, where
\eqn{ns} is the number of states in the \link{alphabet} of the
sequence object. The element \eqn{(i,j)} of the matrix is the cost of
substituting state \eqn{i} with state \eqn{j}.

 With the \code{"CONSTANT"} method, the substitution costs are the
same for all the states, with a default value of 2. An alternative
value can be provided by the user. When the \code{"TRATE"}
(transition rates) method is chosen, the transition rates between all
states are computed using the \link{seqtrate} function. The
substitution cost between states \eqn{i} and \eqn{j} is obtained with
the formula

 \deqn{SC(i,j) = cval -P(i,j) -P(j,i)}

where \eqn{P(i,j)} is the transition rate from state \eqn{i} to
\eqn{j}.
 }

\usage{
 seqsubm(seqdata, method, cval=NULL, with.missing=FALSE,
         miss.cost=NULL, time.varying=FALSE, weighted=TRUE,
		 transition="both", lag=1, missing.trate=FALSE)
}

 \arguments{
  \item{seqdata}{a sequence object as returned by the \link{seqdef} function.}

  \item{method}{method to compute transition rates. At this time, the methods available are constant
      value (\code{method="CONSTANT"}) or substitution costs using transition rates
      (\code{method="TRATE"})}
  \item{cval}{the constant substitution cost if method \code{"CONSTANT"} is chosen. For method \code{"TRATE"},
  the base value from which transition probabilities are subtracted. If \code{NULL}, \code{cval=2}, unless \code{transition}
  is set to \code{"both"} and \code{time.varying} is \code{TRUE} in which case \code{cval=4}.}
	
	\item{with.missing}{if \code{TRUE}, an additional entry is added in the matrix for the missing states.
Hence, a new "missing" state is added to the list of "valid" states. Use this if you want to compute
distances with missing values inside the sequences. See \cite{Gabadinho et al. (2010)} for more details on the
options for handling missing values when computing distances between sequences.}

	\item{miss.cost}{the substitution cost for the missing state. The default set it to \code{cval}}
	\item{time.varying}{Logical. If \code{TRUE} return an \code{\link{array}} containing a distinct
matrix for each time unit. The time is the third dimension (subscript).}
	\item{weighted}{Logical. If \code{TRUE} compute transition rates using weights specified in
\code{seqdata}.}
	\item{transition}{Only used if \code{time.varying=TRUE}. If \code{transition="both"}, it uses the transition
	rates from previous and next state. It can also be set to \code{"previous"} or \code{"next"}.
	}
	\item{lag}{Integer. Only used with (\code{method="TRATE"}). Time between the two states considered to compute
		transition rates (one by default).}
  \item{missing.trate}{Logical.  Only used with (\code{method="TRATE"}).
		If \code{TRUE}, substitution costs with missing state are also based on transition rates. If \code{FALSE} (default value),
		the substitution cost for the missing state are set to \code{miss.cost}.}
}
\seealso{
 \code{\link{seqtrate}}, \code{\link{seqdef}}, \code{\link{seqdist}}.
}

\references{
Gabadinho, A., G. Ritschard, N. S. Müller and M. Studer (2011). Analyzing and Visualizing State Sequences in R with TraMineR. \emph{Journal of Statistical Software} \bold{40}(4), 1-37.

Gabadinho, A., G. Ritschard, M. Studer and N. S. Müller (2010). Mining Sequence Data in
\code{R} with the \code{TraMineR} package: A user's guide. Department of Econometrics and
Laboratory of Demography, University of Geneva.
}

\examples{
## Defining a sequence object with columns 10 to 25
## in the 'biofam' example data set
data(biofam)
biofam.seq <- seqdef(biofam,10:25)

## Optimal matching using transition rates based substitution-cost matrix
## and insertion/deletion costs of 3
trcost <- seqsubm(biofam.seq, method="TRATE")
biofam.om <- seqdist(biofam.seq,method="OM",indel=3,sm=trcost)

## Optimal matching using constant value (2) substitution-cost matrix
## and insertion/deletion costs of 3
ccost <- seqsubm(biofam.seq, method="CONSTANT", cval=2)
biofam.om.c2 <- seqdist(biofam.seq, method="OM",indel=3,sm=ccost)

## Displaying the distance matrix for the first 10 sequences
biofam.om.c2[1:10,1:10]

## =================================
## Example with weights and missings
## =================================
data(ex1)
ex1.seq <- seqdef(ex1,1:13, weights=ex1$weights)

## Unweighted
subm <- seqsubm(ex1.seq, method="TRATE", with.missing=TRUE, weighted=FALSE)
ex1.om <- seqdist(ex1.seq, method="OM", sm=subm, with.missing=TRUE)

## Weighted
subm.w <- seqsubm(ex1.seq, method="TRATE", with.missing=TRUE, weighted=TRUE)
ex1.omw <- seqdist(ex1.seq, method="OM", sm=subm.w, with.missing=TRUE)

ex1.om == ex1.omw
}

\author{Matthias Studer and Alexis Gabadinho (first version) (with Gilbert Ritschard for the help page)}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Dissimilarity measures}
