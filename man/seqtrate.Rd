\name{seqtrate}
\alias{seqtrate}
\title{Compute transition rates between states}
\description{
Returns a matrix with transition rates between states, computed from a set of sequences.
}
\details{
Transition rates are the probabilities of transition from one state to another observed in the sequence data. Substitution costs based on transition rates can be used when computing distances between sequences with the optimal matching method (see \link{seqdist}).
}
\usage{
 seqtrate(seqdata, statl=NULL, time.varying=FALSE, weighted=TRUE,
		  lag=1, with.missing=FALSE)
}
\arguments{
  \item{seqdata}{a sequence object as defined by the \code{\link{seqdef}} function.}
  \item{statl}{a list of states or events for which the transition rates will be computed. If omitted (default), transition rates are computed between the distinct states in seqdata (obtained with the \code{\link{alphabet}} function).}
  \item{time.varying}{Logical. If \code{TRUE}, return an \code{\link{array}} containing a distinct matrix for each time unit. The time is the third dimension (subscript).}
  \item{weighted}{Logical. If \code{TRUE}, compute transition rates using weights specified in \code{seqdata}.}
  \item{lag}{Integer. Time between the two states considered to compute transition rates (one by default).}
  \item{with.missing}{Logical. If \code{FALSE} (default value), returned transition rates ignore missing values.}
}

\value{a matrix of dimension \eqn{ns*ns}, where \eqn{ns} is the number of states in the \link{alphabet} of the sequence object.}

\references{
Gabadinho, A., G. Ritschard, N. S. Müller and M. Studer (2011). Analyzing and Visualizing State Sequences in R with TraMineR. \emph{Journal of Statistical Software} \bold{40}(4), 1-37.
}

\seealso{
 \code{\link{seqdist}} \code{\link{seqsubm}} \code{\link{alphabet}}.
}
\examples{
  ## Loading the 'actcal' example data set
  data(actcal)

  ## Defining a sequence object with data in columns 13 to 24
  ## (activity status from january to december 2000)
  actcal.seq <- seqdef(actcal,13:24,informat='STS')

  ## Computing transition rates
  seqtrate(actcal.seq)

  ## Computing transition rates between states "A" and "B" only
  seqtrate(actcal.seq, c("A","B"))

## ====================
## Example with weights
## ====================
data(ex1)
ex1.seq <- seqdef(ex1,1:13, weights=ex1$weights)

seqtrate(ex1.seq, weighted=FALSE)
seqtrate(ex1.seq, weighted=TRUE)
}
\keyword{nonparametric}
