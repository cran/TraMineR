\name{seqtrate}
\alias{seqtrate}
\title{Transition rates between states.}
\description{
Returns a matrix with transition rates between states, computed from a set of sequences.
}
\details{
Transition rates are the probabilities of transition from one state to another in the sequence data. Substitution costs based on transition rates can be used when computing distances between sequences with the optimal matching method.
}
\usage{
 seqtrate(seqdata, statl=NULL)
}
\arguments{
  \item{seqdata}{a sequence object (see \code{\link{seqdef}} function).}
  \item{statl}{a list of states or events for which the transition rates will be computed. If omitted (default), transition rates are computed between the distinct states in seqdata (obtained with the \code{\link{alphabet}} function).}
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
  seqtrate(actcal.seq,c("A","B"))	
}
\keyword{nonparametric}
