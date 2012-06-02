\name{str.seqe}
\alias{str.seqe}
\alias{str.seqelist}
\alias{as.character.seqe}
\alias{as.character.seqelist}

\title{String representation of sequence object}
\description{
  Return a string representation of an event sequence object. This work on sequence or subsequences as well.
  This is an overiding of the usual corresponding function (\code{str} and \code{as.character})
}
\details{
   The string have the following form:
  time-(e1,e2,...)-elapsedtime-(e2,...)...
  Where time is the time elapsed from previous set of events and (e1,e2,...) a non empty list of simultaneous events.
  This function also work on list.
  You will mainly used \code{str.seqelist} wich a return a single string for a whole list, and \code{as.character.seqelist} which return a vector of string.
  \code{str.seqe} and \code{as.character.seqe} have mainly an internal purpose (there are called by print and more).
}

\seealso{\code{\link{seqecreate}} for a full example }
\examples{
data(actcal.tse)
actcal.seqe <- seqecreate(actcal.tse)
as.character(actcal.seqe)
}
\keyword{internal}
