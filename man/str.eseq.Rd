\name{str.eseq}
\alias{str.eseq}
\alias{str.seqelist}
\alias{as.character.eseq}
\alias{as.character.seqelist}

\title{String representation of an event sequence object}
\description{
  Return a string representation of an event sequence object. Works on event sequences as well as on subsequences.
  This is an overriding of the usual corresponding functions (\code{str} and \code{as.character})
}
\details{
   The string have the following form:
  time-(e1,e2,...)-elapsedtime-(e2,...)...
  Where time is the time elapsed from previous set of events and (e1,e2,...) a non empty list of simultaneous events.
  The function also works on list.
  Use \code{str.seqelist} to get a single string for a whole list, and \code{as.character.seqelist} to get a vector of strings.
  \code{str.eseq} and \code{as.character.eseq} have mainly an internal purpose (they are called by print methods for example).
}

\seealso{\code{\link{seqecreate}} for a full example }
\examples{
data(actcal.tse)
actcal.eseq <- seqecreate(actcal.tse)
head(as.character(actcal.eseq))
}

\author{Matthias Studer (with Gilbert Ritschard for the help page)}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

\keyword{internal}
\keyword{Event Sequences}
