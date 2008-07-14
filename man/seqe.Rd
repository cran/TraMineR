\name{is.seqe}
\alias{is.seqe}
\alias{is.seqelist}
\alias{print.seqe}
\alias{print.seqelist}
\alias{levels.seqelist}
\alias{[.seqelist}
\alias{levels.seqe}
\alias{Math.seqe}
\alias{Math.seqelist}
\alias{Ops.seqe}
\alias{Ops.seqelist}
\alias{Summary.seqe}
\alias{Summary.seqelist}


\title{Event sequences object}
\description{
  TraMineR use an internal sequence object to mine frequent event subsequences. is.seqe return true if the argument is an event sequence object.
  \code{is.seqelist} can be used to check for list of internal object.
  \code{as.seqelist} can be used to set a list as list of internal sequence object.

  There is a print method to handle this kind of object.
}
\usage{
is.seqe(s)
is.seqelist(s)
}
\arguments{
\item{s}{Sequence}
}
\seealso{\code{\link{str.seqelist}}, \code{\link{str.seqe}},\code{\link{as.character.seqelist}} and \code{\link{as.character.seqe}} to get character string of event object}
%\examples{}
\keyword{internal}
