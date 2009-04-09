\name{seqelength}
\alias{seqelength}
\alias{seqesetlength}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Length of event sequences}
\description{
The length of an event sequence is its time span, i.e. the total time of observation. This information is optional but may be useful to perform for instance a survival analysis.
\code{seqelength} retrieves the length the given sequences.
\code{seqesetlength} sets the length of the sequences.

}
\usage{
seqelength(s)
seqesetlength(s,len)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s}{An event sequence object (\code{seqelist}).}
  \item{len}{A list of sequence lengths.}
}
\value{
  \code{seqelength} returns a numeric vector with the length of each sequences.
}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~}

%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
data(actcal.tse)
actcal.seqe <- seqecreate(actcal.tse)
##time to end is added
sl <- numeric()
sl[1:2000] <- 12
##All sequences with same length
seqesetlength(actcal.seqe, sl)
actcal.seqe[1:10]
##Retrieve length
seqelength(actcal.seqe)
}
\keyword{attribute}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
