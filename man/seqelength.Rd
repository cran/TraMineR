\name{seqelength}
\alias{seqelength}
\alias{seqelength<-}
\alias{seqesetlength}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Lengths of event sequences}
\description{
The length of an event sequence is its time span, i.e., the total time of observation. This information is useful to perform for instance a survival analysis. The function
\code{seqelength} retrieves the lengths of the given sequences, while
\code{seqelength <-} sets the length of the sequences.

\code{seqesetlength} is deprecated.

}
\usage{
seqelength(s)
seqelength(s) <- value
seqesetlength(s, len)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s}{An event sequence object (\code{seqelist}).}
  \item{len}{A list of sequence lengths.}
  \item{value}{A list of sequence lengths.}
}
\value{
  \code{seqelength} returns a numeric vector with the lengths of the sequences.
}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~}

%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
data(actcal.tse)
actcal.seqe <- seqecreate(actcal.tse)
## Since endEvent is not specified, contains no sequence lengths
## We set them manually as 12 for all sequences
sl <- numeric()
sl[1:2000] <- 12
seqelength(actcal.seqe) <- sl
actcal.seqe[1:10]
## Retrieve lengths
seqelength(actcal.seqe)
}
\keyword{attribute}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
