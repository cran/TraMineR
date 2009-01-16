\name{seqelength}
\alias{seqelength}
\alias{seqesetlength}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Length of event sequences}
\description{
The length of event sequence correspond to the total time of observation of the corresponding individual. This information is optional but can be used to perform survival analysis.
\code{seqelength} can be used to retrieve the length of some sequences.
\code{seqesetlength} can be used to set the length of these sequences.

}
\usage{
seqelength(s)
seqesetlength(s,len)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s}{A list of sequences.}
  \item{len}{A list of time corresponding to each sequences.}
}
\value{
  \code{seqelength} return a numeric vector with the length of each sequences.
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
