\name{seqecontain}
\alias{seqecontain}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Check if sequence contains events}
\description{
  Check if a sequence or a subsequence contains given events
}
\usage{
seqecontain(seq, eventList, exclude = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seq}{A event sequence object (\code{seqelist}) or a an event subsequence object (\code{subseqelist})}
  \item{eventList}{A list of events}
  \item{exclude}{if \code{TRUE} the search is exclusive and returns \code{FALSE} for any subsequence containing an event that is not in \code{eventList}}
}
\details{
  Checks, for each provided event sequence, if it contains one of the events in \code{eventList}.
  If \code{exclude} is \code{TRUE}, \code{seqecontain} looks if all events of the subsequence are in \code{eventList}.
}
\value{
  A logical vector.
}
\seealso{
\code{\link{seqecreate}} for creating event sequence objects and \code{\link{seqefsub}} for creating event subsequence objects.
}
\examples{
data(actcal.tse)
actcal.seqe <- seqecreate(actcal.tse)

##Searching for frequent subsequences, that is appearing at least 20 times
fsubseq <- seqefsub(actcal.seqe,minSupport=20)

##looking for subsequence with FullTime
seqecontain(fsubseq,c("FullTime"))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{attribute}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
