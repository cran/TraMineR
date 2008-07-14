\name{seqecontain}
\alias{seqecontain}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Event sequence contain event}
\description{
  Check if a sequence or a subsequence contain a specific event
}
\usage{
seqecontain(seq, eventList, exclude = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seq}{A list of event sequence}
  \item{eventList}{A list of event}
  \item{exclude}{if TRUE the search is exclusive, return false for if the subsequence contain an event that is not in eventList}
}
\details{
  Check every event sequence and check if contain an event in \code{eventList}.
  if \code{exclude} is TRUE the search is reversed, i.e. all event of the subsequences should be in \code{eventList}.
}
\value{
  A logical vector
}
\examples{
data(actcal.tse)
actcal.seqe <- seqecreate(actcal.tse$id,actcal.tse$time,actcal.tse$event)

#Searching for frequent subsequences, that is, appearing at least 20 times
fsubseq <- seqefsub(actcal.seqe,minSupport=20)

##looking for subsequence with FullTime
seqecontain(fsubseq$subseq,c("FullTime"))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{attribute}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
