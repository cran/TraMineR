\name{seqistatd}
\alias{seqistatd}
\title{States frequency for each individual sequence}
\description{
  Returns the state frequencies for each sequence in the sequence object.
}
\usage{
 seqistatd(seqdata, with.miss=FALSE)
}
\arguments{
  \item{seqdata}{a sequence object (see \code{\link{seqdef}} function).}
  \item{with.miss}{if set to TRUE, cumulated durations are also computed for the missing status (gaps in the sequences). See \link{seqdef} on options for handling missing values when creating sequence objects.} 
}
%\seealso{}
\examples{
data(actcal)
actcal.seq <- seqdef(actcal,13:24)
seqistatd(actcal.seq[1:10,])
}
\keyword{univar}
