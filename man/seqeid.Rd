\name{seqeid}
\alias{seqeid}
\title{Retrieve id of an event sequence object.}
\description{
  Retrieve id of an event sequence or a list of event sequence object.
}
\usage{
seqeid(s)
}
\arguments{
\item{s}{A sequence or a list of sequence}
}
%\seealso{}
\examples{
data(actcal.tse)
actcal.seqe<-seqecreate(actcal.tse$id,actcal.tse$time,actcal.tse$event)
seqeid(actcal.seqe)
}
\keyword{attribute}
