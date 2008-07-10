\name{seqistatd}
\alias{seqistatd}
\title{States frequency for each individual sequence}
\description{
  Returns the state frequencies for each sequence in the sequence object.
}
\usage{
 seqistatd(seqdata)
}
\arguments{
  \item{seqdata}{a sequence object (see \code{\link{seqdef}} function).}
}
%\seealso{}
\examples{
data(actcal)
actcal.seq <- seqdef(actcal,13:24)
seqistatd(actcal.seq[1:10,])
}
\keyword{univar}
