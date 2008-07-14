\name{seqdss}
\alias{seqdss}
\title{Extract distinct states sequence from a sequence object.}
\description{
Extract distinct states sequence from a sequence object. Returns a sequence object containing the distinct states sequences, ie the durations are not taken into account. The DSS contained in 'D-D-D-D-A-A-A-A-A-A-A-D' is 'D-A-D'. Durations can be extracted with the 'seqdur' function.
 }
\usage{
seqdss(seqdata)
}
\arguments{
   \item{seqdata}{a sequence object (see \code{\link{seqdef}} function).}
}
\value{a sequence object containing the distinct state sequence (DSS) for each sequence in the object given as argument.}
\seealso{
 \code{\link{seqdur}}.
}
\examples{
## Creating a sequence object with the columns 13 to 24 
## in the 'actcal' example data set
data(actcal)
actcal.seq <- seqdef(actcal,13:24)

## Retrieving the DSS
actcal.dss <- seqdss(actcal.seq)
  
## Displaying the DSS for the first 10 sequences
actcal.dss[1:10,]
}
\keyword{attribute}
