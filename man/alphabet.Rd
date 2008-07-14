\name{alphabet}
\alias{alphabet}
\title{Retrieve the alphabet of a sequence object}
\description{
Retrieve the alphabet of a sequence object, that is, the list of all possible states (some states may not appear in the data).
}
\usage{
alphabet(seqdata)
}
\arguments{
  \item{seqdata}{a sequence object (see \code{\link{seqdef}} function).}
 }
\value{a vector containing the alphabet (states).} 
\seealso{
 \code{\link{seqdef}}
}
\examples{
## Creating a sequence object with the columns 13 to 24 
## in the 'actcal' example data set
data(actcal)
actcal.seq <- seqdef(actcal,13:24)

## Retrieving the alphabet
alphabet(actcal.seq)
}
\keyword{attribute}
