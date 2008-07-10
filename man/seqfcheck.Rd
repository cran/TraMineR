\name{seqfcheck}
\alias{seqfcheck}
\title{Check if sequences are in compressed format.}
\description{
  Wether the sequence(s) are in compressed format is checked by counting the number of columns and searching for the '-' or ':' separator. The function returns the separator if it has been found in the data. If the data contains more than one column, the data is supposed to be in the extended format, and 'X' is returned, unless some states codes contain the '-' character (eg. states coded with negative integer values), in which case '-X' is returned.
}
\usage{
 seqfcheck(seqdata)
}
\arguments{
 \item{seqdata}{a vector, data frame or matrix containing sequence data.}
 }
\seealso{
 \code{\link{seqconc} \link{seqdecomp}}.
}
\examples{
data(actcal)
seqfcheck(actcal[,13:24])
}
\keyword{internal}
