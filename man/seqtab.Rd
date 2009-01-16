\name{seqtab}
\alias{seqtab}
\title{Sequences frequency table}
\description{
Returns a frequency table for the sequences (frequency and percent).
}
\usage{
 seqtab(seqdata, tlim=0)
}
\arguments{
  \item{seqdata}{a sequence object as defined by the \code{\link{seqdef}} function.}
  \item{tlim}{if tlim>0, shows only the 'tlim' most frequent sequences.}
}
\seealso{
 \code{\link{seqdplot}}.
}
\examples{
## Loading the 'actcal' example data set
data(actcal)

## Defining a sequence object with data in columns 13 to 24
## (activity status from january to december 2000)
actcal.seq <- seqdef(actcal,13:24,informat='STS')

## 10 most frequent sequences in the data
seqtab(actcal.seq,10)
}
\keyword{univar}
