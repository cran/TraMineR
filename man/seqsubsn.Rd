\name{seqsubsn}
\alias{seqsubsn}
\title{Number of distinct subsequences in a sequence.}
\description{
Computes the number of distinct subsequences in a sequence using Elzinga's algorithm.
}
\usage{
 seqsubsn(seqdata, DSS=TRUE)
}
\arguments{
  \item{seqdata}{a sequence object as defined by the \code{\link{seqdef}} function.}
\item{DSS}{if TRUE, the Distinct State Sequences (DSS, see \code{\link{seqdss}}) are first extracted, eg. the DSS contained in 'D-D-D-D-A-A-A-A-A-A-A-D' is 'D-A-D', and the number of distinct subsequences in the DSS is computed. If FALSE, the number of distinct subsequences is computed from sequences as they appear in the input sequence object. Hence the number of distinct subsequences is in most cases much higher with the DSS=FALSE option.}  
}
\value{a vector containing the number of distinct subsequences for each sequence in the input sequence object.}
\seealso{
 \code{\link{seqdss}}.
}
\examples{
data(actcal)
actcal.seq <- seqdef(actcal,13:24)

## Number of subsequences with DSS=TRUE 
seqsubsn(actcal.seq[1:10,])

## Number of subsequences with DSS=FALSE 
seqsubsn(actcal.seq[1:10,],DSS=FALSE)
}
\keyword{attribute}
