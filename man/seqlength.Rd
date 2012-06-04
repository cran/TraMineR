\name{seqlength}
\alias{seqlength}
\title{Sequence length}
\description{
  Returns the length of sequences.
}
\details{
The length of a sequence is computed by eliminating the missing values at the end (right) and counting the number of states or events. The seqlength function returns a vector containing the length of each sequence in the sequence object given as argument.
}
\usage{
 seqlength(seqdata)
}
\arguments{
  \item{seqdata}{a sequence object created with the \code{\link{seqdef}} function.}
}
%\seealso{}
\examples{
## Loading the 'famform' example data set
data(famform)

## Defining a sequence object with the 'famform' data set
ff.seq <- seqdef(famform)
 
## Retrieving the length of the first 10 sequences 
## in the ff.seq sequence object
seqlength(ff.seq)
}
\keyword{attribute}
