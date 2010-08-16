\name{seqdss}
\alias{seqdss}
\title{Extract distinct states sequence from a sequence object}
\description{
Extract distinct states sequence from a sequence object. 
}
\details{
Returns a sequence object containing the distinct states sequences, ie the durations are not taken into account. The DSS contained in 'D-D-D-D-A-A-A-A-A-A-A-D' is 'D-A-D'. Associated durations can be extracted with the \code{\link{seqdur}} function.

If called with the \code{{with.missing=TRUE}} argument, a missing state in a sequence is considered as the occurence of an additional symbol of the alphabet, and two or more consecutive missing states are considered as two or more occurences of the same state. Hence the DSS of \code{A-A-*-*-*-B-B-C-C-D} is \code{A-*-B-C-D}.
}
\usage{
seqdss(seqdata, with.missing=FALSE)
}
\arguments{
   \item{seqdata}{a sequence object as defined by the \code{\link{seqdef}} function.}
	\item{with.missing}{if set to TRUE, missing statuses (gaps in sequences) also appear in the DSS. See \link{seqdef} on options for handling missing values when creating sequence objects.}
}
\value{a sequence object containing the distinct state sequence (DSS) for each sequence in the object given as argument. }
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

## Example with with.missing argument
data(ex1)
ex1.seq <- seqdef(ex1, 1:13)

seqdss(ex1.seq)
seqdss(ex1.seq, with.missing=TRUE)
}
\keyword{attribute}
