\name{seqpm}
\alias{seqpm}
\title{Find patterns in sequences}
\description{
  Search for a pattern (subsequence) into sequences.
}
\usage{
 seqpm(seqdata, pattern)
}
\arguments{
  \item{seqdata}{a sequence object as defined by the \code{\link{seqdef}} function.}
  \item{pattern}{a character string representing the pattern (subsequence) to search for, without separator between the states.}
}
\details{
	This function searches a pattern (a character string) into a set of sequences and returns the results as a list with two elements: '\code{Nbmatch}' the number of occurrences of the pattern and '\code{MatchesIndex}' the vector of indexes (row numbers) of the sequences that match the pattern (see examples below).
	}

\value{a list with two elements (see details).}

%\seealso{}
\examples{
data(actcal)
actcal.seq <- seqdef(actcal,13:24)

## search for pattern "DAAD"
## (no work-full time work-full time work-no work)
## results are stored in the 'daad' object
daad <- seqpm(actcal.seq,"DAAD")

## Looking at the sequences
## containing the pattern
actcal.seq[daad$MIndex,]

## search for pattern "AD"
## (full time work-no work)
seqpm(actcal.seq,"AD")
}

\author{Alexis Gabadinho}
\keyword{State sequences}
