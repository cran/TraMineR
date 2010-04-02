\name{seqnum}
\alias{seqnum}
\title{Translate a sequence object's alphabet into numerical alphabet, ranging 0-(nbstates-1).}
\description{
  If the alphabet (the list of possible states or events in a set of sequences) is composed of characters, this function converts the sequence data using a numerical alphabet. The first state (for exemple 'A') is coded with the value '0', the second state (for exemple 'B') is coded with the value '1', etc... The function returns a sequence object containing the original sequences coded with the new numerical alphabet. 
}
\usage{
 seqnum(seqdata, with.missing=FALSE)
}
\arguments{
	\item{seqdata}{a sequence object as defined by the \code{\link{seqdef}} function.}

	\item{with.missing}{if TRUE, missing elements in the sequences are turned into numerical values as well. The code for missing values in the sequences is retrieved as the 'nr' attribute of \code{seqdata}.}
}
%\seealso{}

\examples{
data(actcal)
actcal.seq <- seqdef(actcal,13:24)

## The first 10 sequences in the actcal.seq
## sequence object
actcal.seq[1:10,]
alphabet(actcal.seq)

## The first 10 sequences in the actcal.seq
## sequence object with numerical alphabet
seqnum(actcal.seq[1:10,])

## states A,B,C,D are now coded 0,1,2,3
alphabet(seqnum(actcal.seq))
}
\keyword{manip}
