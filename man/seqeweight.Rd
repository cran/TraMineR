\name{seqeweight}
\alias{seqeweight}
\alias{seqeweight<-}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Setting or retrieving weights of an event sequence object.
}
\description{
	Event sequence objects can be weighted. Weights may be used by other functions such as \code{\link{seqefsub}} or \code{\link{seqecmpgroup}} to compute weighted statistics.
}
\usage{
seqeweight(s)
seqeweight(s) <- value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s}{A list of event sequences}
  \item{value}{Numerical vector containing weights}
}
\value{
\code{seqeweight} returns a numerical vector containing the weights associated to each event sequence.
}

\examples{
##Starting with states sequences
##Loading data
data(biofam)
## Creating state sequences
biofam.seq <- seqdef(biofam,10:25,informat='STS')

## Creating event sequences from biofam
biofam.seqe <- seqecreate(biofam.seq, weighted=FALSE)

## Using the weights
seqeweight(biofam.seqe) <- biofam$wp00tbgs

## Now seqefsub accoounts for weights unless weighted is set to FALSE
fsubseq <- seqefsub(biofam.seqe, pMinSupport=0.01)

## Searching for weighted susbsequences which best 
## discriminate the birth cohort
discr <- seqecmpgroup(fsubseq, group=biofam$birthyr>=1940)
plot(discr[1:15])

}
