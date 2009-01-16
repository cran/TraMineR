\name{seqeapplysub}
\alias{seqeapplysub}
\title{Applying Subsequences to Event Sequences}
\description{
Apply every subsequences \code{subseq} to each event sequences and compute the result of the given \code{method}.
}
\usage{
seqeapplysub(subseq, method = "count", constraint = NULL)
}
\arguments{
 \item{subseq}{a list of subsequences}
 \item{method}{type of the result expected, should be one of "count", "presence" or "age"}
 \item{constraint}{
	Other time constraint overriding those used to compute subseq. See \code{\link{seqeconstraint}}
	}
}
\details{
  They are three methods implemented:
  'count' counts the number of occurrence of a given subsequences in an event sequence;
  'presence' returns one if the subsequence is present, 0 otherwise;
  'age' returns the age of appearance of a subsequence in an event sequence. 
In case of multiple possiblities, the age of the first occurences is returned. When the subsequences is not in the sequence, -1 is returned.
}

\value{
The return value is a matrix where every row correspond to a sequence (row names are set accordingly) and each column correspond to a subsequence (col names are set accordingly).
The matrix store the results of the count method.
}

\seealso{See Also as \code{\link{seqecreate}} for more information on how to use event sequence analysis module.}
\examples{
## Loading data
data(actcal.tse)

## Creating sequences
actcal.seqe <- seqecreate(actcal.tse)

## Printing sequences
actcal.seqe[1:10]

## Looking for frequent subsequences
fsubseq <- seqefsub(actcal.seqe,pMinSupport=0.01)

## Counting the number of occurence of each subsequence
msubcount <- seqeapplysub(fsubseq,method="count")
## First lines...
msubcount[1:10,1:10]
## Presence-absence of each subsequence
msubpres <- seqeapplysub(fsubseq,method="presence")
## First lines...
msubpres[1:10,1:10]

## Age of first appearance of each subsequence
msubage <- seqeapplysub(fsubseq,method="age")

## First lines...
msubage[1:10,1:10]
}
\keyword{nonparametric}

