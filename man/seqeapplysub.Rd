\encoding{latin1}
\name{seqeapplysub}
\alias{seqeapplysub}
\title{Checking if event sequences contain given subsequences}
\description{
Checks occurrences of the subsequences \code{subseq} among the event sequences and returns the result according to the selected \code{method}.
}
\usage{
seqeapplysub(subseq, method = "count", constraint = NULL, rules=FALSE)
}
\arguments{
 \item{subseq}{list of subsequences (an event subsequence object) such as created by \code{\link{seqefsub}}}
 \item{method}{type of result, should be one of "count", "presence" or "age"}
 \item{constraint}{
	Time constraints overriding those used to compute subseq. See \code{\link{seqeconstraint}}
	}
 \item{rules}{If set to TRUE, instead of checking occurences of the subsequences among the event sequences, check the occurence of the subsequences inside the subsequences (internally used by \code{seqerules})}
}
\details{
  There are three methods implemented:
  'count' counts the number of occurrence of each given subsequence in each event sequence;
  'presence' returns 1 if the subsequence is present, 0 otherwise;
  'age' returns the age of appearance of each subsequence in each event sequence. 
In case of multiple possibilities, the age of the first occurrence is returned. When the subsequence is not in the sequence, -1 is returned.
}

\value{
The return value is a matrix where each row corresponds to a sequence (row names are set accordingly) and each column corresponds to a subsequence (col names are set accordingly).
The cells of the matrix contain the requested values (count, presence-absence indicator or age).
}

\seealso{\code{\link{seqecreate}} for more information on event sequence object and \cite{Gabadinho et al. (2009)} on how to use the event sequence analysis module.}

\references{
Gabadinho, A., G. Ritschard, M. Studer and N. S. Müller (2009). Mining Sequence Data in \code{R} with \code{TraMineR}: A user's guide for version 1.1. \emph{Department of Econometrics and Laboratory of Demography, University of Geneva}.
}


\examples{
## Loading data
data(actcal.tse)

## Creating the event sequence object
actcal.seqe <- seqecreate(actcal.tse)

## Printing sequences
actcal.seqe[1:10]

## Looking for frequent subsequences
fsubseq <- seqefsub(actcal.seqe,pMinSupport=0.01)

## Counting the number of occurrences of each subsequence
msubcount <- seqeapplysub(fsubseq,method="count")
## First lines...
msubcount[1:10,1:10]
## Presence-absence of each subsequence
msubpres <- seqeapplysub(fsubseq,method="presence")
## First lines...
msubpres[1:10,1:10]

## Age at first appearance of each subsequence
msubage <- seqeapplysub(fsubseq,method="age")

## First lines...
msubage[1:10,1:10]
}
\keyword{nonparametric}

