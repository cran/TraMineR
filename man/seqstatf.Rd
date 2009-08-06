\name{seqstatf}
\alias{seqstatf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{State frequencies in the all sequence data set}
\description{
Frequency of each state of the alphabet in the all sequence data set.
}
\usage{
seqstatf(seqdata, weighted = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seqdata}{a sequence object as defined by the \code{\link{seqdef}} function.}
  \item{weighted}{if TRUE, frequencies account for the weights assigned to the state sequence object 
  (see \code{\link{seqdef}}). Set to FALSE if you want ignore the weights. If no weights were assigned during the creation of the sequence object, \code{weighted=TRUE} will yield the same result as \code{weighted=FALSE} since each sequence is allowed a weight of 1.}
}
\details{
The \code{seqstatf} function computes the (weighted) raw and percentage frequency of each state of the alphabet in \code{seqdata}, i.e the (weighted) sum of the occurences of a state in \code{seqdata}.  
}
\value{
a data.frame with as many rows as the number of states in the alphabet and two columns, one for the raw frequencies (Freq) and one for the percentage frequencies.
}

\seealso{
\code{\link{seqstatd}} for the state distribution by time point (position), \code{\link{seqistatd}} for the state distribution within each sequence.
}

\examples{
## Creating a sequence object from the actcal data set
data(actcal)
actcal.lab <- c("> 37 hours", "19-36 hours", "1-18 hours", "no work")
actcal.seq <- seqdef(actcal, 13:24, labels=actcal.lab)

## States frequencies
seqstatf(actcal.seq)

## Example with weights
data(ex1)
ex1.seq <- seqdef(ex1, 1:13, weights=ex1$weights)

## Unweighted
seqstatf(ex1.seq, weighted=FALSE)

## Weighted
seqstatf(ex1.seq, weighted=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{univar}

