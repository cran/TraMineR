\name{seqtab}
\alias{seqtab}
\title{Frequency table of the sequences}
\description{
Computes the frequency table of the sequences (count and percent of each sequence).
}
\usage{
 seqtab(seqdata, tlim=10, weighted=TRUE, format="SPS")
}
\arguments{
  \item{seqdata}{a sequence object as defined by the \code{\link{seqdef}} function.}
  \item{tlim}{if tlim>0, return frequencies only for the 'tlim' most frequent sequences. Default to 10.}
  \item{weighted}{if TRUE, frequencies account for the weights assigned to the state sequence object 
  (see \code{\link{seqdef}}). Set to FALSE if you want ignore the weights. If no weights were assigned during the creation of the sequence object, \code{weighted=TRUE} will yield the same result as \code{weighted=FALSE} since each sequence is allowed a weight of 1.}
  \item{format}{format used for displaying sequences as rownames in the output table. 
  Default is SPS format, which yields shorter and more readable sequence representations. 
  Alternatively, "STS" may be specified.}
}

\value{An object of class \code{stslist.freq}. This is actually a state sequence object (containing a list of state sequences) with added attributes, among others the \code{freq} attribute containing the frequency table. There are \code{print} and \code{plot} methods for such objects. More sophisticated plots can be produced with the \code{seqplot} function.}

\seealso{
 \code{\link{seqplot}, \link{plot.stslist.freq}}.
}
\examples{
## Creating a sequence object from the actcal data set
data(actcal)
actcal.lab <- c("> 37 hours", "19-36 hours", "1-18 hours", "no work")
actcal.seq <- seqdef(actcal, 13:24, labels=actcal.lab)

## 10 most frequent sequences in the data
seqtab(actcal.seq, tlim=10)

## With tlim=0, we get all distinct sequences in the data set 
## sorted according to their frequency
seqtab(actcal.seq, tlim=0)

## Example with weights
## from biofam data set using weigths
data(ex1)
ex1.seq <-  seqdef(ex1, 1:13, weights=ex1$weights)

## Unweighted frequencies
seqtab(ex1.seq, weighted=FALSE)

## Weighted frequencies
seqtab(ex1.seq, weighted=TRUE)
}
\keyword{univar}
