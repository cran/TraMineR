\name{seqient}
\alias{seqient}
\title{Sequence entropy}
\description{
  Returns the within sequences entropy.
}
\details{
Computes the entropy within each sequence (row) in the dataframe or matrix. The entropy of a sequence equals 0 if all states in the sequence are the same. The function returns a vector with the same number of rows (= the number of sequences) as it's argument. Other measures of entropy are available: 'seqstatd' returns the entropy of the distribution of states for each time unit.
}  
\usage{
 seqient(seqdata, norm=TRUE)
}
\arguments{
  \item{seqdata}{a sequence object (see \code{\link{seqdef}} function).}
  \item{norm}{by default (TRUE), entropy is normalized, ie divided by the maximum entropy. The maximum entropy is computed as the entropy of the alphabet, ie an hypothetic sequence having all the states in the alphabet with equal length. Note that if for example the sequence length is uneven and the number of states in the alphabet is even, the theoretical maximum cannot be observed in the data.}
}
\seealso{
 \code{\link{seqstatd}}.
}
\examples{
data(actcal)
actcal.seq <- seqdef(actcal,13:24)

## Summarize and plots an histogram 
## of the within sequence entropy 
actcal.ient <- seqient(actcal.seq)
summary(actcal.ient)
hist(actcal.ient)
}
\keyword{univar}
