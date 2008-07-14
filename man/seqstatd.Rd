\name{seqstatd}
\alias{seqstatd}
\title{States frequency table and entropy}
\description{
  Returns the state frequencies, the number of valid states and the entropy of the states distribution for each time unit.
}
\usage{
 seqstatd(seqdata, digits=2, norm=TRUE)
}
\arguments{
  \item{seqdata}{a sequence object (see \code{\link{seqdef}} function).}
  \item{digits}{number of digits for the rounding of the results. Default to 2, set to NULL if you don't want any rounding.}
  \item{norm}{if TRUE (the default value), entropy is normalized, ie divided by the entropy of the alphabet. Set to FALSE if you want the entropy without normalization.}
}
\seealso{
 \code{\link{seqdplot}}.
}
\examples{
data(biofam)
biofam.seq <- seqdef(biofam,10:25)
sd <- seqstatd(biofam.seq)
barplot(sd$Entropy, main="Entropy of the states distribution, by age",
        ,xlab="Age",ylab="Entropy",col="green")
}
\keyword{univar}
