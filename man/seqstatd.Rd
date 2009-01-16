\name{seqstatd}
\alias{seqstatd}
\title{Sequence of the states frequency tables and entropy of the states distributions}
\description{
Returns the state frequencies, the number of valid states and the entropy of the states distribution for each time unit.
}
\usage{
 seqstatd(seqdata, digits=2, norm=TRUE)
}
\arguments{
  \item{seqdata}{a sequence object as defined by the \code{\link{seqdef}} function.}

  \item{digits}{number of digits for the rounding of the results. Default to 2, set to NULL if you don't want any rounding.}

  \item{norm}{if TRUE (the default value), entropy is normalized, ie divided by the entropy of the alphabet. Set to FALSE if you want the entropy without normalization.}
}

\details{
In addition to the state distribution at each time point, the \code{seqstatd} function provides also for each time point the number of valid states and the Shannon entropy of the observed state distribution. Letting \eqn{p_i} denote the proportion of cases in state \eqn{i} at the considered time point, the entropy is
\deqn{
h(p_1,\ldots,p_s) = -\sum_{i=1}^s p_i \log_2(p_i)
}
where \eqn{s} is the size of the alphabet. The entropy is 0 when all cases are in the same state and is maximal when the same proportion of cases are in each state. The entropy can be seen as a measure of the diversity of states observed at the considered time point. An application of such a measure (but with aggregated transversal data) can be seen in \cite{Billari (2001)} and \cite{Fussell (2005)}. 
}

\references{
Billari, F. C. (2001). The analysis of early life courses: complex descriptions of the transition to adulthood. \emph{Journal of Population Research} 18 (2), 119-24.

Fussell, E. (2005). Measuring the early adult life course in Mexico: An application of the entropy index. In R. Macmillan (Ed.), \emph{The Structure of the Life Course: Standardized? Individualized? Differentiated?}, Advances in Life Course Research, Vol. 9, pp. 91-122. Amsterdam: Elsevier.
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
