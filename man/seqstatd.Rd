\name{seqstatd}
\alias{seqstatd}
\title{Sequence of transversal state distributions and their entropies}
\description{
Returns the state frequencies, the number of valid states and the entropy of the state distribution at each position in the
sequence.
}
\usage{
 seqstatd(seqdata, weighted=TRUE, with.missing=FALSE, norm=TRUE)
}
\arguments{
  \item{seqdata}{a state sequence object as defined by the \code{\link{seqdef}} function.}

  \item{weighted}{if \code{TRUE}, distributions account for the weights assigned to the state sequence object
  (see \code{\link{seqdef}}).
  Set as \code{FALSE} if you want ignore the weights.}

  \item{with.missing}{If \code{FALSE} (default value), returned distributions ignore missing values.}

  \item{norm}{if \code{TRUE} (default value), entropy is normalized, ie divided by the entropy of the alphabet. Set as \code{FALSE} if
      you want the entropy without normalization.}
}

\details{
In addition to the state distribution at each position in the sequence, the \code{seqstatd} function provides also for each time
point the number of valid states and the Shannon entropy of the observed state distribution. Letting \eqn{p_i} denote the
proportion of cases in state \eqn{i} at the considered time point, the entropy is
\deqn{
h(p_1,\ldots,p_s) = -\sum_{i=1}^{s} p_i \log(p_i)
}
where \eqn{s} is the size of the alphabet. The log is here the natural (base e) logarithm. The entropy is 0 when all cases are in the same state and is maximal when the same proportion of cases are in each state. The entropy can be seen as a measure of the diversity of states observed at the considered time point. An application of such a measure (but with aggregated transversal data) can be seen in \cite{Billari
(2001)} and \cite{Fussell (2005)}.
}


\references{
Billari, F. C. (2001). The analysis of early life courses: complex descriptions of the transition to adulthood.
\emph{Journal of Population Research} 18 (2), 119-24.

Fussell, E. (2005). Measuring the early adult life course in Mexico: An application of the entropy index. In R. Macmillan
(Ed.), \emph{The Structure of the Life Course: Standardized? Individualized? Differentiated?}, Advances in Life Course
Research, Vol. 9, pp. 91-122. Amsterdam: Elsevier.
}

\seealso{
 \code{\link{plot.stslist.statd}} the plot method for objects of class \code{stslist.statd}, \cr
 \code{\link{seqdplot}} for higher level plot of transversal distributions and\cr
 \code{\link{seqHtplot}} for plotting the transversal entropy over sequence positions.
}
\examples{
data(biofam)
biofam.seq <- seqdef(biofam,10:25)
sd <- seqstatd(biofam.seq)
## Plotting the state distribution
plot(sd, type="d")

## Plotting the entropy indexes
plot(sd, type="Ht")

## ====================
## example with weights
## ====================
data(ex1)
ex1.seq <- seqdef(ex1, 1:13, weights=ex1$weights)

## Unweighted
seqstatd(ex1.seq, weighted=FALSE)

seqstatd(ex1.seq, weighted=TRUE)
}

\author{Alexis Gabadinho (with Gilbert Ritschard for the help page)}

\keyword{Transversal characteristics}
