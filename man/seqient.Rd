\name{seqient}
\alias{seqient}
\title{Within sequences entropy}

\description{Within sequences entropy}

\details{
The seqient function returns the Shannon entropy of each sequence in \code{seqdata}. The entropy of a sequence is computed using the formula
\deqn{h(\pi_1,\ldots,\pi_s)=-\sum_{i=1}^{s}\pi_i\log \pi_i}{h(p_1,...,p_s) = - sum_{i=1}^{s} p_i log(p_i)}
%
where \eqn{s} is the size of the alphabet and \eqn{\pi_i}{p_i} the proportion of occurrences of the \eqn{i}th state in the considered sequence. The log is here the natural (base e) logarithm. The entropy can be interpreted as the `uncertainty' of predicting the states in a given sequence. If all states in the sequence are the same, the entropy is equal to 0. The maximum entropy for a sequence of length 12 with an alphabet of 4 states is 1.386294 and is attained when each of the four states appears 3 times. 
}  
\usage{
 seqient(seqdata, norm=TRUE, with.missing=FALSE)
}
\arguments{
  \item{seqdata}{a sequence object as returned by the the \code{\link{seqdef}} function.}

  \item{norm}{by default (TRUE), entropy is normalized, ie divided by the maximum entropy. The maximum entropy is computed as the entropy of the alphabet, ie an hypothetic sequence having all the states in the alphabet with equal length. Note that if for example the sequence length is uneven and the number of states in the alphabet is even, the theoretical maximum cannot be observed in the data.}

	\item{with.missing}{if set to TRUE, missing status (gaps in sequences) is handled as an additional state when computing the state distribution in the sequence.}
}

\value{
a vector whose number of elements is the number of sequences in \code{seqdata}, containing the entropy value of each sequence.
}

\references{
Gabadinho, A., G. Ritschard, M. Studer and N. S. Muller (2009). Mining Sequence Data in \code{R} with the \code{TraMineR} package: A user's guide. \emph{Department of Econometrics and Laboratory of Demography, University of Geneva}.
}

\seealso{
 \code{\link{seqstatd}} for the entropy of the transversal state distributions by positions in the sequence.
}
\examples{
data(actcal)
actcal.seq <- seqdef(actcal,13:24)

## Summarize and plots an histogram 
## of the within sequence entropy 
actcal.ient <- seqient(actcal.seq)
summary(actcal.ient)
hist(actcal.ient)

## Examples using with.missing argument
data(ex1)
ex1.seq <- seqdef(ex1, 1:13, weights=ex1$weights)

seqient(ex1.seq)
seqient(ex1.seq, with.missing=TRUE)
}
\keyword{univar}
