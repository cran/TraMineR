\name{seqient}
\alias{seqient}
\title{Within sequence entropies}

\description{Computes normalized or non-normalized within sequence entropies}

\details{
The seqient function returns the Shannon entropy of each sequence in \code{seqdata}. The entropy of a sequence is computed using the formula
%
\deqn{h(\pi_1,\ldots,\pi_s)=-\sum_{i=1}^{s}\pi_i\log \pi_i}{h(p_1,...,p_s) = - sum_{i=1}^{s} p_i log(p_i)}
%
where \eqn{s} is the size of the alphabet and \eqn{\pi_i}{p_i} the proportion of occurrences of the \eqn{i}th state in the considered sequence. The log is here the natural logarithm, i.e., the logarithm in base \eqn{e}. The entropy can be interpreted as the `uncertainty' of predicting the states in a given sequence. If all states in the sequence are the same, the entropy is equal to 0. The maximum entropy for a sequence of length 12 with an alphabet of 4 states is 1.386294 and is attained when each of the four states appears 3 times.

Normalization can be requested with the \code{norm=TRUE} option, in which case the returned value is the entropy divided by the entropy of the alphabet. The later is an upper bound for the entropy of sequences made from this alphabet. It exactly is the maximal entropy when the sequence length is a multiple of the alphabet size. The value of the normalized entropy is independent of the chosen logarithm base.
}
\usage{
 seqient(seqdata, norm=TRUE, base=exp(1), with.missing=FALSE)
}
\arguments{
  \item{seqdata}{a sequence object as returned by the the \code{\link{seqdef}} function.}

  \item{norm}{logical: should the entropy be normalized? \code{TRUE} by default. (see details)}
  \item{base}{real positive value: base of the logarithm used in the entropy formula (see details). If entropy is normalized (\code{norm=TRUE}), its value is the same whatever the base. Default is exp(1), i.e., the natural logarithm is used.}

  \item{with.missing}{logical: if \code{TRUE}, the missing state (gap in sequences) is handled as an additional state when computing the state distribution in the sequence.}
}

\value{
a vector with an entropy value for each sequence in \code{seqdata}; the vector length is equal to the number of sequences.
}

\references{
Gabadinho, A., G. Ritschard, N. S. Müller and M. Studer (2011). Analyzing and Visualizing State Sequences in R with TraMineR. \emph{Journal of Statistical Software} \bold{40}(4), 1-37.

Gabadinho, A., G. Ritschard, M. Studer and N. S. Müller (2009). Mining Sequence Data in \code{R} with the \code{TraMineR} package: A user's guide. \emph{Department of Econometrics and Laboratory of Demography, University of Geneva}.
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

\author{Alexis Gabadinho}

\keyword{Longitudinal characteristics}
\keyword{State sequences}
