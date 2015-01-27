\name{seqlogp}
%\Rdversion{1.1}
\alias{seqlogp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Logarithm of the probabilities of state sequences
}
\description{
Compute the logarithm of the probability of each state sequence obtained from a state transition model. The probability of a sequence is equal to the product of each state probability of the sequence. There are several methods to compute a state probability.
}
\usage{
seqlogp(seqdata, prob="trate", time.varying=TRUE,
        begin="freq", weighted=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{seqdata}{The sequence to compute the probabilities.}
	
	\item{prob}{either the name (\code{"trate"} or \code{"freq"$} of the probability model to use to compute the state probabilities, or an \code{\link{array}} specifying the transition probabilities at each position \eqn{t} (see details).}
	
	\item{time.varying}{Logical. If \code{TRUE}, the probabilities (transitions or frequencies) are computed separately for each time \eqn{t} point.}
	
	\item{begin}{Model used to compute the probability of the first state. Either \code{"freq"} to use the observed frequencies on the first period or a vector specifying the probability of each state of the alphabet.}
	
	\item{weighted}{Logical. If \code{TRUE}, uses the weights specified in \code{seqdata} when computing the observed transition rates.}
}
\details{
The sequence likelihood \eqn{P(s)} is defined as the product of the probability with which each of its observed successive state is supposed to occur at its position.
Let \eqn{s=s_{1}s_{2} \cdots s_{\ell}}{s=s_1s_2 ... s_l} be a sequence of length \eqn{\ell}{l}. Then
\deqn{
	P(s)=P(s_{1},1) \cdot P(s_{2},2) \cdots P(s_{\ell},\ell)
}{
	P(s)=P(s_1, 1) * P(s_2, 2) * ... * P(s_l, l)
}
with \eqn{P(s_{t},t)}{P(s_t,t)} the probability to observe state \eqn{s_t} at position \eqn{t}.

The question is how to determinate the state probabilities \eqn{P(s_t,t)}. Several methods are available and can be set using the \code{prob} argument.

One commonly used method for computing them is to postulate a Markov model, which can be of various order. We can consider probabilities derived from the first order Markov model, that is, each \eqn{P(s_t,t)}, \eqn{t>1} is set as the transition rate \eqn{p(s_t|s_{t-1})}{p(s_t|s_(t-1))}. This is available in \code{seqlogp} by setting \code{prob="trate"}. \cr
The transition rates may be considered constant over time/positions (\code{time.varying=FALSE}), that is estimated across sequences from the observations at positions \eqn{t} and \eqn{t-1} for all \eqn{t} together. Time varying transition rates may also be considered (\code{time.varying=TRUE}), in which case they are computed separately for each position, that is estimated across sequences from the observations at positions \eqn{t} and \eqn{t-1} for each \eqn{t}, yielding an array of transition matrices. The user may also specify his own transition rates array or matrix.

Another method is to use the frequency of a state at each position to set \eqn{P(s_{t},t)}{P(s_t,t)} (\code{prob="freq"}). In the latter case, the probability of a sequence is independent of the probability of the transitions. Here again, the frequencies can be computed all together (\code{time.varying=FALSE}) or separately for each position \eqn{t} (\code{time.varying=TRUE}).
%
For \eqn{t=1}, we set \eqn{P(s_1,1)} to the observed frequency of the state \eqn{s_1} at position 1. Alternatively, the \code{begin} argument allows to specify the probability of the first state.

The likelihood \eqn{P(s)} being generally very small, \code{seqlogp} return \eqn{-\log P(s)}{-log(P(s))}. The latter quantity is minimal when \eqn{P(s)} is equal to \eqn{1}.
}
\value{
	A vector containing the logarithm of each sequence probability.
}

\examples{
## Creating the sequence objects using weigths
data(biofam)
biofam.seq <-  seqdef(biofam, 10:25, weights=biofam$wp00tbgs)

## Computing sequence probabilities
biofam.prob <- seqlogp(biofam.seq)
## Comparing the probability of each cohort
cohort <- biofam$birthyr>1940
boxplot(biofam.prob~cohort)

}

\author{Matthias Studer and Alexis Gabadinho (with Gilbert Ritschard for the help page)}

\keyword{State sequences}
\keyword{Longitudinal characteristics}  