\name{seqST}
\alias{seqST}
\title{Sequences turbulence}
\description{
  Computes the turbulence for each sequence in a sequence data set using the measure proposed by Elzinga.
}

\details{
 Sequence turbulence is a measure proposed by \cite{Elzinga
& Liefbroer (2007)}. It is based on the number \eqn{\phi(x)}{phi(x)}
of distinct subsequences that can be extracted from the distinct
state sequence and the variance of the consecutive times \eqn{t_i}
spent in the distinct states. For a sequence \eqn{x}, the formula is

  \deqn{T(x)=\log_{2}(\phi(x)\,\frac{s_{t,max}^2(x) + 1}{s_t^2(x) + 1})}{%
        T(x) = log_2 [phi(x) (s_{t,max}^2(x) + 1)/(s_t^2(x) + 1)]}

where \eqn{s_t^2(x)} is the variance of the successive state
durations in sequence \eqn{x} and \eqn{s_{t,max}^2(x)} is the maximum
value that this variance can take given the total duration of the
sequence. This maximum is computed as

\deqn{s_{t,max}^2 =(d-1)(1-\bar{t})^2}{%
      s_{t,max}^2 = (d-1)(1- bar{t})^2}

 where \eqn{\bar{t}}{bar{t}} is the mean consecutive time spent in the
distinct states, i.e. the sequence duration divided by the number
\eqn{d} of distinct states in the sequence.

}
\usage{
 seqST(seqdata)
}
\arguments{
  \item{seqdata}{a state sequence object as returned by the the \code{\link{seqdef}} function.}
}

\value{a vector of length equal to the number of sequences in
\code{seqdata} containing the turbulence value of each sequence. }

\seealso{
 \code{\link{seqient}} for computing the within sequence entropy.
}
\references{
  Elzinga, Cees H. and Liefbroer, Aart C. (2007).
  De-standardization of Family-Life Trajectories of Young Adults: A Cross-National Comparison Using Sequence Analysis.
  \emph{European Journal of Population}, 23, 225-250.
}
\examples{
  ## Loading the 'actcal' example data set
  data(actcal)

  ## Defining a sequence object with data in columns 13 to 24
  ## (activity status from january to december 2000)
  actcal.seq <- seqdef(actcal,13:24, informat='STS')

  ## Computing the sequences turbulence
  turb <- seqST(actcal.seq)

  ## Histogram for the turbulence
  hist(turb)
}
\keyword{univar}
