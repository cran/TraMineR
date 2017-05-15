\name{seqST}
\alias{seqST}
\title{Sequences turbulence}
\description{
  Computes Elzinga's turbulence for each sequence in a sequence data set.
}

\details{
 Sequence turbulence is a measure proposed by \cite{Elzinga
& Liefbroer (2007)}. It is based on the number \eqn{\phi(x)}{phi(x)}
of distinct subsequences that can be extracted from the distinct
successive state sequence and the variance of the consecutive times \eqn{t_i}
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

The function searches for missing states in the sequences and if found, adds the missing state to the alphabet for the computation of the turbulence. In this case the \code{\link{seqdss}} and \code{\link{seqdur}} functions for extracting the distinct successive state sequences and the associated durations are called with the \code{{with.missing=TRUE}} argument. A missing state in a sequence is considered as the occurrence of an additional symbol of the alphabet, and two or more consecutive missing states are considered as two or more occurences of the same state. Hence the DSS of \code{A-A-*-*-*-B-B-C-C-D} is \code{A-*-B-C-D} and the associated durations are \code{2-3-2-2-1}.

The normalized value is obtained by subtracting 1 to the index and dividing by the turbulence value of a sequence made by repeating successively the alphabet up to the maximal length in \code{seqdata}.
}
\usage{
 seqST(seqdata, norm=FALSE)
}
\arguments{
  \item{seqdata}{a state sequence object as returned by the the \code{\link{seqdef}} function.}
  
  \item{norm}{logical: Should the turbulence index be normalized?}
}

\value{a vector of length equal to the number of sequences in
\code{seqdata} containing the turbulence value of each sequence. Normalized values are returned when \code{norm=TRUE}.}

\seealso{
 \code{\link{seqdss}}, \code{\link{seqdur}}. For another composite measure of sequence complexity see and \code{\link{seqici}}.
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
  actcal.seq <- seqdef(actcal[,13:24], informat='STS')

  ## Computing the sequences turbulence
  turb <- seqST(actcal.seq)

  ## Histogram for the turbulence
  hist(turb)
  
  ## Normalized turbulence
  turb.norm <- seqST(actcal.seq, norm=TRUE)

}

\author{Alexis Gabadinho and Gilbert Ritschard}

\keyword{Longitudinal characteristics}
