\name{seqST}
\alias{seqST}
\title{Sequences turbulence}
\description{
  Computes the turbulence for each sequence in a sequence data set, using the measure proposed by Elzinga.
}
\details{
	Returns a vector whose number of elements is the number of sequences in the input data set.
}
\usage{
 seqST(seqdata)
}
\arguments{
  \item{seqdata}{a sequence object (see \code{\link{seqdef}} function).}
}
\seealso{
 \code{}.
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
  actcal.seq <- seqdef(actcal,13:24,informat='STS')

  ## Computing the sequences turbulence
  turb <- seqST(actcal.seq)

  ## Histogram for the turbulence
  hist(turb)
}
\keyword{univar}
