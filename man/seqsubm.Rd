\name{seqsubm}
\alias{seqsubm}
\title{Create a substitution-cost matrix}
\description{
The substitution-cost matrix is used when computing distances between sequences by the method of optimal matching. The function creates the substitution matrix using either a constant or the transition rates computed from the sequence data or other methods to be implemented in the future.
}
\details{
	The substitution-cost matrix has dimension ns*ns, where ns is the number of distinct states found in the sequence dataset. The element (i,j) of the matrix is the cost of substituting state i whith state j.
	In the "constant" method, the substitution costs are the same for all the states, the value is provided by the user. In the
	"transition rates" method, the transition rates between all states are computed using the seqtrate function. The substitution cost between states 'Si' and 'Sj' is obtained with the formula SC(i,j) = 2 -P(Si,Sj) -P(Sj,Si) where P(Si,Sj) is the transition rate between states i and j.  
}
\usage{
 seqsubm(seqdata, method, cval)
}
\arguments{
  \item{seqdata}{a sequence object created with the \code{\link{seqdef}} function.}
  \item{method}{method to compute transition rates. At this time, the methods available are constant value (method="CONSTANT") or substitution costs using transition rates (method="TRATE")}
  \item{cval}{the constant substitution cost if method "CONSTANT" is choosen. Otherwise, do not specify.}
}
\seealso{
 \code{\link{seqtrate}} \code{\link{seqdef}}.
}
\examples{
  ## Defining a sequence object with columns 10 to 25 
  ## in the 'biofam' example data set 
  data(biofam)
  biofam.seq <- seqdef(biofam,10:25)

  ## Optimal matching using transition rates based substitution-cost matrix
  ## and insertion/deletion costs of 3
  trcost <- seqsubm(biofam.seq, method="TRATE")
  biofam.om <- seqdist(biofam.seq,method="OM",indel=3,sm=trcost)

  ## Optimal matching using constant value (2) substitution-cost matrix
  ## and insertion/deletion costs of 3
  ccost <- seqsubm(biofam.seq, method="CONSTANT", cval=2)
  biofam.om.c2 <- seqdist(biofam.seq, method="OM",indel=3,sm=ccost)

  ## Displaying the distance matrix for the first 10 sequences
  biofam.om.c2[1:10,1:10]
}
\keyword{misc}
